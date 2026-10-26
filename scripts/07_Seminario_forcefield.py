###################################################################################################################
#   Automated Force Fields for Metals / easyPARM
#   EP2 provenance-enabled variant of 07_Seminario_forcefield.py
#
#   Default EP2 policy: keep Hessian/Seminario BOND and ANGLE terms unaveraged.
#   Rationale: averaging by atom-type key can erase chemically distinct local terms.
#
#   Optional environment control:
#     EP2_HESSIAN_BOND_ANGLE_POLICY=all_unaveraged       # default
#     EP2_HESSIAN_BOND_ANGLE_POLICY=selected_unaveraged  # preserve S/s, ce, os, metal terms; average others
#     EP2_HESSIAN_BOND_ANGLE_POLICY=average_all          # original easyPARM-style averaging, provenance-aware
###################################################################################################################

from collections import defaultdict
from pathlib import Path
import os
import re

PROVENANCE_FILE = "ep2_parameter_provenance.tsv"
POLICY = os.environ.get("EP2_HESSIAN_BOND_ANGLE_POLICY", "all_unaveraged").strip().lower()
VALID_POLICIES = {"all_unaveraged", "selected_unaveraged", "average_all"}
if POLICY not in VALID_POLICIES:
    raise RuntimeError(
        f"Invalid EP2_HESSIAN_BOND_ANGLE_POLICY={POLICY!r}. "
        f"Valid: {', '.join(sorted(VALID_POLICIES))}"
    )


SECTION_HEADERS = {"MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"}


def append_provenance(section, key, atom_ids, value_1, value_2, source, source_file, notes=""):
    path = Path(PROVENANCE_FILE)
    write_header = not path.exists()
    with path.open("a") as handle:
        if write_header:
            handle.write("section\tatom_type_key\tatom_ids\tvalue_1\tvalue_2\tsource\tsource_file\tnotes\n")
        handle.write(
            f"{section}\t{key}\t{atom_ids}\t{value_1}\t{value_2}\t"
            f"{source}\t{source_file}\t{notes}\n"
        )


def split_comment(line):
    if "#" in line:
        before, comment = line.split("#", 1)
        return before.rstrip(), comment.strip()
    return line.rstrip(), ""


def extract_atoms_from_comment(comment):
    match = re.search(r"(?:^|\s)atoms=([^\s#]+)", comment)
    return match.group(1) if match else ""


def parse_parameter_line(line):
    """Return (key, value1, value2, comment) or None for empty/unparseable lines."""
    data_part, comment = split_comment(line)
    if not data_part.strip():
        return None
    fields = data_part.split()
    if len(fields) < 3:
        return None
    try:
        value1 = float(fields[1])
        value2 = float(fields[2])
    except ValueError:
        return None
    return fields[0], value1, value2, comment


def read_metal_atom_types():
    """Read metal atom types from metal_number.dat and NEW_COMPLEX.mol2 when available."""
    metal_ids = []
    try:
        with open("metal_number.dat", "r") as handle:
            metal_ids = [int(x.strip()) for x in handle if x.strip()]
    except FileNotFoundError:
        return set()
    except ValueError:
        return set()

    if not metal_ids:
        return set()

    atom_types = {}
    try:
        with open("NEW_COMPLEX.mol2", "r") as handle:
            in_atom = False
            for line in handle:
                if line.startswith("@<TRIPOS>ATOM"):
                    in_atom = True
                    continue
                if line.startswith("@<TRIPOS>") and in_atom:
                    break
                if in_atom:
                    parts = line.split()
                    if len(parts) >= 6:
                        try:
                            atom_types[int(parts[0])] = parts[5]
                        except ValueError:
                            pass
    except FileNotFoundError:
        return set()

    return {atom_types[i] for i in metal_ids if i in atom_types}


METAL_TYPES = read_metal_atom_types()


def key_parts(key):
    return [part.strip() for part in key.split("-") if part.strip()]


def is_sulfur_like(atom_type):
    """Preserve atom types that visibly represent sulfur in this EP2/easyPARM context."""
    t = atom_type.strip()
    tl = t.lower()
    if not t:
        return False
    # Explicitly avoid common non-sulfur element-like atom types that also start with s.
    if tl in {"se", "si", "sn", "sb", "sc", "sr", "sm"}:
        return False
    return t.startswith("S") or t.startswith("s")


def is_selected_unaveraged_key(key):
    parts = key_parts(key)
    for part in parts:
        pl = part.lower()
        if is_sulfur_like(part):
            return True
        if pl in {"ce", "os"}:
            return True
        if part in METAL_TYPES:
            return True
    return False


def mark_unaveraged_line(line, section_name, reason="policy=all_unaveraged"):
    parsed = parse_parameter_line(line)
    if parsed is None:
        return None
    key, value1, value2, comment = parsed
    atoms = extract_atoms_from_comment(comment)
    base = f"{key} {value1:.2f} {value2:.3f}"
    source = "HESSIAN_UNAVERAGED"
    old_comment = f"; upstream_comment={comment}" if comment else ""
    atom_text = f" atoms={atoms}" if atoms else ""
    new_comment = (
        f"# ep2_source={source} n=1{atom_text} source_file=forcefield.dat "
        f"{reason}{old_comment}"
    )
    append_provenance(section_name, key, atoms, f"{value1:.2f}", f"{value2:.3f}", source, "forcefield.dat", reason)
    return f"{base}    {new_comment}\n"


def group_for_averaging(lines):
    grouped = defaultdict(lambda: {"sum1": 0.0, "sum2": 0.0, "count": 0, "atoms": []})
    for line in lines:
        parsed = parse_parameter_line(line)
        if parsed is None:
            continue
        key, value1, value2, comment = parsed
        grouped[key]["sum1"] += value1
        grouped[key]["sum2"] += value2
        grouped[key]["count"] += 1
        atoms = extract_atoms_from_comment(comment)
        if atoms:
            grouped[key]["atoms"].append(atoms)
    return grouped


def write_averaged_group(file_handle, key, data, section_name, reason):
    count = data["count"]
    avg1 = data["sum1"] / count
    avg2 = data["sum2"] / count
    atoms_joined = "|".join(data["atoms"]) if data["atoms"] else ""
    source = "HESSIAN_UNAVERAGED" if count == 1 else "HESSIAN_AVERAGED"
    atom_text = f" atoms={atoms_joined}" if atoms_joined else ""
    comment = f"# ep2_source={source} n={count}{atom_text} source_file=forcefield.dat {reason}"
    file_handle.write(f"{key} {avg1:.2f} {avg2:.3f}    {comment}\n")
    append_provenance(section_name, key, atoms_joined, f"{avg1:.2f}", f"{avg2:.3f}", source, "forcefield.dat", reason)


def write_processed_section(file_handle, lines, section_name):
    if POLICY == "all_unaveraged":
        for line in lines:
            marked = mark_unaveraged_line(line, section_name, "policy=all_unaveraged")
            if marked is not None:
                file_handle.write(marked)
        return

    if POLICY == "average_all":
        grouped = group_for_averaging(lines)
        for key, data in grouped.items():
            write_averaged_group(file_handle, key, data, section_name, "policy=average_all")
        return

    # selected_unaveraged: pass through sulfur/ce/os/metal terms, average the rest.
    average_bucket = []
    for line in lines:
        parsed = parse_parameter_line(line)
        if parsed is None:
            continue
        key = parsed[0]
        if is_selected_unaveraged_key(key):
            marked = mark_unaveraged_line(line, section_name, "policy=selected_unaveraged; selector=s_or_ce_or_os_or_metal")
            if marked is not None:
                file_handle.write(marked)
        else:
            average_bucket.append(line)

    grouped = group_for_averaging(average_bucket)
    for key, data in grouped.items():
        write_averaged_group(file_handle, key, data, section_name, "policy=selected_unaveraged; selector=nonselected_averaged")


with open("forcefield.dat", "r") as file:
    content = file.readlines()

bond_start = bond_end = angle_start = angle_end = dihe_start = dihe_end = None

for i, line in enumerate(content):
    if line.strip() == "BOND":
        bond_start = i + 1
    elif line.strip() == "ANGLE":
        bond_end = i
        angle_start = i + 1
    elif line.strip() == "DIHE":
        angle_end = i
        dihe_start = i
        for j, dihe_line in enumerate(content[i:]):
            if dihe_line.strip() == "IMPROPER":
                dihe_end = i + j
                break

if None in (bond_start, bond_end, angle_start, angle_end, dihe_start, dihe_end):
    raise RuntimeError("Could not locate expected BOND/ANGLE/DIHE/IMPROPER sections in forcefield.dat")

with open("forcefield2.dat", "w") as file:
    for line in content[:bond_start]:
        file.write(line)

    write_processed_section(file, content[bond_start:bond_end], "BOND")
    file.write("\n")

    for line in content[bond_end:angle_start]:
        file.write(line)

    write_processed_section(file, content[angle_start:angle_end], "ANGLE")
    file.write("\n")

    # DIHE section: preserve existing lines. In the current EP2 chain these are
    # generic placeholder torsions from 06_get_atom_type.py, not Hessian torsion barriers.
    for line in content[dihe_start:dihe_end]:
        file.write(line)
    file.write("\n")

    for line in content[dihe_end:]:
        file.write(line)

print(f"ep2: 07_Seminario_forcefield.py completed with EP2_HESSIAN_BOND_ANGLE_POLICY={POLICY}")
