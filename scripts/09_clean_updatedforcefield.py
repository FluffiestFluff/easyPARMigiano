###################################################################################################################
#   Automated Force Fields for Metals / easyPARM
#   EP2 provenance-enabled variant of 09_clean_updatedforcefield.py
#
#   This version prefers the Hessian/Seminario BOND and ANGLE stream from
#   forcefield2.dat over parmchk/GAFF lines whenever the same atom-type key is
#   present. It also removes parmchk metal terms and reinserts Hessian metal
#   DIHE/IMPROPER terms as in the original workflow.
###################################################################################################################

import periodictable
import re
from pathlib import Path

PROVENANCE_FILE = "ep2_parameter_provenance.tsv"
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


def parse_ep2_source(line):
    match = re.search(r"ep2_source=([A-Za-z0-9_]+)", line)
    return match.group(1) if match else "HESSIAN_INSERTED_FROM_FORCEFIELD2"


def parse_atoms(line):
    match = re.search(r"atoms=([^\s#;]+)", line)
    return match.group(1) if match else ""


def parse_key_values(line):
    before_comment = line.split("#", 1)[0]
    fields = before_comment.split()
    if len(fields) < 3:
        return "", "", ""
    return fields[0], fields[1], fields[2]


def key_parts_from_line(line):
    key, _, _ = parse_key_values(line)
    if not key or "-" not in key:
        return []
    return [part.strip() for part in key.split("-") if part.strip()]


def normalize_key(key):
    parts = [part.strip() for part in key.split("-") if part.strip()]
    if len(parts) == 2:
        return "-".join(sorted(parts, key=str.lower))
    if len(parts) == 3:
        # Angles A-B-C are equivalent to C-B-A; preserve the central atom.
        outer = sorted([parts[0], parts[2]], key=str.lower)
        return f"{outer[0]}-{parts[1]}-{outer[1]}"
    return "-".join(parts)


def read_metal_numbers(file_path):
    try:
        with open(file_path, "r") as file:
            return [int(line.strip()) for line in file if line.strip()]
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return []


def get_element_from_atom_name(atom_name):
    match = re.match(r"([A-Za-z]+)", atom_name)
    if match:
        element_name = match.group(1)
        if len(element_name) >= 2:
            possible_element = element_name[:2].capitalize()
            if possible_element in {elem.symbol for elem in periodictable.elements}:
                return possible_element
        possible_element = element_name[0].upper()
        if possible_element in {elem.symbol for elem in periodictable.elements}:
            return possible_element
    return None


def extract_atom_types_from_mol2(file_path):
    try:
        with open(file_path, "r") as file:
            lines = file.readlines()
        atom_info = {}
        is_atom_section = False

        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                is_atom_section = True
                continue
            if line.startswith("@<TRIPOS>"):
                is_atom_section = False
                continue
            if is_atom_section:
                parts = line.split()
                if len(parts) >= 6:
                    atom_id = int(parts[0])
                    full_atom_name = parts[1]
                    atom_type = parts[5]
                    atom_info[atom_id] = {"name": full_atom_name, "type": atom_type}

        return atom_info
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return {}


def create_atom_type_to_element_mapping(atom_info):
    atom_type_to_element = {}
    atom_type_to_name = {}

    for atom_data in atom_info.values():
        atom_type = atom_data["type"]
        atom_name = atom_data["name"]
        atom_type_to_name.setdefault(atom_type, set()).add(atom_name)

    for atom_type, atom_names in atom_type_to_name.items():
        elements = {get_element_from_atom_name(name) for name in atom_names if get_element_from_atom_name(name)}
        if len(elements) == 1:
            element = elements.pop()
            if element:
                atom_type_to_element[atom_type] = element

    return atom_type_to_element


def make_atom_type_patterns(atom_types):
    patterns = []
    for atom_type in atom_types:
        pattern = rf"(^|\s|\-|\.)({re.escape(atom_type)})($|\s|\-|\.)"
        patterns.append(pattern)
    return patterns


def line_contains_any_atom_type(line, atom_types):
    stripped = line.strip()
    for pattern in make_atom_type_patterns(atom_types):
        if re.search(pattern, stripped):
            return True
    return False


def remove_metal_lines(input_file, output_file, metal_atom_types):
    """Remove parmchk-derived metal terms from BOND/ANGLE/DIHE/IMPROPER."""
    try:
        with open(input_file, "r") as infile, open(output_file, "w") as outfile:
            current_section = None
            removed = 0
            for line in infile:
                stripped_line = line.strip()

                if stripped_line in SECTION_HEADERS:
                    current_section = stripped_line
                    outfile.write(line)
                    continue

                if not stripped_line:
                    outfile.write(line)
                    continue

                if current_section in {"BOND", "ANGLE", "DIHE", "IMPROPER"} and line_contains_any_atom_type(line, metal_atom_types):
                    removed += 1
                    continue

                outfile.write(line)
        print(f"ep2: removed {removed} parmchk/GAFF metal-containing bonded term(s).")
    except Exception as e:
        print(f"Error updating FRCMOD file: {str(e)}")
        raise


def update_mass_section(input_file, output_file, metal_atom_types, atom_info):
    try:
        atom_type_to_element = create_atom_type_to_element_mapping(atom_info)
        mol2_atom_types = {data["type"] for data in atom_info.values()}

        with open(input_file, "r") as infile, open(output_file, "w") as outfile:
            in_mass_section = False

            for line in infile:
                stripped_line = line.strip()

                if stripped_line == "MASS":
                    in_mass_section = True
                    outfile.write(line)
                    continue

                if in_mass_section:
                    if stripped_line in SECTION_HEADERS - {"MASS"}:
                        in_mass_section = False
                        outfile.write(line)
                        continue

                    if not stripped_line:
                        outfile.write(line)
                        continue

                    parts = line.split()
                    if len(parts) >= 3:
                        atom_type = parts[0]
                        if atom_type not in mol2_atom_types:
                            continue
                        try:
                            current_mass = float(parts[1])
                        except ValueError:
                            outfile.write(line)
                            continue

                        element_symbol = atom_type_to_element.get(atom_type)
                        if element_symbol:
                            element = getattr(periodictable, element_symbol, None)
                            if element:
                                correct_mass = element.mass
                                if abs(current_mass - correct_mass) > 0.1:
                                    mass_start_pos = line.find(str(current_mass))
                                    if mass_start_pos != -1:
                                        mass_end_pos = mass_start_pos + len(str(current_mass))
                                        new_line = (
                                            f"{line[:mass_start_pos]}"
                                            f"{correct_mass:8.3f}"
                                            f"{line[mass_end_pos:]}"
                                        )
                                        outfile.write(new_line)
                                        continue

                    outfile.write(line)
                else:
                    outfile.write(line)

    except Exception as e:
        print(f"Error updating MASS section: {str(e)}")
        raise


def collect_hessian_bond_angle_lines(forcefield_file):
    """Collect all Hessian/Seminario BOND and ANGLE lines from forcefield2.dat."""
    section_lines = {"BOND": [], "ANGLE": []}
    replacement_keys = {"BOND": set(), "ANGLE": set()}
    current_section = None

    with open(forcefield_file, "r") as infile:
        for line in infile:
            stripped = line.strip()
            if stripped in SECTION_HEADERS:
                current_section = stripped
                continue
            if current_section not in {"BOND", "ANGLE"} or not stripped:
                continue

            key, v1, v2 = parse_key_values(line)
            if not key:
                continue

            # forcefield2.dat is the Hessian/Seminario reference stream. Prefer explicit
            # HESSIAN labels, but accept unlabeled lines for compatibility with old easyPARM.
            if "ep2_source=" in line and "ep2_source=HESSIAN" not in line:
                continue

            section_lines[current_section].append(line if line.endswith("\n") else line + "\n")
            replacement_keys[current_section].add(normalize_key(key))

            append_provenance(
                current_section,
                key,
                parse_atoms(line),
                v1,
                v2,
                parse_ep2_source(line),
                forcefield_file,
                "collected_by=09_clean_updatedforcefield.py; hessian_preferred_over_parmchk=yes",
            )

    return section_lines, replacement_keys


def replace_bond_angle_lines_with_hessian(forcefield_file, frcmod_file):
    """
    Remove parmchk/GAFF BOND/ANGLE lines whose keys are present in forcefield2.dat,
    then insert the forcefield2.dat Hessian/Seminario BOND/ANGLE lines.
    """
    hessian_lines, replacement_keys = collect_hessian_bond_angle_lines(forcefield_file)
    if not hessian_lines["BOND"] and not hessian_lines["ANGLE"]:
        print("ep2: no Hessian/Seminario BOND/ANGLE replacements found in forcefield2.dat; keeping parmchk terms.")
        return

    with open(frcmod_file, "r") as infile:
        lines = infile.readlines()

    output = []
    current_section = None
    inserted = {"BOND": False, "ANGLE": False}
    removed = {"BOND": 0, "ANGLE": 0}

    for line in lines:
        stripped = line.strip()

        if stripped in SECTION_HEADERS:
            current_section = stripped
            output.append(line)
            if stripped in {"BOND", "ANGLE"} and not inserted[stripped]:
                for repl in hessian_lines[stripped]:
                    output.append(repl)
                inserted[stripped] = True
            continue

        if current_section in {"BOND", "ANGLE"}:
            key, _, _ = parse_key_values(line)
            if key and normalize_key(key) in replacement_keys[current_section]:
                removed[current_section] += 1
                continue

        output.append(line)

    with open(frcmod_file, "w") as outfile:
        outfile.writelines(output)

    print(
        "ep2: Hessian/Seminario override: "
        f"removed {removed['BOND']} parmchk BOND line(s), inserted {len(hessian_lines['BOND'])}; "
        f"removed {removed['ANGLE']} parmchk ANGLE line(s), inserted {len(hessian_lines['ANGLE'])}."
    )


def copy_metal_dihe_improper_lines(forcefield_file, output_file, metal_atom_types):
    """Copy metal-containing DIHE/IMPROPER terms from forcefield2.dat after BOND/ANGLE override."""
    try:
        with open(forcefield_file, "r") as infile:
            forcefield_lines = infile.readlines()
        with open(output_file, "r") as outfile:
            updated_lines = outfile.readlines()

        section_lines = {"DIHE": [], "IMPROPER": []}
        current_section = None
        seen = {"DIHE": set(), "IMPROPER": set()}

        for line in forcefield_lines:
            stripped_line = line.strip()
            if stripped_line in SECTION_HEADERS:
                current_section = stripped_line
                continue
            if current_section not in section_lines or not stripped_line:
                continue
            if line_contains_any_atom_type(line, metal_atom_types):
                key, v1, v2 = parse_key_values(line)
                dedup = line.split("#", 1)[0].strip()
                if dedup in seen[current_section]:
                    continue
                seen[current_section].add(dedup)
                section_lines[current_section].append(line if line.endswith("\n") else line + "\n")
                append_provenance(
                    current_section,
                    key,
                    parse_atoms(line),
                    v1,
                    v2,
                    parse_ep2_source(line),
                    forcefield_file,
                    "inserted_by=09_clean_updatedforcefield.py; replaced_parmchk_metal_dihe_or_improper=yes",
                )

        final_lines = []
        inserted = {"DIHE": False, "IMPROPER": False}
        for line in updated_lines:
            stripped_line = line.strip()
            final_lines.append(line)
            if stripped_line in section_lines and not inserted[stripped_line]:
                final_lines.extend(section_lines[stripped_line])
                inserted[stripped_line] = True

        with open(output_file, "w") as outfile:
            outfile.writelines(final_lines)

    except Exception as e:
        print(f"Error inserting metal DIHE/IMPROPER lines into the correct sections: {str(e)}")
        raise


def main():
    try:
        atom_info = extract_atom_types_from_mol2("NEW_COMPLEX.mol2")
        if not atom_info:
            print("Error: No atom information extracted from mol2 file.")
            return

        metal_ids = read_metal_numbers("metal_number.dat")
        if not metal_ids:
            print("Error: No metal atom IDs found.")
            return

        metal_atom_types = {atom_info[metal_id]["type"] for metal_id in metal_ids if metal_id in atom_info}
        if not metal_atom_types:
            print("Error: Could not determine metal atom types.")
            return

        remove_metal_lines("updated_COMPLEX_modified.frcmod", "temp_COMPLEX_modified.frcmod", metal_atom_types)
        update_mass_section("temp_COMPLEX_modified.frcmod", "updated_COMPLEX_modified2.frcmod", metal_atom_types, atom_info)

        # Prefer Hessian/Seminario BOND and ANGLE terms globally.
        replace_bond_angle_lines_with_hessian("forcefield2.dat", "updated_COMPLEX_modified2.frcmod")

        # Preserve original metal DIHE/IMPROPER override behavior without duplicating BOND/ANGLE.
        copy_metal_dihe_improper_lines("forcefield2.dat", "updated_COMPLEX_modified2.frcmod", metal_atom_types)

    except Exception as e:
        print(f"An error occurred in the main function: {str(e)}")
        print("Traceback:")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
