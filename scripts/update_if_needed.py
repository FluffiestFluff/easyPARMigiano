###################################################################################################################
#   Automated Force Fields for Metals / easyPARM update_if_needed.py-derived helper
#   EP2-local copy: same scaling rules, with explicit provenance logging.
###################################################################################################################

import os
import re
import shutil
import sys
from pathlib import Path

PROVENANCE_FILE = "ep2_parameter_provenance.tsv"

BOND_SCALE_RULES = [
    (20, 4.599),
]

ANGLE_SCALE_RULES = [
    (5, 11.599),
    (10, 7.799),
    (20, 3.599),
    (29, 2.699),
]

LINE_RE = re.compile(
    r"^([A-Za-z0-9*']+"
    r"(?:\s*-\s*[A-Za-z0-9*']+)+"
    r")"
    r"(\s+)([-0-9.]+)"
    r"(\s+)([-0-9.]+)"
    r"(.*)"
)


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


def apply_scale(force, rules):
    """Return (scaled_force, factor) according to the first matching rule, or unchanged."""
    for threshold, factor in rules:
        if force < threshold:
            return force * factor, factor
    return force, 1.0


def add_or_extend_ep2_comment(trailing, payload):
    if "#" in trailing:
        return trailing.rstrip() + " ; " + payload
    if trailing.strip():
        return trailing.rstrip() + " # " + payload
    return " # " + payload


def process_frcmod(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()

    section = None
    output = []
    updated = 0
    unchanged = 0
    skipped = 0

    for line in lines:
        stripped = line.rstrip("\n")

        if re.match(r"^BOND", stripped):
            section = "BOND"
            output.append(line)
            continue
        elif re.match(r"^ANGLE", stripped):
            section = "ANGLE"
            output.append(line)
            continue
        elif re.match(r"^DIHE|^IMPROPER|^NONBON", stripped):
            section = stripped[:6]
            output.append(line)
            continue

        if not stripped or stripped.startswith("#"):
            output.append(line)
            continue

        if section in ("BOND", "ANGLE") and re.match(r"^[A-Za-z]", stripped):
            match = LINE_RE.match(stripped)
            if match:
                atom_tuple = match.group(1)
                sep1 = match.group(2)
                force_str = match.group(3)
                sep2 = match.group(4)
                eq_val = match.group(5)
                trailing = match.group(6)

                force = float(force_str)
                if section == "BOND":
                    new_force, factor = apply_scale(force, BOND_SCALE_RULES)
                else:
                    new_force, factor = apply_scale(force, ANGLE_SCALE_RULES)

                formatted_force = f"{new_force:.2f}"

                if new_force != force:
                    payload = (
                        "ep2_scaled=EASYPARM_WEAK_FORCE "
                        f"old_k={force_str} new_k={formatted_force} factor={factor:.3f} "
                        f"scaling_script={os.path.basename(__file__)}"
                    )
                    trailing = add_or_extend_ep2_comment(trailing, payload)
                    print(
                        "  [%s] %-24s force: %s -> %s"
                        % (section, atom_tuple.strip(), force_str, formatted_force)
                    )
                    append_provenance(
                        section,
                        atom_tuple.strip(),
                        "NA",
                        formatted_force,
                        eq_val,
                        "EASYPARM_WEAK_FORCE_SCALED",
                        os.path.basename(filepath),
                        f"old_k={force_str};factor={factor:.3f}",
                    )
                    updated += 1
                else:
                    unchanged += 1

                output.append(
                    f"{atom_tuple}{sep1}{formatted_force}{sep2}{eq_val}{trailing}\n"
                )
                continue
            skipped += 1

        output.append(line)

    backup = filepath + ".bak"
    shutil.copy2(filepath, backup)
    with open(filepath, "w") as f:
        f.writelines(output)

    print(f"Updated force constants: {updated}")
    print(f"Unchanged force constants: {unchanged}")
    print(f"Skipped BOND/ANGLE lines not matching parser: {skipped}")
    print(f"Backup written to: {backup}")


def main():
    if len(sys.argv) != 2:
        print("Usage: python3 update_if_needed.py <path/to/COMPLEX.frcmod>", file=sys.stderr)
        sys.exit(1)

    frcmod_path = sys.argv[1]
    if not os.path.exists(frcmod_path):
        print(f"Error: File not found: {frcmod_path}", file=sys.stderr)
        sys.exit(1)

    process_frcmod(frcmod_path)


if __name__ == "__main__":
    main()
