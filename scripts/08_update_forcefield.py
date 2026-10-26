###################################################################################################################
#   Automated Force Fields for Metals / easyPARM
#   EP2 provenance-enabled variant of 08_update_forcefield.py
###################################################################################################################

import itertools
import re
from pathlib import Path

PROVENANCE_FILE = "ep2_parameter_provenance.tsv"


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


def add_or_extend_ep2_comment(line, comment_payload):
    if "ep2_source=" in line:
        return line
    if "#" in line:
        return line + " ; " + comment_payload
    return line + " # " + comment_payload


def classify_parmchk_line(line):
    lower = line.lower()
    if "same as" in lower and "penalty score" in lower:
        source = "PARMCHK_GAFF_ANALOGY"
    elif "attn" in lower or "need revision" in lower:
        source = "PARMCHK_ATTENTION_NEEDS_REVISION"
    else:
        source = "PARMCHK_GAFF_OR_FORCEFIELD"

    same_as = ""
    penalty = ""

    same_match = re.search(r"same\s+as\s+([^,]+)", line, flags=re.IGNORECASE)
    if same_match:
        same_as = same_match.group(1).strip()

    penalty_match = re.search(r"penalty\s+score\s*=\s*([-+0-9.]+)", line, flags=re.IGNORECASE)
    if penalty_match:
        penalty = penalty_match.group(1).strip()

    notes = []
    if same_as:
        notes.append(f"same_as={same_as.replace(' ', '')}")
    if penalty:
        notes.append(f"penalty={penalty}")

    return source, ";".join(notes)


def parse_key_values(line):
    before_comment = line.split("#", 1)[0]
    fields = before_comment.split()
    if len(fields) < 3:
        return "", "", ""
    return fields[0], fields[1], fields[2]


# Extract atom type of mol2.
def extract_atom_types_from_mol2(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        atom_types = {}
        is_atom_section = False
        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                is_atom_section = True
                continue
            if line.startswith("@<TRIPOS>BOND"):
                is_atom_section = False
                break
            if is_atom_section:
                parts = line.split()
                if len(parts) >= 6:
                    atom_id = int(parts[0])
                    atom_type = parts[5]
                    atom_types[atom_id] = atom_type
        return atom_types
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return {}


# Function to read metal numbers from a file.
def read_metal_numbers(file_path):
    try:
        with open(file_path, 'r') as file:
            return [int(line.strip()) for line in file if line.strip()]
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return []


# Function to read bond distances and create a dictionary of bonds.
def read_distance_file(file_path):
    bonds = {}
    try:
        with open(file_path, 'r') as file:
            for line in file:
                parts = line.split()
                if len(parts) == 3:
                    atom1, atom2, distance = map(float, parts)
                    bonds[int(atom1), int(atom2)] = distance
        return bonds
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return {}


# Replace the old atom type with the new atom type.
def replace_atom_type(line, old_to_new_types):
    match = re.match(r'^(\S+(?:\s*-\s*\S+)*)(.*)$', line)
    if not match:
        return [line]

    atom_types, rest_of_line = match.groups()
    atom_parts = [part.strip() for part in re.split(r'\s*-\s*', atom_types)]

    combinations = []
    for part in atom_parts:
        if part in old_to_new_types:
            combinations.append(old_to_new_types[part])
        else:
            combinations.append([part])

    all_combinations = list(itertools.product(*combinations))

    new_lines = []
    for comb in all_combinations:
        new_atom_types = '-'.join(comb)
        new_line = f"{new_atom_types}{rest_of_line}"
        new_lines.append(new_line)

    return new_lines


# Update the frcmod with the new atom type.
def update_frcmod_file(input_file, output_file, old_to_new_types, metal_specific_types):
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            current_section = None
            for line in infile:
                original_line = line.rstrip('\n')

                if original_line in ["MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"]:
                    current_section = original_line
                    outfile.write(original_line + '\n')
                    continue

                if not original_line or original_line.startswith('Remark') or original_line.startswith('#'):
                    outfile.write(original_line + '\n')
                    continue

                if current_section == "NONBON":
                    parts = original_line.split()
                    if parts and parts[0] in old_to_new_types:
                        for new_type in old_to_new_types[parts[0]]:
                            new_line = original_line.replace(parts[0], new_type, 1)
                            outfile.write(new_line + '\n')
                    else:
                        outfile.write(original_line + '\n')
                else:
                    new_lines = replace_atom_type(original_line, old_to_new_types)
                    for new_line in new_lines:
                        if current_section in {"BOND", "ANGLE", "DIHE"}:
                            source, notes = classify_parmchk_line(original_line)
                            key, v1, v2 = parse_key_values(new_line)
                            payload = f"ep2_source={source} source_file={input_file}"
                            if notes:
                                payload += f" {notes}"
                            new_line = add_or_extend_ep2_comment(new_line, payload)
                            append_provenance(current_section, key, "NA", v1, v2, source, input_file, notes)
                        outfile.write(new_line + '\n')

    except Exception as e:
        print(f"Error updating FRCMOD file: {str(e)}")


# Main function to run the replacement.
def main():
    old_atom_types = extract_atom_types_from_mol2('COMPLEX.mol2')
    new_atom_types = extract_atom_types_from_mol2('NEW_COMPLEX.mol2')
    metal_ids = read_metal_numbers('metal_number.dat')
    bonds = read_distance_file('distance.dat')

    metal_bonded_atoms = {metal_id: set() for metal_id in metal_ids}
    for (atom1, atom2), distance in bonds.items():
        if atom1 in metal_ids:
            metal_bonded_atoms[atom1].add(atom2)
        if atom2 in metal_ids:
            metal_bonded_atoms[atom2].add(atom1)

    old_to_new_types = {}
    metal_specific_types = {metal_id: set() for metal_id in metal_ids}

    # Create a mapping of old atom types to all possible new atom types.
    for atom_id, old_type in old_atom_types.items():
        if atom_id in new_atom_types:
            new_type = new_atom_types[atom_id]
            if old_type not in old_to_new_types:
                old_to_new_types[old_type] = set()
            old_to_new_types[old_type].add(new_type)

            for metal_id, bonded_atoms in metal_bonded_atoms.items():
                if atom_id in bonded_atoms:
                    metal_specific_types[metal_id].add(new_type)

    # Include original types in the replacement set.
    for old_type in list(old_to_new_types.keys()):
        old_to_new_types[old_type].add(old_type)

    old_to_new_types = {k: list(v) for k, v in old_to_new_types.items()}

    update_frcmod_file('COMPLEX_modified.frcmod', 'updated_COMPLEX_modified.frcmod', old_to_new_types, metal_specific_types)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"An error occurred in the main function: {str(e)}")
        print("Traceback:")
        import traceback
        traceback.print_exc()
