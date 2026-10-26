###################################################################################################################
#   Automated Force Fields for Metals     /$$$$$$$   /$$$$$$  /$$$$$$$  /$$      /$$                              #
#                                        | $$__  $$ /$$__  $$| $$__  $$| $$$    /$$$                              #
#   /$$$$$$   /$$$$$$   /$$$$$$$ /$$   /$$| $$  \ $$| $$  \ $$| $$  \ $$| $$$$  /$$$$                             #
#  /$$__  $$ |____  $$ /$$_____/| $$  | $$| $$$$$$$/| $$$$$$$$| $$$$$$$/| $$ $$/$$ $$                             #
# | $$$$$$$$  /$$$$$$$|  $$$$$$ | $$  | $$| $$____/ | $$__  $$| $$__  $$| $$  $$$| $$                             #
# | $$_____/ /$$__  $$ \____  $$| $$  | $$| $$      | $$  | $$| $$  \ $$| $$\  $ | $$                             #
# |  $$$$$$$|  $$$$$$$ /$$$$$$$/|  $$$$$$$| $$      | $$  | $$| $$  | $$| $$ \/  | $$                             #
#  \_______/ \_______/|_______/  \____  $$|__/      |__/  |__/|__/     |__/|__/     |__/                          #
#                               /$$  | $$                                                                         #
#                              |  $$$$$$/              ep2 modified postclean                                      #
#                               \______/                                                                          #
###################################################################################################################
#
# ep2 modification:
#   Preserve W-centered ANGLE entries instead of averaging them by atom-type key.
#
# Why:
#   The original postclean script averages repeated ANGLE entries that share the same
#   atom-type key, e.g. s1-W-s2. That is acceptable for generic equivalent organic
#   terms, but problematic for anisotropic metal sites. Different W-centered atom-index
#   angles can collapse to the same type key and lose their Java/XYZViewer-like
#   Seminario values.
#
# Policy implemented here:
#   BOND  section: same as original; duplicate keys are averaged.
#   ANGLE section:
#       - if central atom type is W, preserve the original line and do NOT average it.
#       - otherwise, duplicate keys are averaged as in the original script.
#   DIHE  section: same as original; duplicate keys are averaged.
#
# Note:
#   Preserving W-centered duplicate type keys may create repeated ANGLE records with
#   the same key. For production Amber parameters, the better long-term solution is
#   to split atom types upstream so distinct W-centered angles have distinct keys.
###################################################################################################################

import re


PRESERVE_W_CENTERED_ANGLES = True


def process_frcmod_file(file_path, atomtype_file):
    """
    Process an frcmod-like file and write updated_<file_path>.

    Original behavior:
        average duplicate BOND, ANGLE, and DIHE entries by atom-type key.

    ep2 behavior:
        average duplicate BOND and DIHE entries;
        average duplicate ANGLE entries only when they are not W-centered.

    Parameters
    ----------
    file_path : str
        Input frcmod file, e.g. 'updated_COMPLEX_modified2.frcmod'

    atomtype_file : str
        Atom-type file. Read for compatibility with the original script.
        The current script does not actively use the atom_types set.
    """
    with open(atomtype_file, 'r') as at_file:
        atom_types = set(line.strip() for line in at_file)

    with open(file_path, 'r') as file:
        lines = file.readlines()

    bond_section = False
    angle_section = False
    dihe_section = False

    data_bond = {}
    data_angle = {}
    data_dihe = {}

    # Accept integer or decimal numeric fields.
    num = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)"
    pattern_bond = re.compile(rf"^(\S+)-(\S+)\s+({num})\s+({num})$")
    pattern_angle = re.compile(rf"^(\S+)-(\S+)-(\S+)\s+({num})\s+({num})$")
    pattern_dihe = re.compile(rf"^(\S+)-(\S+)-(\S+)-(\S+)\s+(\d+)\s+({num})\s+({num})\s+({num})$")

    # First pass: collect averages.
    # W-centered angles are intentionally excluded from averaging.
    for raw_line in lines:
        line = raw_line.strip()

        if line == "BOND":
            bond_section = True
            angle_section = False
            dihe_section = False
            continue
        elif line == "ANGLE":
            bond_section = False
            angle_section = True
            dihe_section = False
            continue
        elif line == "DIHE":
            bond_section = False
            angle_section = False
            dihe_section = True
            continue
        elif line in {"IMPROPER", "NONBON", "MASS"}:
            bond_section = False
            angle_section = False
            dihe_section = False
            continue

        if bond_section:
            match = pattern_bond.match(line)
            if match:
                col1, col2, col3, col4 = match.groups()
                col3, col4 = float(col3), float(col4)

                key = f"{col1}-{col2}"
                if key not in data_bond:
                    data_bond[key] = {'sum_col3': 0.0, 'sum_col4': 0.0, 'count': 0}
                data_bond[key]['sum_col3'] += col3
                data_bond[key]['sum_col4'] += col4
                data_bond[key]['count'] += 1

        elif angle_section:
            match = pattern_angle.match(line)
            if match:
                col1, col2, col3, col4, col5 = match.groups()
                col4, col5 = float(col4), float(col5)

                is_w_centered = (col2 == "W")
                if PRESERVE_W_CENTERED_ANGLES and is_w_centered:
                    # Do not collect W-centered angles for averaging.
                    continue

                key = f"{col1}-{col2}-{col3}"
                if key not in data_angle:
                    data_angle[key] = {'sum_col4': 0.0, 'sum_col5': 0.0, 'count': 0}
                data_angle[key]['sum_col4'] += col4
                data_angle[key]['sum_col5'] += col5
                data_angle[key]['count'] += 1

        elif dihe_section:
            match = pattern_dihe.match(line)
            if match:
                col1, col2, col3, col4, col5, col6, col7, col8 = match.groups()
                col5 = int(col5)
                col6, col7, col8 = float(col6), float(col7), float(col8)

                key = f"{col1}-{col2}-{col3}-{col4}"
                if key not in data_dihe:
                    data_dihe[key] = {
                        'sum_col5': 0,
                        'sum_col6': 0.0,
                        'sum_col7': 0.0,
                        'sum_col8': 0.0,
                        'count': 0
                    }
                data_dihe[key]['sum_col5'] += col5
                data_dihe[key]['sum_col6'] += col6
                data_dihe[key]['sum_col7'] += col7
                data_dihe[key]['sum_col8'] += col8
                data_dihe[key]['count'] += 1

    output_file_path = 'updated_' + file_path

    with open(output_file_path, 'w') as output_file:
        bond_section = False
        angle_section = False
        dihe_section = False

        for raw_line in lines:
            line = raw_line.strip()

            if line == "BOND":
                bond_section = True
                angle_section = False
                dihe_section = False
                output_file.write(line + "\n")
                continue
            elif line == "ANGLE":
                bond_section = False
                angle_section = True
                dihe_section = False
                output_file.write(line + "\n")
                continue
            elif line == "DIHE":
                bond_section = False
                angle_section = False
                dihe_section = True
                output_file.write(line + "\n")
                continue
            elif line in {"IMPROPER", "NONBON", "MASS"}:
                bond_section = False
                angle_section = False
                dihe_section = False
                output_file.write(line + "\n")
                continue

            if bond_section:
                match = pattern_bond.match(line)
                if match:
                    col1, col2, _, _ = match.groups()
                    key = f"{col1}-{col2}"

                    if key in data_bond:
                        sum_col3 = data_bond[key]['sum_col3']
                        sum_col4 = data_bond[key]['sum_col4']
                        count = data_bond[key]['count']

                        avg_col3 = sum_col3 / count
                        avg_col4 = sum_col4 / count

                        output_file.write(f"{col1}-{col2} {avg_col3:.2f} {avg_col4:.3f}\n")
                        data_bond.pop(key)
                else:
                    output_file.write(line + "\n")

            elif angle_section:
                match = pattern_angle.match(line)
                if match:
                    col1, col2, col3, col4, col5 = match.groups()
                    key = f"{col1}-{col2}-{col3}"

                    is_w_centered = (col2 == "W")

                    if PRESERVE_W_CENTERED_ANGLES and is_w_centered:
                        # ep2 logic:
                        # Preserve the original W-centered angle term and do not average.
                        output_file.write(
                            f"{col1}-{col2}-{col3} {float(col4):.3f} {float(col5):.3f}\n"
                        )
                        continue

                    if key in data_angle:
                        sum_col4 = data_angle[key]['sum_col4']
                        sum_col5 = data_angle[key]['sum_col5']
                        count = data_angle[key]['count']

                        avg_col4 = sum_col4 / count
                        avg_col5 = sum_col5 / count

                        output_file.write(f"{col1}-{col2}-{col3} {avg_col4:.3f} {avg_col5:.3f}\n")
                        data_angle.pop(key)
                else:
                    output_file.write(line + "\n")

            elif dihe_section:
                match = pattern_dihe.match(line)
                if match:
                    col1, col2, col3, col4, _, _, _, _ = match.groups()
                    key = f"{col1}-{col2}-{col3}-{col4}"

                    if key in data_dihe:
                        sum_col5 = data_dihe[key]['sum_col5']
                        sum_col6 = data_dihe[key]['sum_col6']
                        sum_col7 = data_dihe[key]['sum_col7']
                        sum_col8 = data_dihe[key]['sum_col8']
                        count = data_dihe[key]['count']

                        avg_col5 = round(sum_col5 / count)
                        avg_col6 = sum_col6 / count
                        avg_col7 = sum_col7 / count
                        avg_col8 = sum_col8 / count

                        output_file.write(
                            f"{col1}-{col2}-{col3}-{col4}   "
                            f"{avg_col5:d}    {avg_col6:.3f}       "
                            f"{avg_col7:.4f}           {avg_col8:.3f}\n"
                        )
                        data_dihe.pop(key)
                else:
                    output_file.write(line + "\n")

            else:
                output_file.write(line + "\n")

    print(f"Wrote {output_file_path}")


if __name__ == "__main__":
    process_frcmod_file('updated_COMPLEX_modified2.frcmod', 'new_atomtype.dat')
