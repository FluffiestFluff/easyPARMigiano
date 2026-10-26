#!/usr/bin/env python3
"""
apply_custom_atom_types.py

Apply user-defined atom-type overrides to a MOL2 file.

Default use inside EP2:
    python3 apply_custom_atom_types.py custom_atom_types.dat NEW_COMPLEX.mol2 NEW_COMPLEX.mol2

Supported override formats, whitespace-separated:

    # safest: atom id from MOL2 atom column 1
    7 sD1
    atom_id 7 sD1
    id 7 sD1

    # atom name from MOL2 atom column 2
    atom_name S7 sD1
    name S7 sD1

    # global current atom-type replacement from MOL2 atom column 6
    atom_type s1 sD1
    type s1 sD1
    old_type s1 sD1

Lines beginning with # are ignored. Inline comments are allowed.

Notes:
- Atom-id overrides are recommended because atom names/types can repeat.
- This script edits only the MOL2 atom type column. Coordinates, charges, atom
  names, residue IDs, and bonds are preserved.
- A backup of the output target is written as <output>.before_custom_atom_types.bak
  when input and output are the same existing file.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Tuple


VALID_SELECTOR_KEYS = {"atom_id", "id", "atom", "atom_name", "name", "atom_type", "type", "old_type"}


def strip_inline_comment(line: str) -> str:
    return line.split("#", 1)[0].strip()


def read_overrides(path: Path):
    by_id: Dict[int, str] = {}
    by_name: Dict[str, str] = {}
    by_type: Dict[str, str] = {}

    with path.open() as handle:
        for line_no, raw in enumerate(handle, start=1):
            line = strip_inline_comment(raw)
            if not line:
                continue

            fields = line.split()
            if len(fields) == 2:
                selector, new_type = fields
                try:
                    atom_id = int(selector)
                except ValueError:
                    raise ValueError(
                        f"{path}:{line_no}: two-column format requires integer atom ID: {raw.rstrip()}"
                    )
                by_id[atom_id] = new_type
                continue

            if len(fields) < 3:
                raise ValueError(f"{path}:{line_no}: expected 2 or 3+ fields: {raw.rstrip()}")

            key = fields[0].lower()
            if key not in VALID_SELECTOR_KEYS:
                raise ValueError(
                    f"{path}:{line_no}: unknown selector '{fields[0]}'. "
                    f"Use atom_id/id, atom_name/name, or atom_type/type/old_type."
                )

            value = fields[1]
            new_type = fields[2]

            if key in {"atom_id", "id", "atom"}:
                try:
                    atom_id = int(value)
                except ValueError:
                    raise ValueError(f"{path}:{line_no}: atom_id must be an integer: {raw.rstrip()}")
                by_id[atom_id] = new_type
            elif key in {"atom_name", "name"}:
                by_name[value] = new_type
            elif key in {"atom_type", "type", "old_type"}:
                by_type[value] = new_type

    return by_id, by_name, by_type


def apply_overrides_to_mol2(input_mol2: Path, output_mol2: Path, override_file: Path) -> Tuple[int, List[str]]:
    by_id, by_name, by_type = read_overrides(override_file)
    lines = input_mol2.read_text().splitlines(keepends=True)

    output: List[str] = []
    in_atom_section = False
    changed = 0
    log_lines: List[str] = []

    for raw in lines:
        if raw.startswith("@<TRIPOS>ATOM"):
            in_atom_section = True
            output.append(raw)
            continue
        if raw.startswith("@<TRIPOS>") and in_atom_section:
            in_atom_section = False
            output.append(raw)
            continue

        if not in_atom_section or not raw.strip():
            output.append(raw)
            continue

        parts = raw.split()
        if len(parts) < 9:
            output.append(raw)
            continue

        try:
            atom_id = int(parts[0])
        except ValueError:
            output.append(raw)
            continue

        atom_name = parts[1]
        old_type = parts[5]
        new_type = None
        reason = None

        # Priority: explicit atom id > atom name > current atom type.
        if atom_id in by_id:
            new_type = by_id[atom_id]
            reason = f"atom_id={atom_id}"
        elif atom_name in by_name:
            new_type = by_name[atom_name]
            reason = f"atom_name={atom_name}"
        elif old_type in by_type:
            new_type = by_type[old_type]
            reason = f"atom_type={old_type}"

        if new_type is None or new_type == old_type:
            output.append(raw)
            continue

        if len(new_type) > 4:
            log_lines.append(
                f"WARNING atom {atom_id} {atom_name}: new type '{new_type}' is longer than 4 characters; "
                "Amber formats may become fragile."
            )

        parts[5] = new_type
        changed += 1
        log_lines.append(f"atom_id={atom_id}\tatom_name={atom_name}\told_type={old_type}\tnew_type={new_type}\tselector={reason}")

        # Reconstruct a conventional MOL2 atom line. This is intentionally simple
        # and preserves the actual values, not exact original spacing.
        try:
            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
            subst_id = int(parts[6])
            subst_name = parts[7]
            charge = float(parts[8])
            line_out = (
                f"{atom_id:7d} {atom_name:<8s} {x:10.4f} {y:10.4f} {z:10.4f} "
                f"{new_type:<6s} {subst_id:3d} {subst_name:<8s} {charge:12.6f}\n"
            )
        except Exception:
            # Fallback if unexpected MOL2 columns are present.
            line_out = " ".join(parts) + "\n"
        output.append(line_out)

    if input_mol2.resolve() == output_mol2.resolve() and output_mol2.exists():
        backup = output_mol2.with_suffix(output_mol2.suffix + ".before_custom_atom_types.bak")
        shutil.copy2(output_mol2, backup)
        log_lines.append(f"backup={backup}")

    output_mol2.write_text("".join(output))
    return changed, log_lines


def main() -> None:
    parser = argparse.ArgumentParser(description="Apply custom atom-type overrides to a MOL2 file.")
    parser.add_argument("override_file", help="Custom atom-type override file.")
    parser.add_argument("input_mol2", help="Input MOL2 file, usually NEW_COMPLEX.mol2.")
    parser.add_argument(
        "output_mol2",
        nargs="?",
        help="Output MOL2 file. Default: overwrite input_mol2.",
    )
    parser.add_argument(
        "--log",
        default="custom_atom_types_applied.tsv",
        help="Log file. Default: custom_atom_types_applied.tsv",
    )

    args = parser.parse_args()
    override_file = Path(args.override_file)
    input_mol2 = Path(args.input_mol2)
    output_mol2 = Path(args.output_mol2) if args.output_mol2 else input_mol2

    if not override_file.exists():
        raise FileNotFoundError(f"Override file not found: {override_file}")
    if not input_mol2.exists():
        raise FileNotFoundError(f"MOL2 file not found: {input_mol2}")

    changed, log_lines = apply_overrides_to_mol2(input_mol2, output_mol2, override_file)

    log_path = Path(args.log)
    with log_path.open("w") as log:
        log.write("# custom atom type overrides applied by apply_custom_atom_types.py\n")
        log.write("# atom_id\tatom_name\told_type\tnew_type\tselector\n")
        for item in log_lines:
            log.write(item + "\n")

    print(f"Applied custom atom types: {changed}")
    print(f"MOL2 written: {output_mol2}")
    print(f"Log written: {log_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
