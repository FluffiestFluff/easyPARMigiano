#!/usr/bin/env python3
"""
extract_minimal_frcmod.py

Create minimal.frcmod containing only parameters that mention custom atom types.

Input custom atom-type file formats:
    7        sD1
    atom_id 7 sD1
    atom_name S7 sD1
    atom_type s1 sD1

The script extracts the final/new atom type from each active non-comment line.

Policy:
    MASS: keep entries whose atom type is custom.
    NONBON: keep entries whose atom type is custom.
    BOND/ANGLE/DIHE/IMPROPER: keep entries where at least one atom type is custom.

Comment policy for kept parameter lines:
    Keep only:
        # atoms=... ep2_source=...
    Drop all other provenance fields/comments in minimal.frcmod.

The full provenance remains in COMPLEX.frcmod and ep2_parameter_provenance.tsv.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set


SECTIONS = {"MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"}


def strip_inline_comment(line: str) -> str:
    return line.split("#", 1)[0].strip()


def get_inline_comment(line: str) -> str:
    if "#" not in line:
        return ""
    return line.split("#", 1)[1].strip()


def read_custom_atom_types(path: Path) -> Set[str]:
    custom_types: Set[str] = set()

    with path.open() as handle:
        for raw in handle:
            line = strip_inline_comment(raw)
            if not line:
                continue

            parts = line.split()
            if len(parts) < 2:
                continue

            if parts[0] in {"atom_id", "atom_name", "atom_type"}:
                if len(parts) >= 3:
                    custom_types.add(parts[2])
            else:
                custom_types.add(parts[1])

    return custom_types


def split_key_atom_types(line: str) -> List[str]:
    body = strip_inline_comment(line)
    if not body:
        return []

    match = re.match(r"^(.+?)(?=\s+[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)", body)
    if match:
        key = match.group(1).strip()
    else:
        fields = body.split()
        key = fields[0].strip() if fields else ""

    if not key:
        return []

    if "-" not in key:
        return [key]

    return [part.strip() for part in re.split(r"\s*-\s*", key) if part.strip()]


def first_token(line: str) -> str:
    body = strip_inline_comment(line)
    parts = body.split()
    return parts[0] if parts else ""


def line_contains_custom_type(line: str, custom_types: Set[str], section: str) -> bool:
    if section in {"BOND", "ANGLE", "DIHE", "IMPROPER"}:
        atom_types = split_key_atom_types(line)
        return any(atom_type in custom_types for atom_type in atom_types)

    if section in {"MASS", "NONBON"}:
        return first_token(line) in custom_types

    return False


def extract_comment_field(comment: str, key: str) -> Optional[str]:
    """
    Extract field from comments like:
        atoms=1,2 ep2_source=HESSIAN_RAW source_file=...
        ep2_source=HESSIAN_UNAVERAGED n=1 source_members="... atoms=1,2 ..."
    """
    if key == "ep2_source":
        m = re.search(r"(?:^|\s)ep2_source=([^\s;#]+)", comment)
        return m.group(1) if m else None

    if key == "atoms":
        # Prefer direct atoms=...
        m = re.search(r"(?:^|\s)atoms=([^\s;#]+)", comment)
        if m:
            return m.group(1)

        # In averaged members/comments, atoms may be embedded in quoted strings.
        m = re.search(r"atoms=([0-9,;:\-]+)", comment)
        return m.group(1) if m else None

    return None


def minimal_comment(line: str) -> str:
    comment = get_inline_comment(line)

    atoms = extract_comment_field(comment, "atoms")
    source = extract_comment_field(comment, "ep2_source")

    fields = []
    if atoms:
        fields.append(f"atoms={atoms}")
    if source:
        fields.append(f"ep2_source={source}")

    if not fields:
        return ""

    return " # " + " ".join(fields)


def clean_parameter_line(line: str) -> str:
    body = strip_inline_comment(line)
    if not body:
        return ""
    return body.rstrip() + minimal_comment(line) + "\n"


def sectionize_frcmod(lines: List[str]) -> Dict[str, List[str]]:
    sections: Dict[str, List[str]] = {section: [] for section in SECTIONS}
    current = None

    for line in lines:
        stripped = line.strip()
        if stripped in SECTIONS:
            current = stripped
            continue

        if current in SECTIONS:
            sections[current].append(line)

    return sections


def filter_section(section: str, lines: Iterable[str], custom_types: Set[str]) -> List[str]:
    kept: List[str] = []

    for line in lines:
        stripped = line.strip()

        if not stripped:
            continue

        if stripped.startswith("#") or stripped.startswith("Remark"):
            continue

        if line_contains_custom_type(line, custom_types, section):
            cleaned = clean_parameter_line(line)
            if cleaned:
                kept.append(cleaned)

    return kept


def write_minimal_frcmod(input_frcmod: Path, output_frcmod: Path, custom_types: Set[str]) -> None:
    lines = input_frcmod.read_text().splitlines(keepends=True)
    sections = sectionize_frcmod(lines)

    with output_frcmod.open("w") as out:
        out.write("Remark minimal frcmod generated by EP2\n")
        out.write("Remark contains only parameters mentioning custom atom types\n")
        out.write("Remark comments reduced to atoms= and ep2_source= only\n")
        out.write("Remark custom_atom_types: " + ",".join(sorted(custom_types)) + "\n\n")

        for section in ["MASS", "BOND", "ANGLE", "DIHE", "IMPROPER", "NONBON"]:
            kept = filter_section(section, sections.get(section, []), custom_types)
            out.write(section + "\n")
            for line in kept:
                out.write(line)
            out.write("\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate minimal.frcmod containing only parameters with custom atom types."
    )
    parser.add_argument("--custom-types", required=True, help="custom_atom_types.dat file.")
    parser.add_argument("--frcmod", default="COMPLEX.frcmod", help="Input frcmod.")
    parser.add_argument("--output", default="minimal.frcmod", help="Output frcmod.")
    args = parser.parse_args()

    custom_path = Path(args.custom_types).expanduser().resolve()
    frcmod_path = Path(args.frcmod).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()

    if not custom_path.exists():
        raise FileNotFoundError(f"Custom atom-type file not found: {custom_path}")

    if not frcmod_path.exists():
        raise FileNotFoundError(f"Input frcmod not found: {frcmod_path}")

    custom_types = read_custom_atom_types(custom_path)
    if not custom_types:
        raise RuntimeError(f"No custom atom types parsed from: {custom_path}")

    write_minimal_frcmod(frcmod_path, output_path, custom_types)

    print(f"Custom atom types: {', '.join(sorted(custom_types))}")
    print(f"Wrote: {output_path}")


if __name__ == "__main__":
    main()
