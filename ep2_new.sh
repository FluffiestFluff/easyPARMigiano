#!/usr/bin/env bash
set -euo pipefail

###################################################################################################################
# ep2_new.sh — lightweight EP2 launcher for easyPARM-derived workflow
#
# Positional interface:
#   ./ep2_new.sh
#       Interactive mode. Answers are typed manually.
#
#   ./ep2_new.sh choices.txt
#       Non-interactive mode. Answers are read line-by-line from choices.txt.
#
#   ./ep2_new.sh choices.txt custom_atom_types.dat
#       Non-interactive mode plus custom atom-type override file.
#
# Explicit interface:
#   ./ep2_new.sh --choices choices.txt
#   ./ep2_new.sh --choices choices.txt --custom-atom-types custom_atom_types.dat
#   ./ep2_new.sh --run-dir /path/to/run --choices choices.txt --custom-atom-types custom_atom_types.dat
#
# Notes:
#   - The choices file is passed directly to 01_ep2.sh through stdin.
#   - Therefore the file must contain answers in exactly the same order as the prompts.
#   - Blank lines are allowed only when a prompt should receive an empty answer.
#   - Comments are not stripped by this launcher; do not add comments inside choices.txt unless 01_ep2.sh expects them.
#   - The custom atom-type file is exported as EP2_CUSTOM_ATOM_TYPES, then used by 01_ep2.sh.
###################################################################################################################

SCRIPT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -d "$SCRIPT_ROOT/scripts" ]; then
    SCRIPT_DIR="$SCRIPT_ROOT/scripts"
else
    SCRIPT_DIR="$SCRIPT_ROOT"
fi

EP2_DRIVER="$SCRIPT_DIR/01_ep2.sh"

if [ ! -f "$EP2_DRIVER" ]; then
    echo "Error: 01_ep2.sh not found in: $SCRIPT_DIR" >&2
    exit 1
fi

usage() {
    cat <<USAGE
Usage:
  $(basename "$0")
      Run EP2 interactively.

  $(basename "$0") choices.txt
      Run EP2 non-interactively using answers from choices.txt.

  $(basename "$0") choices.txt custom_atom_types.dat
      Run EP2 non-interactively and apply custom atom-type overrides.

  $(basename "$0") --choices choices.txt
      Same as above, explicit choices-file form.

  $(basename "$0") --choices choices.txt --custom-atom-types custom_atom_types.dat
      Explicit form with custom atom-type overrides.

  $(basename "$0") --run-dir /path/to/run choices.txt custom_atom_types.dat
      Use /path/to/run as the working directory.

Environment flags are still supported, for example:
  EP2_KEEP_INTERMEDIATES=0 EP2_APPLY_EASYPARM_SCALING=1 $(basename "$0") choices.txt custom_atom_types.dat
USAGE
}

RUN_DIR="$PWD"
CHOICES_FILE=""
CUSTOM_ATOM_TYPES_FILE=""

while [ "$#" -gt 0 ]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        --run-dir)
            if [ "$#" -lt 2 ]; then
                echo "Error: --run-dir requires a path." >&2
                exit 1
            fi
            RUN_DIR="$2"
            shift 2
            ;;
        --choices|--input|--answers)
            if [ "$#" -lt 2 ]; then
                echo "Error: $1 requires a choices file." >&2
                exit 1
            fi
            CHOICES_FILE="$2"
            shift 2
            ;;
        --custom-atom-types|--atom-types|--custom-types)
            if [ "$#" -lt 2 ]; then
                echo "Error: $1 requires a custom atom-type file." >&2
                exit 1
            fi
            CUSTOM_ATOM_TYPES_FILE="$2"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "Error: unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
        *)
            if [ -z "$CHOICES_FILE" ]; then
                CHOICES_FILE="$1"
                shift
            elif [ -z "$CUSTOM_ATOM_TYPES_FILE" ]; then
                CUSTOM_ATOM_TYPES_FILE="$1"
                shift
            else
                echo "Error: unexpected extra positional argument: $1" >&2
                usage >&2
                exit 1
            fi
            ;;
    esac
done

RUN_DIR="$(cd "$RUN_DIR" && pwd)"

if [ -n "$CHOICES_FILE" ]; then
    if [ ! -f "$CHOICES_FILE" ]; then
        echo "Error: choices file not found: $CHOICES_FILE" >&2
        exit 1
    fi
    CHOICES_FILE="$(cd "$(dirname "$CHOICES_FILE")" && pwd)/$(basename "$CHOICES_FILE")"
fi

if [ -n "$CUSTOM_ATOM_TYPES_FILE" ]; then
    if [ ! -f "$CUSTOM_ATOM_TYPES_FILE" ]; then
        echo "Error: custom atom-type file not found: $CUSTOM_ATOM_TYPES_FILE" >&2
        exit 1
    fi
    CUSTOM_ATOM_TYPES_FILE="$(cd "$(dirname "$CUSTOM_ATOM_TYPES_FILE")" && pwd)/$(basename "$CUSTOM_ATOM_TYPES_FILE")"
    export EP2_CUSTOM_ATOM_TYPES="$CUSTOM_ATOM_TYPES_FILE"
fi

echo " "
echo "════════════════════════════════════════════════════════════════════════"
echo " easyPARM ep2 launcher"
echo "════════════════════════════════════════════════════════════════════════"
echo "Run directory:       $RUN_DIR"
echo "Script directory:    $SCRIPT_DIR"

if [ -n "$CUSTOM_ATOM_TYPES_FILE" ]; then
    echo "Custom atom types:   $CUSTOM_ATOM_TYPES_FILE"
else
    echo "Custom atom types:   none supplied as arg2"
fi

if [ -n "$CHOICES_FILE" ]; then
    echo "Input mode:          non-interactive"
    echo "Choices file:        $CHOICES_FILE"
    echo " "
    bash "$EP2_DRIVER" "$RUN_DIR" < "$CHOICES_FILE"
else
    echo "Input mode:          interactive"
    echo " "
    bash "$EP2_DRIVER" "$RUN_DIR"
fi
