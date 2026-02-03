#!/bin/bash
# -----------------------------------------------------------------------------
# Uncrustify wrapper for VS Code
# -----------------------------------------------------------------------------
# This script wraps uncrustify to add post-processing for CUDA files.
# Set this as the "Uncrustify › Executable Path" in VS Code settings.
#
# It passes all arguments to uncrustify, then runs cuda-after-uncrustify.sh
# if the file being formatted is a .cu file.
#
# Works with:
#   - Format Document
#   - Format Selection (uncrustify handles this via --frag flag)
#
# How to set up in VS Code:
# 
# Navigate to .vscode settings in this repository and open .vscode/settings.json.
# Then, set the "uncrustify.executablePath" setting to point, e.g.
# 
#   "uncrustify.configPath.linux": "${workspaceFolder}/.github/uncrustify.cfg",
#   "uncrustify.executablePath.linux": "/home/maxwell-rosen/gkeyll/.github/uncrustify-wrapper.sh",
# To enable formatting of .h and .cu files, add these lines to .vscode/settings.json:
#   "files.associations": {
#      "*.h": "cpp",
#      "*.cu": "cpp",
# -----------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CUDA_POST_SCRIPT="$SCRIPT_DIR/cuda-after-uncrustify.sh"
C_H_POST_SCRIPT="$SCRIPT_DIR/c-h-after-uncrustify.sh"

# Find the actual uncrustify binary
UNCRUSTIFY_BIN=$(which uncrustify)
if [ -z "$UNCRUSTIFY_BIN" ]; then
  echo "Error: uncrustify not found in PATH" >&2
  exit 1
fi

# Run uncrustify with all passed arguments
"$UNCRUSTIFY_BIN" "$@"
UNCRUSTIFY_EXIT=$?

# If uncrustify failed, exit with its error code
if [ $UNCRUSTIFY_EXIT -ne 0 ]; then
  exit $UNCRUSTIFY_EXIT
fi

# Check if any argument is a .cu, .c, or .h file and if --replace or -o was used
# (meaning the file was modified in place or output was written)
CUDA_FILE=""
C_H_FILE=""
HAS_REPLACE=0

for arg in "$@"; do
  case "$arg" in
    --replace|--no-backup)
      HAS_REPLACE=1
      ;;
    *.cu)
      CUDA_FILE="$arg"
      ;;
    *.c|*.h)
      C_H_FILE="$arg"
      ;;
  esac
done

# Postprocessing scripts
# .c and .h files
if [ -n "$C_H_FILE" ] && [ $HAS_REPLACE -eq 1 ] && [ -x "$C_H_POST_SCRIPT" ]; then
  "$C_H_POST_SCRIPT" "$C_H_FILE" >/dev/null 2>&1
fi

# CUDA files
if [ -n "$CUDA_FILE" ] && [ $HAS_REPLACE -eq 1 ] && [ -x "$CUDA_POST_SCRIPT" ]; then
  "$CUDA_POST_SCRIPT" "$CUDA_FILE" >/dev/null 2>&1
fi

exit 0
