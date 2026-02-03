#!/bin/bash
# -----------------------------------------------------------------------------
# Post-process C/H files after uncrustify formatting
# -----------------------------------------------------------------------------
# Fixes: Function return type pointer style
#   FROM: struct foo *
#         function_name(...)
#   TO:   struct foo*
#         function_name(...)
#
# Usage:
#   .github/c-after-uncrustify.sh <file.c> [file2.h ...]
#   .github/c-after-uncrustify.sh core/zero/*.h
# -----------------------------------------------------------------------------

if [ $# -eq 0 ]; then
  echo "Usage: $0 <file.c> [file2.h ...]"
  echo "Post-processes C/H files after uncrustify formatting."
  exit 1
fi

for file in "$@"; do
  if [ ! -f "$file" ]; then
    echo "Warning: File not found: $file"
    continue
  fi

  echo "Processing: $file"

  # Fix function return type pointer style:
  # When a line ends with " *" and the next line is a function name,
  # change " *" to "*" (remove the space before the asterisk)
  # 
  # Pattern: "word *\n" followed by a line starting with a function name
  # This handles:
  #   struct gkyl_dg_array_mask *
  #   gkyl_dg_array_mask_acquire(...)
  # Becomes:
  #   struct gkyl_dg_array_mask*
  #   gkyl_dg_array_mask_acquire(...)

  perl -i -0pe '
    # Match: (type) space asterisk newline (function_name)
    # Replace with: (type) asterisk newline (function_name)
    s/(\w) \*\n(\w+\s*\()/$1*\n$2/g;
  ' "$file"

done

echo "Done."
