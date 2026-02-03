#!/bin/bash
# -----------------------------------------------------------------------------
# Post-process CUDA files after uncrustify formatting
# -----------------------------------------------------------------------------
# Fixes issues caused by uncrustify not understanding CUDA syntax:
# 1. Removes spaces in/around kernel launch operators: << < -> <<<, >> > -> >>>
# 2. Fixes indentation of arguments after kernel launch calls
#
# Usage:
#   .github/cuda-after-uncrustify.sh <file.cu> [file2.cu ...]
#   .github/cuda-after-uncrustify.sh core/zero/*.cu
# -----------------------------------------------------------------------------

if [ $# -eq 0 ]; then
  echo "Usage: $0 <file.cu> [file2.cu ...]"
  echo "Post-processes CUDA files after uncrustify formatting."
  exit 1
fi

for file in "$@"; do
  if [ ! -f "$file" ]; then
    echo "Warning: File not found: $file"
    continue
  fi

  echo "Processing: $file"

  # Fix 1: Remove spaces in and around <<< and >>> operators
  # Handles: << < -> <<<, >> > -> >>>, and removes surrounding spaces
  sed -i 's/<< </<<</g; s/>> >/>>>/g' "$file"
  sed -i 's/ <<<\s*/<<</' "$file"
  sed -i 's/\s*>>> />>>/' "$file"

  # Fix 2: Fix indentation of continuation lines after kernel launches
  # Find lines with <<<...>>>( and fix indentation of following lines
  # until we hit the closing );
  perl -i -pe '
    BEGIN { $in_kernel_call = 0; $base_indent = ""; }
    
    # Detect kernel launch line ending with (
    if (/^(\s*)(\S+)\s*<<<.*>>>\s*\(\s*$/) {
      $in_kernel_call = 1;
      $base_indent = $1 . "  ";  # 2 spaces from function start
      next;
    }
    
    # If we are in a kernel call, fix indentation
    if ($in_kernel_call) {
      # Check if this line ends the call
      if (/\);\s*$/) {
        s/^(\s*)/$base_indent/;
        $in_kernel_call = 0;
      } else {
        # Fix indentation of continuation lines
        s/^(\s*)/$base_indent/;
      }
    }
  ' "$file"

	# Fix 3: Function return type pointer style:
	# When a line ends with " *" and the next line is a function name,
	# change " *" to "*" (remove the space before the asterisk)
  perl -i -0pe '
    # Match: (type) space asterisk newline (function_name)
    # Replace with: (type) asterisk newline (function_name)
    s/(\w) \*\n(\w+\s*\()/$1*\n$2/g;
  ' "$file"

done

echo "Done."
