#!/bin/bash
# Script to build a gkeyll module and check for warnings
# Usage: ./make-module.sh <module> <target> [jobs]
# Example: ./make-module.sh gyrokinetic regression 8
# Example: ./make-module.sh core unit 8

set -o pipefail

MODULE="$1"
TARGET="$2"
JOBS="${3:-8}"  # Default to 8 parallel jobs

if [ -z "$MODULE" ] || [ -z "$TARGET" ]; then
    echo "Usage: $0 <module> <target> [jobs]"
    echo "Available modules: core, vlasov, gyrokinetic, moments, pkpm"
    echo "Available targets: unit, regression"
    echo "Example: $0 gyrokinetic regression 8"
    exit 1
fi

# Validate target
if [[ ! "$TARGET" =~ ^(unit|regression)$ ]]; then
    echo "Error: Invalid target '$TARGET'. Must be 'unit' or 'regression'"
    exit 1
fi

BUILD_LOG="build_${MODULE}_${TARGET}.log"
WARNINGS_LOG="${MODULE}_${TARGET}_warnings.log"
ERRORS_LOG="${MODULE}_${TARGET}_errors.log"

echo "========================================"
echo "Building ${MODULE}-${TARGET} with ${JOBS} jobs"
echo "========================================"

# Run make and capture output silently (only save to log)
# Don't exit on error yet - we want to check the log first
MODULE_UPPER=$(echo "$MODULE" | tr '[:lower:]' '[:upper:]')
if ! make -j"${JOBS}" "${MODULE}-${TARGET}" > "$BUILD_LOG" 2>&1; then
    # Build failed - show the errors
    echo ""
    echo "::error::Build failed for ${MODULE} ${TARGET}"
    
    # Try to extract specific error lines
    if grep -E "(error:|undefined reference|fatal error)" "$BUILD_LOG" > "$ERRORS_LOG" 2>/dev/null; then
        echo "=== COMPILER ERRORS IN ${MODULE_UPPER} ==="
        cat "$ERRORS_LOG"
        echo "================================="
    else
        # If no specific errors found, show the last 50 lines of the log
        echo "=== BUILD OUTPUT (last 50 lines) ==="
        tail -50 "$BUILD_LOG"
        echo "================================="
    fi
    
    exit 1
fi

# Build succeeded - check for warnings
if grep -E "warning:" "$BUILD_LOG" > "$WARNINGS_LOG" 2>/dev/null; then
    echo ""
    echo "::warning::Compiler warnings found in ${MODULE} module"
    echo "=== COMPILER WARNINGS IN ${MODULE_UPPER} ==="
    cat "$WARNINGS_LOG"
    echo "================================="
fi

echo "========================================"
echo "${MODULE} ${TARGET} build completed successfully"
echo "========================================"

# Remove only if variables are set and non-empty
# I'm suspicious that we might exceed our CI storage limits if we do not delete these files
for log in "$BUILD_LOG" "$WARNINGS_LOG" "$ERRORS_LOG"; do
    [[ -n "$log" ]] && rm -f "$log"
done