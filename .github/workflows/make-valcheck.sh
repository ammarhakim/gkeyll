#!/bin/bash
# Script to run valgrind checks for a gkeyll module and analyze results
# Usage: ./make-valcheck.sh <module> [jobs]
# Example: ./make-valcheck.sh core 1

set -o pipefail

MODULE="$1"
JOBS="${2:-1}"  # Default to 1 parallel job

if [ -z "$MODULE" ]; then
    echo "Usage: $0 <module> [jobs]"
    echo "Available modules: core, vlasov, gyrokinetic, moments, pkpm"
    exit 1
fi

VALCHECK_LOG="${MODULE}_valcheck.log"

echo "========================================"
echo "Running ${MODULE}-valcheck with ${JOBS} jobs"
echo "========================================"

# Run make valcheck and capture output
if ! make -j"${JOBS}" "${MODULE}-valcheck" > "$VALCHECK_LOG" 2>&1; then
    echo ""
    echo "::error::Valgrind check execution failed for ${MODULE} module"
    # Show the last 50 lines if make itself failed
    echo "=== VALCHECK OUTPUT (last 50 lines) ==="
    tail -50 "$VALCHECK_LOG"
    echo "================================="
    exit 1
fi

# Convert module name to uppercase for display (compatible with older bash)
MODULE_UPPER=$(echo "$MODULE" | tr '[:lower:]' '[:upper:]')

# Count valgrind results
CLEAN=$(grep -c "is valgrind clean" "$VALCHECK_LOG" 2>/dev/null || echo "0")
LEAKY=$(grep -c "has error issues or memory leaks" "$VALCHECK_LOG" 2>/dev/null || echo "0")

# Ensure we have single integer values
CLEAN=$(echo "$CLEAN" | tr -d '\n' | head -1)
LEAKY=$(echo "$LEAKY" | tr -d '\n' | head -1)

# Extract tests with issues
if [ "$LEAKY" -gt 0 ]; then
    echo ""
    echo "::error::Memory leaks or errors detected in ${MODULE} module"
    echo "=== VALGRIND ISSUES IN ${MODULE_UPPER} ==="
    
    # Show which tests have issues
    grep "has error issues or memory leaks" "$VALCHECK_LOG"
    
    echo ""
    echo "=== VALGRIND ERROR DETAILS ==="
    
    # For each leaky test, try to show the valgrind error file contents
    grep "has error issues or memory leaks" "$VALCHECK_LOG" | while read -r line; do
        # Extract the test path from the line (e.g., ../build/core/unit/ctest_array)
        TEST_PATH=$(echo "$line" | grep -oE '\.\./build/[^ ]+')
        if [ -n "$TEST_PATH" ]; then
            # Convert ../build to ./build for running from root directory
            VAL_ERR_FILE=$(echo "${TEST_PATH}_val_err" | sed 's|^\.\./|./|')
            TEST_NAME=$(basename "$TEST_PATH")
            
            echo ""
            echo "--- ${TEST_NAME} ---"
            if [ -f "$VAL_ERR_FILE" ]; then
                # Show the error summary from valgrind output
                cat "$VAL_ERR_FILE"
            else
                echo "  (valgrind error file not found: $VAL_ERR_FILE)"
            fi
        fi
    done
    
    echo "================================="
fi

echo ""
echo "========================================"
echo "${MODULE_UPPER} VALGRIND SUMMARY"
echo "========================================"
echo "Clean:   $CLEAN"
echo "Leaky:   $LEAKY"
echo "========================================"

# Fail if any tests have memory issues
if [ "$LEAKY" -gt 0 ]; then
    echo "::error::${MODULE} valgrind check failed: $LEAKY tests with memory issues"
    exit 1
fi

echo "========================================"
echo "${MODULE} valgrind checks completed successfully"
echo "========================================"
