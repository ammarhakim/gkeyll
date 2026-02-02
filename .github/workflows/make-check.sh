#!/bin/bash
# Script to run tests for a gkeyll module and analyze results
# Usage: ./make-check.sh <module> [jobs]
# Example: ./make-check.sh core 8

set -o pipefail

MODULE="$1"
JOBS="${2:-8}"  # Default to 8 parallel jobs

if [ -z "$MODULE" ]; then
    echo "Usage: $0 <module> [jobs]"
    echo "Available modules: core, vlasov, gyrokinetic, moments, pkpm"
    exit 1
fi

TEST_LOG="${MODULE}_test.log"

echo "========================================"
echo "Running ${MODULE}-check with ${JOBS} jobs"
echo "========================================"

# Run make and capture output silently
if ! make -j"${JOBS}" "${MODULE}-check" > "$TEST_LOG" 2>&1; then
    echo ""
    echo "::error::Test execution failed for ${MODULE} module"
    # Show the last 50 lines if make itself failed
    echo "=== TEST OUTPUT (last 50 lines) ==="
    tail -50 "$TEST_LOG"
    echo "================================="
    exit 1
fi

# Count test results
PASSED=$(grep -c "\[ OK \]" "$TEST_LOG" 2>/dev/null || echo "0")
FAILED=$(grep -c "\[ FAILED \]" "$TEST_LOG" 2>/dev/null || echo "0")

# Ensure we have single integer values (remove any whitespace/newlines)
PASSED=$(echo "$PASSED" | tr -d '\n' | head -1)
FAILED=$(echo "$FAILED" | tr -d '\n' | head -1)

# Check for segfaults or crashes
CRASHED=0
if grep -qiE "(segmentation fault|segfault|core dumped|signal|abort)" "$TEST_LOG"; then
    echo "::error::${MODULE} tests crashed with segmentation fault or signal"
    CRASHED=1
fi

# Extract failed test names if any
if [ "$FAILED" -gt 0 ]; then
    echo ""
    echo "::error::Failed tests detected"
    echo "=== FAILED TESTS IN ${MODULE^^} ==="
    
    # Extract test executable paths and failed test names
    # Look for lines like: ../build/core/unit/ctest_array
    # followed by lines with [ FAILED ]
    awk '
        /^\.\.\/build\// { 
            # Extract the test path (e.g., ../build/core/unit/ctest_array)
            gsub(/^\.\.\/build\//, "")  # Remove ../build/ prefix
            current_test = $0
        }
        /\[ FAILED \]/ {
            if (current_test != "") {
                print current_test
                current_test = ""
            }
            print $0
        }
    ' "$TEST_LOG"
    
    # Also show the FAILED summary line
    grep "^FAILED:" "$TEST_LOG" || true
    echo "================================="
fi

# Check for unexpected output (lines that aren't test results, echoes, or known patterns)
# This catches extra print statements
UNEXPECTED=$(grep -vE "(^\.\./build/|Test |SUCCESS:|FAILED:|warning:|^ *$|make\[|Entering|Leaving|\[ OK \]|\[ FAILED \]|^cd |^echo )" "$TEST_LOG" | head -20 || true)
if [ -n "$UNEXPECTED" ]; then
    echo ""
    echo "::warning::Unexpected output detected (possible debug prints):"
    echo "$UNEXPECTED"
fi

echo ""
echo "========================================"
echo "${MODULE^^} TEST SUMMARY"
echo "========================================"
echo "Passed:  $PASSED"
echo "Failed:  $FAILED"
echo "========================================"

# Fail if any tests failed or crashed
if [ "$FAILED" -gt 0 ] || [ "$CRASHED" -eq 1 ]; then
    echo "::error::${MODULE} tests failed: $FAILED failures, $CRASHED crashes"
    exit 1
fi

# Check for SUCCESS message from test framework
if ! grep -q "SUCCESS: All unit tests have passed" "$TEST_LOG"; then
    echo "::error::Not all ${MODULE} test suites completed successfully"
    exit 1
fi

echo "========================================"
echo "${MODULE} tests completed successfully"
echo "========================================"
