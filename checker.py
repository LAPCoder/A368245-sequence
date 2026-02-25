#!/usr/bin/env python3
# This script was 100% AI-generated
# Do like me, you shouldn't read it
"""
Verification script for base^n computations with digit sum checking.

FORMULA VERIFIED:
    n-th_root(big_number) == digit_sum_in_base(big_number) - n

USAGE:
    ./verify_hex_sum.py [--base BASE] < input.txt
    
    (./your_c_program 1000 30 > out_hex.pgm >/dev/null) | python checker.py --base 16
    
    Example:
        ./verify_hex_sum.py --base 16 < results.txt
        cat results.txt | ./verify_hex_sum.py --base 16
        ./your_c_program | ./verify_hex_sum.py --base 16

INPUT FORMAT:
    Each line should be: BIG_NUMBER - N
    
    Example:
        4 - 2
        196 - 2
        481890304 - 6

OPTIONS:
    --base BASE    : Set the digit base (default: 16 for hexadecimal)
                     Supported: 2 (binary), 8 (octal), 10 (decimal), 16 (hex)
    --help         : Show this help message

OUTPUT:
    ✓ indicates verification passed
    ✗ indicates verification failed
"""

import sys
import argparse

def digit_sum(n, base=16):
    """Sum the digits of a number in a given base"""
    if base == 16:
        return sum(int(d, 16) for d in hex(n)[2:])
    elif base == 10:
        return sum(int(d) for d in str(n))
    elif base == 8:
        return sum(int(d, 8) for d in oct(n)[2:])
    elif base == 2:
        return sum(int(d) for d in bin(n)[2:])
    else:
        raise ValueError(f"Unsupported base: {base}. Use 2, 8, 10, or 16.")

def integer_root(n, k):
    """Compute the integer k-th root of n using binary search"""
    if n == 0:
        return 0
    if k == 1:
        return n
    
    low, high = 0, n
    while low <= high:
        mid = (low + high) // 2
        mid_pow = mid ** k
        
        if mid_pow == n:
            return mid
        elif mid_pow < n:
            low = mid + 1
        else:
            high = mid - 1
    
    return high  # Return floor of root

def verify_line(big_number, n, base=16):
    """
    Verify: n-th_root(big_number) == digit_sum_in_base(big_number) - n
    
    Args:
        big_number: The large number (result of root^n)
        n: The exponent/root degree
        base: The numerical base for digit sum (default: 16)
    
    Returns:
        True if verification passed, False otherwise
    """
    # Compute the root
    root = integer_root(big_number, n)
    
    # Verify root^n == big_number
    if root ** n != big_number:
        print(f"✗ {root}^{n} ≠ {big_number}")
        return False
    
    # Compute digit sum in the specified base
    digit_sum_value = digit_sum(big_number, base)
    
    # Expected root value
    expected_root = digit_sum_value - n
    
    # Verify the formula
    if root == expected_root:
        base_name = {2: "binary", 8: "octal", 10: "decimal", 16: "hexadecimal"}.get(base, f"base-{base}")
        print(f"✓ {big_number} (n={n}): root={root}, {base_name}_sum={digit_sum_value}, {root} == {digit_sum_value}-{n}")
        return True
    else:
        base_name = {2: "binary", 8: "octal", 10: "decimal", 16: "hexadecimal"}.get(base, f"base-{base}")
        print(f"✗ {big_number} (n={n}): root={root}, {base_name}_sum={digit_sum_value}, {root} ≠ {digit_sum_value}-{n}")
        return False

def main():
    """Read stdin and verify each line"""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Verify: n-th_root(big_number) == digit_sum(big_number) - n",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  ./verify_hex_sum.py --base 16 < results.txt
  cat results.txt | ./verify_hex_sum.py --base 10
  ./your_c_program | ./verify_hex_sum.py --base 16
        """
    )
    parser.add_argument(
        '--base',
        type=int,
        default=16,
        choices=[2, 8, 10, 16],
        help='Numerical base for digit sum (default: 16 for hexadecimal)'
    )
    
    args = parser.parse_args()
    base = args.base
    
    base_names = {2: "binary", 8: "octal", 10: "decimal", 16: "hexadecimal"}
    print(f"Verifying with {base_names[base]} (base-{base}) digit sums...\n")
    
    passed = 0
    failed = 0
    
    # Read and process each line from stdin
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        
        # Parse input: "BIG_NUMBER - N"
        parts = line.split('-')
        if len(parts) != 2:
            print(f"Invalid format: {line}")
            failed += 1
            continue
        
        try:
            big_number = int(parts[0].strip())
            n = int(parts[1].strip())
            
            if verify_line(big_number, n, base):
                passed += 1
            else:
                failed += 1
        except ValueError as e:
            print(f"Parse error: {line} ({e})")
            failed += 1
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"Results: {passed} passed, {failed} failed")
    
    # Exit with error code if any failed
    sys.exit(0 if failed == 0 else 1)

if __name__ == "__main__":
    main()