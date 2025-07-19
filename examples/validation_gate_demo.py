#!/usr/bin/env python3
"""
Validation Gate System Demonstration

This demonstrates how the mandatory validation gate prevents analysis
without proper validation procedures.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from udt.data_loaders import requires_validation, validate_before_analysis, ValidationError


# Example 1: Analysis function that REQUIRES validation
@requires_validation('cmb', 'data/cmb_planck/')
def analyze_cmb_power_spectrum(data_directory):
    """
    This function can only run AFTER validation is completed.
    The decorator ensures validation procedures are followed.
    """
    print(f"✓ CMB analysis running on validated data from {data_directory}")
    print("✓ This means:")
    print("  - Data integrity was verified")
    print("  - Artifact correction was validated")
    print("  - Bias testing was completed")
    print("  - Systematic uncertainties were assessed")
    return "CMB analysis completed with validated procedures"


# Example 2: Analysis without validation (will be BLOCKED)
@requires_validation('supernova', 'data/')
def analyze_supernovae(data_directory):
    """This will be blocked until validation is completed."""
    print(f"Supernova analysis on {data_directory}")
    return "Supernova analysis completed"


def demonstrate_validation_gate():
    """Demonstrate how the validation gate system works."""
    
    print("VALIDATION GATE SYSTEM DEMONSTRATION")
    print("=" * 40)
    
    # Attempt 1: Try to run analysis WITHOUT validation (should fail)
    print("\n1. Attempting analysis WITHOUT validation:")
    print("   analyze_cmb_power_spectrum('data/cmb_planck/')")
    
    try:
        result = analyze_cmb_power_spectrum('data/cmb_planck/')
        print(f"   Result: {result}")
    except ValidationError as e:
        print("   ✗ BLOCKED! Validation gate prevented analysis:")
        print(f"   {str(e)[:200]}...")
        print("   → Analysis cannot proceed without validation")
    
    # Step 2: Complete proper validation
    print("\n2. Completing mandatory validation:")
    
    try:
        validation_result = validate_before_analysis('cmb', 'data/cmb_planck/')
        print(f"   ✓ Validation completed: {validation_result['status']}")
        print(f"   ✓ Checks completed: {validation_result['checks_completed']}")
    except Exception as e:
        print(f"   ✗ Validation failed: {e}")
        return
    
    # Attempt 3: Try analysis AFTER validation (should succeed)
    print("\n3. Attempting analysis AFTER validation:")
    print("   analyze_cmb_power_spectrum('data/cmb_planck/')")
    
    try:
        result = analyze_cmb_power_spectrum('data/cmb_planck/')
        print(f"   ✓ SUCCESS! {result}")
    except Exception as e:
        print(f"   ✗ Unexpected error: {e}")
    
    # Attempt 4: Try different dataset without validation (should fail)
    print("\n4. Attempting different dataset WITHOUT validation:")
    print("   analyze_supernovae('data/')")
    
    try:
        result = analyze_supernovae('data/')
        print(f"   Result: {result}")
    except ValidationError as e:
        print("   ✗ BLOCKED! Each dataset requires separate validation")
        print("   → Must run validate_before_analysis('supernova', 'data/') first")
    
    print("\n" + "=" * 50)
    print("VALIDATION GATE SYSTEM SUMMARY:")
    print("✓ Prevents analysis without proper validation")
    print("✓ Enforces scientific rigor requirements")
    print("✓ Blocks bypass attempts with clear error messages")
    print("✓ Ensures artifact correction and bias testing")
    print("✓ Maintains validation logs for audit trail")
    print("\nThis system prevents the mistakes I made earlier!")


if __name__ == "__main__":
    demonstrate_validation_gate()