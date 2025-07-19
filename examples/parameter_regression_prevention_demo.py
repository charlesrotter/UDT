#!/usr/bin/env python3
"""
Parameter Regression Prevention System Demonstration

Shows how the parameter registry prevents the specific error I made:
using galactic-scale R0 for CMB analysis.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from udt.diagnostics.parameter_registry import (
    get_validated_parameters, 
    validate_analysis_parameters,
    parameter_registry
)
from udt.diagnostics.mandatory_validation_gate import requires_validation, ValidationError


def demonstrate_parameter_prevention():
    """Demonstrate how the system prevents parameter regression."""
    
    print("PARAMETER REGRESSION PREVENTION DEMONSTRATION")
    print("=" * 50)
    print("Preventing the specific error I made: wrong R0 scale for CMB")
    print()
    
    # Show the correct parameters for different analyses
    print("1. VALIDATED PARAMETER REGISTRY:")
    print("-" * 30)
    
    for analysis_type in ['galaxy', 'supernova', 'cmb']:
        params = get_validated_parameters(analysis_type)
        print(f"{analysis_type.upper():>10}: R0 = {params['R0_mpc']:>8.1f} Mpc ({params['description']})")
    
    print()
    
    # Demonstrate parameter validation catching errors
    print("2. PARAMETER VALIDATION CATCHING ERRORS:")
    print("-" * 40)
    
    # Case A: Correct CMB parameters (should pass)
    print("Case A: Using CORRECT CMB parameters")
    correct_cmb_params = {'R0_mpc': 13041.1}
    validation = validate_analysis_parameters('cmb', correct_cmb_params)
    print(f"   Status: {validation['status']}")
    if validation['warnings']:
        print(f"   Warnings: {validation['warnings']}")
    if validation['errors']:
        print(f"   Errors: {validation['errors']}")
    print("   → Analysis would proceed normally")
    print()
    
    # Case B: Wrong CMB parameters (my error - should catch)
    print("Case B: Using WRONG CMB parameters (my error)")
    wrong_cmb_params = {'R0_mpc': 3582}  # Galactic scale used for CMB!
    validation = validate_analysis_parameters('cmb', wrong_cmb_params)
    print(f"   Status: {validation['status']}")
    if validation['warnings']:
        for warning in validation['warnings']:
            print(f"   ⚠ Warning: {warning}")
    if validation['errors']:
        for error in validation['errors']:
            print(f"   ✗ ERROR: {error}")
    if validation['recommendations']:
        for rec in validation['recommendations']:
            print(f"   → {rec}")
    print()
    
    # Case C: Galaxy parameters for galaxy analysis (should pass)
    print("Case C: Using CORRECT galaxy parameters")
    correct_galaxy_params = {'R0_mpc': 0.038}
    validation = validate_analysis_parameters('galaxy', correct_galaxy_params)
    print(f"   Status: {validation['status']}")
    print("   → Analysis would proceed normally")
    print()
    
    # Show automatic parameter selection
    print("3. AUTOMATIC PARAMETER SELECTION:")
    print("-" * 35)
    print("Instead of manually choosing parameters, use registry:")
    print()
    
    for analysis_type in ['galaxy', 'supernova', 'cmb']:
        params = get_validated_parameters(analysis_type)
        R0 = params['R0_mpc']
        scale = params['scale']
        print(f"# For {analysis_type} analysis:")
        print(f"params = get_validated_parameters('{analysis_type}')")
        print(f"R0_{scale} = params['R0_mpc']  # {R0} Mpc")
        print()
    
    # Show validation history
    print("4. VALIDATION SUMMARY:")
    print("-" * 20)
    summary = parameter_registry.get_validation_summary()
    
    print("Scale-specific parameters:")
    for scale, info in summary['scales'].items():
        print(f"  {scale.upper():>13}: R0 = {info['R0_mpc']:>8.1f} Mpc")
        if info['best_chi2_dof']:
            print(f"                Best χ²/dof = {info['best_chi2_dof']:.2f}")
    
    print()
    print("5. INTEGRATION WITH VALIDATION GATE:")
    print("-" * 38)
    print("The validation gate now automatically:")
    print("✓ Loads correct parameters from registry")
    print("✓ Validates provided parameters against registry")  
    print("✓ Warns about parameter deviations")
    print("✓ Blocks analysis with wrong scale parameters")
    print("✓ Tracks validation history and best results")
    
    print()
    print("PREVENTION SUMMARY:")
    print("=" * 20)
    print("✗ BEFORE: Could use any R0 value, leading to wrong scale")
    print("✓ NOW: Registry enforces correct scale-specific parameters")
    print("✓ NOW: Automatic validation catches scale mismatches")
    print("✓ NOW: Clear recommendations for correct parameters")
    print("✓ NOW: Integrated with mandatory validation gate")
    
    print()
    print("This system would have PREVENTED my CMB regression error!")


# Example of protected analysis function
@requires_validation('cmb', 'data/cmb_planck/')
def protected_cmb_analysis(data_directory):
    """
    CMB analysis that's protected by validation gate and parameter registry.
    
    This function will:
    1. Be blocked until validation is completed
    2. Automatically get correct CMB parameters from registry
    3. Validate any provided parameters against registry
    4. Warn about parameter deviations
    """
    # Get validated parameters automatically
    params = get_validated_parameters('cmb')
    R0_cmb = params['R0_mpc']  # Automatically gets 13041.1 Mpc
    
    print(f"✓ Using validated CMB parameters: R0 = {R0_cmb} Mpc")
    print("✓ Analysis proceeding with correct scale parameters")
    
    return f"CMB analysis completed with R0 = {R0_cmb} Mpc"


if __name__ == "__main__":
    demonstrate_parameter_prevention()
    
    print("\n" + "="*60)
    print("TESTING PROTECTED ANALYSIS FUNCTION:")
    print("="*60)
    
    try:
        # This will fail due to validation gate
        result = protected_cmb_analysis('data/cmb_planck/')
        print(f"Result: {result}")
    except ValidationError:
        print("✗ Analysis blocked - validation required first")
        print("→ Must run validate_before_analysis('cmb', 'data/cmb_planck/') first")