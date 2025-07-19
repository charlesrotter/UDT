#!/usr/bin/env python3
"""
Test Quantum Validation Regression Fix
======================================

Test if using ParameterRegistry fixes the LIGO and muon g-2 regression.
"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent))
from src.udt.diagnostics.parameter_registry import ParameterRegistry

def test_parameter_compatibility():
    """Test if quantum validation parameters are compatible with registry."""
    
    print("QUANTUM VALIDATION PARAMETER COMPATIBILITY TEST")
    print("=" * 55)
    
    registry = ParameterRegistry()
    
    # Current hardcoded values in quantum scripts
    current_quantum_params = {
        'ligo_R0_cosmo': 3582.0,  # From ligo_gravitational_wave_udt_analysis.py
        'ligo_R0_gw': 3582.0,    # From ligo_gravitational_wave_udt_analysis.py  
        'muon_R0_cosmo': 3582.0, # From muon_g2_udt_analysis.py
        'muon_R0_lab': 1e-10,    # From muon_g2_udt_analysis.py
    }
    
    # Registry validated values
    try:
        cosmo_params = registry.get_parameters_for_analysis('supernova')
        cmb_params = registry.get_parameters_for_analysis('cmb')
        sparc_params = registry.get_parameters_for_analysis('sparc')
        
        registry_params = {
            'cosmo_R0': cosmo_params['R0_mpc'],      # 3000.0 Mpc
            'cmb_R0': cmb_params['R0_mpc'],          # 13041.1 Mpc  
            'galactic_R0': sparc_params['R0_mpc'],   # 0.038 Mpc
        }
        
        print("Current Quantum Script Parameters:")
        for key, value in current_quantum_params.items():
            print(f"  {key}: {value}")
            
        print("\nValidated Registry Parameters:")
        for key, value in registry_params.items():
            print(f"  {key}: {value}")
            
        # Check compatibility
        print("\nCompatibility Analysis:")
        
        # LIGO should probably use cosmological scale
        ligo_diff = abs(current_quantum_params['ligo_R0_cosmo'] - registry_params['cosmo_R0'])
        ligo_percent_diff = (ligo_diff / registry_params['cosmo_R0']) * 100
        
        if ligo_percent_diff < 20:  # Within 20%
            print(f"  + LIGO R0: {ligo_percent_diff:.1f}% difference - ACCEPTABLE")
        else:
            print(f"  - LIGO R0: {ligo_percent_diff:.1f}% difference - POTENTIAL REGRESSION")
            
        # Muon g-2 lab scale is very different - might be intentional
        print(f"  ? Muon lab R0: {current_quantum_params['muon_R0_lab']:.2e} Mpc - quantum scale")
        print(f"    (vs galactic {registry_params['galactic_R0']:.3f} Mpc - very different scales)")
        
        return ligo_percent_diff < 20
        
    except Exception as e:
        print(f"Registry access failed: {e}")
        return False

def suggest_fixes():
    """Suggest fixes for quantum validation regression."""
    
    print("\n" + "=" * 55)
    print("REGRESSION FIX RECOMMENDATIONS")
    print("=" * 55)
    
    print("\n1. IMMEDIATE FIXES:")
    print("   - Update quantum scripts to use ParameterRegistry")
    print("   - Replace hardcoded R0 values with registry lookups")
    print("   - Ensure scale-appropriate parameter selection")
    
    print("\n2. PARAMETER MAPPING:")
    print("   - LIGO: Use cosmological scale (3000 Mpc) for GW propagation")
    print("   - Muon g-2: Use quantum/lab scale (may need specific calibration)")
    print("   - Both: Import ParameterRegistry and ValidationGate")
    
    print("\n3. VALIDATION APPROACH:")
    print("   - Re-run with registry parameters")
    print("   - Compare results to previous VALIDATED benchmarks")
    print("   - Document any legitimate parameter differences")
    
    print("\n4. RISK ASSESSMENT:")
    print("   - HIGH: Quantum domain validation is critical for UDT credibility")
    print("   - MEDIUM: Parameter inconsistency may have cascading effects")
    print("   - LOW: Fixes should be straightforward registry integration")

if __name__ == "__main__":
    # Test compatibility
    compatible = test_parameter_compatibility()
    
    # Provide recommendations
    suggest_fixes()
    
    # Final assessment
    print("\n" + "=" * 55)
    print("QUANTUM REGRESSION ASSESSMENT")
    print("=" * 55)
    
    if compatible:
        print("MINOR REGRESSION: Parameter differences exist but may be acceptable")
        print("Recommend registry integration for consistency and validation")
    else:
        print("MAJOR REGRESSION: Significant parameter inconsistencies detected")
        print("URGENT: Update quantum scripts to use validated parameters")
        
    print("\nNext steps:")
    print("1. Integrate ParameterRegistry into quantum validation scripts")
    print("2. Test if this resolves the LIGO/muon g-2 failures")
    print("3. Compare restored results to previous VALIDATED benchmarks")