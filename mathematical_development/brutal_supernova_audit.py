#!/usr/bin/env python3
"""
Brutal Supernova Data Audit
===========================

Apply the same brutal honesty to supernova data that we applied to galactic dynamics.
No mercy - expose all problems, contradictions, and weaknesses.

Key questions:
1. Are the R0 values realistic or fitted?
2. How does UDT compare to LCDM properly?
3. Are we cherry-picking good results?
4. What do the chi2 values really mean?

Author: Charles Rotter
Date: 2025-01-17
"""

import numpy as np
import pandas as pd

class BrutalSupernovaAudit:
    """Brutally honest audit of UDT supernova analysis."""
    
    def __init__(self):
        print("BRUTAL SUPERNOVA DATA AUDIT")
        print("=" * 35)
        print("Applying maximum scientific skepticism")
        print("No mercy for wishful thinking")
        print("=" * 35)
        print()
    
    def analyze_reported_results(self):
        """Analyze the claimed supernova results."""
        
        print("SUPERNOVA RESULTS CLAIMED:")
        print("-" * 30)
        
        # From the files we just read
        print("CSP DR3 (RAW data):")
        print("  R0 = 4,645 Mpc")
        print("  M_B = -18.56")  
        print("  RMS = 1.168 mag")
        print("  chi2/dof = 8,627")
        print()
        
        print("Pantheon+ (RAW mB):")
        print("  R0 = 3,153 Mpc")
        print("  M_B = -18.56")
        print("  RMS = 0.360 mag")
        print("  chi2/dof = 64.1")
        print()
        
        print("IMMEDIATE RED FLAGS:")
        print("1. R0 values differ by 47% between datasets")
        print("2. CSP chi2/dof = 8,627 is CATASTROPHICALLY BAD")
        print("3. Pantheon chi2/dof = 64 is also very poor")
        print("4. RMS values differ by factor of 3.2")
        print()
        
    def assess_chi2_disaster(self):
        """Assess what the chi2 values actually mean."""
        
        print("CHI-SQUARED DISASTER ANALYSIS:")
        print("-" * 35)
        
        print("CSP: chi2/dof = 8,627")
        print("- This means UDT fits CSP data ~8,600x WORSE than random chance")
        print("- This is not a 'reasonable fit' - it's a catastrophic failure")
        print("- For comparison, chi2/dof > 5 is considered very poor")
        print("- chi2/dof = 8,627 means the model is fundamentally wrong")
        print()
        
        print("Pantheon+: chi2/dof = 64")
        print("- Still terrible - means model is ~64x worse than expected")
        print("- Any chi2/dof > 3-5 indicates serious model problems")
        print("- This is not evidence FOR UDT, it's evidence AGAINST it")
        print()
        
        print("WHAT GOOD FITS LOOK LIKE:")
        print("- chi2/dof ~ 0.8-1.2: Excellent fit")
        print("- chi2/dof ~ 1.2-2.0: Acceptable fit")
        print("- chi2/dof ~ 2.0-5.0: Poor fit, model likely wrong")
        print("- chi2/dof > 5.0: Model is fundamentally wrong")
        print()
        
        print("VERDICT: Both UDT supernova fits are TERRIBLE")
        print()
        
    def compare_with_lcdm_properly(self):
        """What would LCDM give for these same datasets?"""
        
        print("PROPER LCDM COMPARISON:")
        print("-" * 25)
        
        print("The documentation claims:")
        print("'UDT beats LCDM on supernova data'")
        print()
        
        print("BUT WHERE IS THE ACTUAL LCDM FIT?")
        print("- No LCDM chi2 values reported")
        print("- No direct comparison on same data")
        print("- No error bars or significance tests")
        print()
        
        print("EXPECTED LCDM PERFORMANCE:")
        print("- Pantheon+ designed for LCDM: chi2/dof ~ 1.0-1.5")
        print("- CSP with LCDM: chi2/dof ~ 1.0-2.0")
        print()
        
        print("UDT REALITY CHECK:")
        print("If UDT gives chi2/dof = 64-8627")
        print("and LCDM gives chi2/dof ~ 1-2")
        print("Then LCDM is ~30-4000x BETTER than UDT!")
        print()
        
        print("CONCLUSION: UDT does NOT beat LCDM on supernovae")
        print("The opposite is almost certainly true")
        print()
        
    def examine_parameter_fitting(self):
        """Are the R0 values genuinely predicted or just fitted?"""
        
        print("PARAMETER FITTING ANALYSIS:")
        print("-" * 30)
        
        print("UDT fits TWO parameters:")
        print("1. R0 (characteristic scale)")
        print("2. M_B (absolute magnitude)")
        print()
        
        print("LCDM also fits TWO parameters:")
        print("1. H0 (Hubble constant)")
        print("2. M_B (absolute magnitude)")
        print()
        
        print("FAIR COMPARISON?")
        print("- Both models have same number of free parameters")
        print("- Both can fit magnitude-redshift relation")
        print("- UDT is NOT more fundamental if it requires fitting")
        print()
        
        print("R0 VALUE PROBLEMS:")
        print("- CSP gives R0 = 4,645 Mpc")
        print("- Pantheon gives R0 = 3,153 Mpc")
        print("- 47% difference between datasets")
        print("- If R0 is fundamental, why does it vary?")
        print()
        
        print("GALACTIC COMPARISON:")
        print("- Galactic R0 was 'fitted' as 38 kpc")
        print("- Cosmological R0 is ~3000-4600 Mpc")  
        print("- Scale ratio: ~100,000:1")
        print("- This is the same scale problem we identified!")
        print()
        
    def assess_data_contamination_claims(self):
        """Assess claims about data contamination."""
        
        print("DATA CONTAMINATION ASSESSMENT:")
        print("-" * 35)
        
        print("CONTAMINATION CLAIMS:")
        print("- 'Previous analysis used contaminated m_b_corr'")
        print("- 'Now using clean raw mB data'")
        print("- 'Eliminated unrealistic chi2/dof = 0.54'")
        print()
        
        print("REALITY CHECK:")
        print("The 'fix' made results WORSE, not better:")
        print("- chi2/dof went from 0.54 to 64.1")
        print("- That's ~120x worse performance")
        print("- A 'clean' dataset should give better fits, not worse")
        print()
        
        print("WHAT THIS REALLY MEANS:")
        print("- The 'contaminated' data actually fit UDT well")
        print("- The 'clean' data reveals UDT doesn't work")
        print("- This is evidence AGAINST UDT, not for it")
        print()
        
        print("DATA CONTAMINATION PARADOX:")
        print("If UDT is correct physics:")
        print("- Clean data should fit better than contaminated data")
        print("- Instead, clean data fits 120x worse")
        print("- This suggests UDT is wrong, not that data was contaminated")
        print()
        
    def examine_cosmological_scale_issue(self):
        """Examine if the cosmological R0 values make sense."""
        
        print("COSMOLOGICAL R0 SCALE ANALYSIS:")
        print("-" * 35)
        
        print("CLAIMED R0 VALUES:")
        print("- CSP: 4,645 Mpc")
        print("- Pantheon: 3,153 Mpc")
        print("- Average: ~3,900 Mpc")
        print()
        
        print("COSMOLOGICAL CONTEXT:")
        print("- Observable universe radius: ~14,000 Mpc")
        print("- Hubble radius c/H0: ~4,300 Mpc")
        print("- UDT R0 ~ Hubble radius")
        print()
        
        print("PHYSICAL INTERPRETATION:")
        print("UDT claims temporal dilation becomes important")
        print("at distances comparable to the Hubble radius")
        print()
        
        print("BUT WAIT - HUGE PROBLEM:")
        print("If R0 ~ Hubble radius, then for galaxies:")
        print("- Galaxy scale: 10-50 kpc")
        print("- r/R0 ~ 0.000003-0.000015")
        print("- tau = R0/(R0+r) ≈ 1 - r/R0 ≈ 0.999997")
        print("- UDT effects are NEGLIGIBLE at galactic scales!")
        print()
        
        print("THIS PROVES GALACTIC UDT MUST BE WRONG:")
        print("With cosmologically determined R0 ~ 4000 Mpc,")
        print("galactic UDT effects are 10^-6 level - unmeasurable!")
        print()
        
    def final_brutal_assessment(self):
        """Final brutal scientific assessment."""
        
        print("FINAL BRUTAL ASSESSMENT:")
        print("=" * 30)
        print()
        
        print("SUPERNOVA DATA REVEALS UDT PROBLEMS:")
        print()
        
        print("1. TERRIBLE FITS TO DATA:")
        print("   - chi2/dof = 64-8627 (should be ~1)")
        print("   - UDT performs catastrophically worse than LCDM")
        print("   - This is evidence AGAINST UDT")
        print()
        
        print("2. INCONSISTENT PARAMETERS:")
        print("   - R0 varies by 47% between datasets")
        print("   - If UDT is fundamental physics, R0 should be universal")
        print("   - Parameter fitting, not genuine physics")
        print()
        
        print("3. SCALE PROBLEM CONFIRMED:")
        print("   - Cosmological R0 ~ 4000 Mpc makes galactic UDT negligible")
        print("   - This mathematically proves galactic applications are wrong")
        print("   - Scale mismatch is fatal to unified theory")
        print()
        
        print("4. DATA CONTAMINATION EXCUSE:")
        print("   - 'Clean' data fits 120x worse than 'contaminated' data")
        print("   - If UDT were correct, clean data should fit better")
        print("   - This is evidence UDT is wrong, not that data was bad")
        print()
        
        print("5. NO PROPER LCDM COMPARISON:")
        print("   - Claims 'UDT beats LCDM' without showing LCDM fits")
        print("   - LCDM almost certainly gives chi2/dof ~ 1")
        print("   - UDT is likely ~50-5000x worse than LCDM")
        print()
        
        print("SCIENTIFIC VERDICT:")
        print("===============")
        print("THE SUPERNOVA DATA DOES NOT SUPPORT UDT")
        print()
        print("Evidence:")
        print("- Catastrophically poor fits (chi2/dof >> 1)")
        print("- Inconsistent parameters between datasets")
        print("- Scale problems that doom galactic applications")
        print("- No evidence UDT beats LCDM (likely much worse)")
        print()
        
        print("The supernova analysis, when examined critically,")
        print("provides evidence AGAINST UDT, not for it.")
        print()
        
        return "supernova_evidence_against_udt"
    
    def run_complete_audit(self):
        """Run complete brutal audit."""
        
        self.analyze_reported_results()
        self.assess_chi2_disaster()
        self.compare_with_lcdm_properly()
        self.examine_parameter_fitting()
        self.assess_data_contamination_claims()
        self.examine_cosmological_scale_issue()
        verdict = self.final_brutal_assessment()
        
        return verdict

def main():
    """Run brutal supernova audit."""
    
    auditor = BrutalSupernovaAudit()
    verdict = auditor.run_complete_audit()
    
    print("AUDIT COMPLETE")
    print("=" * 15)
    print("The supernova data does NOT validate UDT.")
    print("Critical examination reveals fundamental problems")
    print("that were obscured by optimistic interpretation.")
    
    return verdict

if __name__ == "__main__":
    main()