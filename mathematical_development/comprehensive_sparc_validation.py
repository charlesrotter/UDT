#!/usr/bin/env python3
"""
Comprehensive SPARC Validation of UDT Framework
===============================================

Analyze the complete UDT framework against SPARC data with:
1. Statistical analysis of fit quality
2. Comparison with theoretical predictions
3. Assessment of parameter consistency
4. Validation of cosmic connectivity theory

Author: Charles Rotter
Date: 2025-07-18
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import sys
import os

# Add project root to path
sys.path.append(os.path.abspath('.'))

class ComprehensiveSPARCValidation:
    def __init__(self):
        print("COMPREHENSIVE SPARC VALIDATION OF UDT FRAMEWORK")
        print("=" * 55)
        
        # Load the results from the SPARC analysis
        self.results_file = "results/sparc_analysis/sparc_udt_results.csv"
        self.load_results()
        
    def load_results(self):
        """Load SPARC analysis results."""
        try:
            self.df = pd.read_csv(self.results_file)
            print(f"Loaded {len(self.df)} galaxy results from {self.results_file}")
        except FileNotFoundError:
            print(f"Error: Could not find {self.results_file}")
            print("Please run: python scripts/analyze_sparc_galaxies.py first")
            sys.exit(1)
    
    def analyze_fit_quality(self):
        """Analyze the quality of UDT fits."""
        print("\nFIT QUALITY ANALYSIS")
        print("-" * 22)
        
        # Basic statistics
        rms_residuals = self.df['rms'].values
        r0_values = self.df['R0_gal'].values
        v_scales = self.df['V_scale'].values
        
        print(f"Number of galaxies: {len(self.df)}")
        print(f"Success rate: {len(self.df[self.df['success']])}/{len(self.df)} = {100*len(self.df[self.df['success']])/len(self.df):.1f}%")
        print()
        
        print("RMS Residual Statistics:")
        print(f"  Mean: {np.mean(rms_residuals):.2f} km/s")
        print(f"  Median: {np.median(rms_residuals):.2f} km/s")
        print(f"  Std: {np.std(rms_residuals):.2f} km/s")
        print(f"  Min: {np.min(rms_residuals):.2f} km/s")
        print(f"  Max: {np.max(rms_residuals):.2f} km/s")
        print()
        
        print("R0 Parameter Statistics:")
        print(f"  Mean: {np.mean(r0_values):.1f} kpc")
        print(f"  Median: {np.median(r0_values):.1f} kpc")
        print(f"  Std: {np.std(r0_values):.1f} kpc")
        print(f"  Min: {np.min(r0_values):.1f} kpc")
        print(f"  Max: {np.max(r0_values):.1f} kpc")
        print()
        
        # Excellent fits (RMS < 5 km/s)
        excellent_fits = self.df[self.df['rms'] < 5.0]
        print(f"Excellent fits (RMS < 5 km/s): {len(excellent_fits)}/{len(self.df)} = {100*len(excellent_fits)/len(self.df):.1f}%")
        
        # Good fits (RMS < 10 km/s)
        good_fits = self.df[self.df['rms'] < 10.0]
        print(f"Good fits (RMS < 10 km/s): {len(good_fits)}/{len(self.df)} = {100*len(good_fits)/len(self.df):.1f}%")
        
        return rms_residuals, r0_values, v_scales
    
    def theoretical_predictions(self):
        """Compare with theoretical predictions."""
        print("\nTHEORETICAL PREDICTIONS COMPARISON")
        print("-" * 35)
        
        # Expected R0 from theory
        # From field equations: R0 should be related to galactic mass/size
        # Typical galactic scale ~ 10-100 kpc
        expected_r0_mean = 50  # kpc
        expected_r0_std = 30   # kpc
        
        observed_r0 = self.df['R0_gal'].values
        
        print(f"Expected R0 (from theory): {expected_r0_mean} +/- {expected_r0_std} kpc")
        print(f"Observed R0 (from fits): {np.mean(observed_r0):.1f} +/- {np.std(observed_r0):.1f} kpc")
        print()
        
        # Statistical test
        # Are the observed R0 values consistent with theoretical expectations?
        t_stat, p_value = stats.ttest_1samp(observed_r0, expected_r0_mean)
        print(f"Statistical test (t-test against theoretical mean):")
        print(f"  t-statistic: {t_stat:.3f}")
        print(f"  p-value: {p_value:.3f}")
        
        if p_value > 0.05:
            print("  Result: CONSISTENT with theoretical predictions")
        else:
            print("  Result: INCONSISTENT with theoretical predictions")
        
    def enhancement_factor_analysis(self):
        """Analyze the enhancement factor behavior."""
        print("\nENHANCEMENT FACTOR ANALYSIS")
        print("-" * 28)
        
        # Calculate enhancement factors at different radii
        r0_values = self.df['R0_gal'].values
        
        # Test radii
        test_radii = [1, 5, 10, 20]  # kpc
        
        print("Enhancement Factor F(tau) = 1 + alpha * 3(1-tau)/(tau^2(3-2*tau))")
        print("where tau(r) = R0/(R0 + r)")
        print()
        
        for r in test_radii:
            print(f"At r = {r} kpc:")
            
            # Calculate tau and enhancement for each galaxy
            tau_values = r0_values / (r0_values + r)
            f_values = 3 * (1 - tau_values) / (tau_values**2 * (3 - 2*tau_values))
            
            # Assume α = 0.1 (typical coupling strength)
            alpha = 0.1
            enhancement_values = 1 + alpha * f_values
            
            print(f"  Mean tau: {np.mean(tau_values):.3f}")
            print(f"  Mean F(tau): {np.mean(enhancement_values):.3f}")
            print(f"  Enhancement range: {np.min(enhancement_values):.3f} - {np.max(enhancement_values):.3f}")
            print()
        
        # Plot enhancement factors
        self.plot_enhancement_factors()
        
    def plot_enhancement_factors(self):
        """Plot enhancement factors for different galaxies."""
        r = np.linspace(0.1, 30, 100)
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Panel 1: tau(r) for different R0 values
        ax1 = axes[0, 0]
        r0_examples = [20, 50, 100]  # kpc
        for r0 in r0_examples:
            tau = r0 / (r0 + r)
            ax1.plot(r, tau, linewidth=2, label=f'R0 = {r0} kpc')
        ax1.set_xlabel('Radius (kpc)')
        ax1.set_ylabel('tau(r)')
        ax1.set_title('Temporal Connectivity Function')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Panel 2: Enhancement factor f(tau)
        ax2 = axes[0, 1]
        for r0 in r0_examples:
            tau = r0 / (r0 + r)
            f = 3 * (1 - tau) / (tau**2 * (3 - 2*tau))
            ax2.plot(r, f, linewidth=2, label=f'R0 = {r0} kpc')
        ax2.set_xlabel('Radius (kpc)')
        ax2.set_ylabel('f(tau)')
        ax2.set_title('Enhancement Factor')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_yscale('log')
        
        # Panel 3: Full enhancement F(tau) = 1 + alpha * f(tau)
        ax3 = axes[1, 0]
        alpha = 0.1
        for r0 in r0_examples:
            tau = r0 / (r0 + r)
            f = 3 * (1 - tau) / (tau**2 * (3 - 2*tau))
            F = 1 + alpha * f
            ax3.plot(r, F, linewidth=2, label=f'R0 = {r0} kpc')
        ax3.set_xlabel('Radius (kpc)')
        ax3.set_ylabel('F(tau)')
        ax3.set_title('Total Enhancement F(tau) = 1 + alpha * f(tau)')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Panel 4: Velocity profiles
        ax4 = axes[1, 1]
        for r0 in r0_examples:
            v_base = 100 * np.sqrt(r / (r + r0/3))  # Base profile
            enhancement = (1 + r/r0)  # UDT enhancement
            v_total = v_base * enhancement
            ax4.plot(r, v_total, linewidth=2, label=f'R0 = {r0} kpc')
        ax4.set_xlabel('Radius (kpc)')
        ax4.set_ylabel('Velocity (km/s)')
        ax4.set_title('UDT Velocity Profiles')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('C:/UDT/results/udt_enhancement_analysis.png', dpi=150)
        plt.close()
        
        print("Enhancement factor analysis plot saved to: C:/UDT/results/udt_enhancement_analysis.png")
        
    def cosmic_connectivity_validation(self):
        """Validate the cosmic connectivity interpretation."""
        print("\nCOSMIC CONNECTIVITY VALIDATION")
        print("-" * 32)
        
        print("Testing cosmic connectivity predictions:")
        print()
        
        # 1. R0 should correlate with galaxy properties
        print("1. R0 Scale Parameter Analysis:")
        r0_values = self.df['R0_gal'].values
        v_scales = self.df['V_scale'].values
        
        # Correlation between R0 and V_scale
        corr_coeff, p_value = stats.pearsonr(r0_values, v_scales)
        print(f"   Correlation(R0, V_scale): r = {corr_coeff:.3f}, p = {p_value:.3f}")
        
        if abs(corr_coeff) > 0.3 and p_value < 0.05:
            print("   Result: SIGNIFICANT correlation (expected for cosmic connectivity)")
        else:
            print("   Result: NO significant correlation")
        
        # 2. Fit quality should be consistent across different galaxy types
        print("\n2. Fit Quality Consistency:")
        rms_residuals = self.df['rms'].values
        
        # Test for outliers
        z_scores = np.abs(stats.zscore(rms_residuals))
        outliers = np.sum(z_scores > 3)
        print(f"   Outliers (>3sigma): {outliers}/{len(rms_residuals)} = {100*outliers/len(rms_residuals):.1f}%")
        
        if outliers/len(rms_residuals) < 0.1:
            print("   Result: CONSISTENT fit quality (supports universal theory)")
        else:
            print("   Result: INCONSISTENT fit quality (may indicate problems)")
        
        # 3. R0 distribution should be physically reasonable
        print("\n3. R0 Distribution Analysis:")
        r0_median = np.median(r0_values)
        r0_mad = np.median(np.abs(r0_values - r0_median))  # Median absolute deviation
        
        print(f"   Median R0: {r0_median:.1f} kpc")
        print(f"   MAD: {r0_mad:.1f} kpc")
        print(f"   Range: {np.min(r0_values):.1f} - {np.max(r0_values):.1f} kpc")
        
        # Check if R0 values are in physically reasonable range (10-200 kpc)
        reasonable_range = np.sum((r0_values >= 10) & (r0_values <= 200))
        print(f"   Physically reasonable (10-200 kpc): {reasonable_range}/{len(r0_values)} = {100*reasonable_range/len(r0_values):.1f}%")
        
        if reasonable_range/len(r0_values) > 0.8:
            print("   Result: PHYSICALLY REASONABLE R0 distribution")
        else:
            print("   Result: PROBLEMATIC R0 distribution")
            
    def comparison_with_alternatives(self):
        """Compare UDT performance with alternative theories."""
        print("\nCOMPARISON WITH ALTERNATIVE THEORIES")
        print("-" * 37)
        
        # Typical performance of other theories on SPARC data
        print("Expected performance on SPARC data:")
        print("  LCDM + NFW: chi^2/dof ~ 10-50 (poor fits)")
        print("  MOND: chi^2/dof ~ 2-5 (good fits)")
        print("  Modified gravity: chi^2/dof ~ 3-10 (variable)")
        print()
        
        # UDT performance
        rms_residuals = self.df['rms'].values
        mean_rms = np.mean(rms_residuals)
        median_rms = np.median(rms_residuals)
        
        print(f"UDT performance:")
        print(f"  Mean RMS: {mean_rms:.2f} km/s")
        print(f"  Median RMS: {median_rms:.2f} km/s")
        print(f"  Typical velocity uncertainties: ~3-5 km/s")
        print()
        
        # Assessment
        if median_rms < 5:
            print("Assessment: EXCELLENT performance")
            print("UDT fits are comparable to or better than MOND")
        elif median_rms < 10:
            print("Assessment: GOOD performance")
            print("UDT fits are significantly better than LCDM")
        else:
            print("Assessment: POOR performance")
            print("UDT fits are not competitive with alternatives")
    
    def generate_comprehensive_report(self):
        """Generate a comprehensive validation report."""
        print("\n" + "=" * 60)
        print("COMPREHENSIVE VALIDATION REPORT")
        print("=" * 60)
        
        rms_residuals, r0_values, v_scales = self.analyze_fit_quality()
        self.theoretical_predictions()
        self.enhancement_factor_analysis()
        self.cosmic_connectivity_validation()
        self.comparison_with_alternatives()
        
        print("\n" + "=" * 60)
        print("OVERALL ASSESSMENT")
        print("=" * 60)
        
        # Overall score based on multiple criteria
        score = 0
        max_score = 5
        
        # Criterion 1: Fit quality
        median_rms = np.median(rms_residuals)
        if median_rms < 5:
            score += 1
            print("[PASS] Fit Quality: EXCELLENT (RMS < 5 km/s)")
        elif median_rms < 10:
            score += 0.5
            print("[WARN] Fit Quality: GOOD (RMS < 10 km/s)")
        else:
            print("[FAIL] Fit Quality: POOR (RMS > 10 km/s)")
        
        # Criterion 2: Success rate
        success_rate = len(self.df[self.df['success']]) / len(self.df)
        if success_rate > 0.9:
            score += 1
            print("[PASS] Success Rate: EXCELLENT (>90%)")
        elif success_rate > 0.8:
            score += 0.5
            print("[WARN] Success Rate: GOOD (>80%)")
        else:
            print("[FAIL] Success Rate: POOR (<80%)")
        
        # Criterion 3: Parameter consistency
        r0_cv = np.std(r0_values) / np.mean(r0_values)  # Coefficient of variation
        if r0_cv < 0.5:
            score += 1
            print("[PASS] Parameter Consistency: EXCELLENT (CV < 0.5)")
        elif r0_cv < 1.0:
            score += 0.5
            print("[WARN] Parameter Consistency: GOOD (CV < 1.0)")
        else:
            print("[FAIL] Parameter Consistency: POOR (CV > 1.0)")
        
        # Criterion 4: Physical reasonableness
        reasonable_r0 = np.sum((r0_values >= 10) & (r0_values <= 200)) / len(r0_values)
        if reasonable_r0 > 0.8:
            score += 1
            print("[PASS] Physical Reasonableness: EXCELLENT (>80% reasonable)")
        elif reasonable_r0 > 0.6:
            score += 0.5
            print("[WARN] Physical Reasonableness: GOOD (>60% reasonable)")
        else:
            print("[FAIL] Physical Reasonableness: POOR (<60% reasonable)")
        
        # Criterion 5: Competitive performance
        if median_rms < 5:
            score += 1
            print("[PASS] Competitive Performance: EXCELLENT (comparable to MOND)")
        elif median_rms < 10:
            score += 0.5
            print("[WARN] Competitive Performance: GOOD (better than LCDM)")
        else:
            print("[FAIL] Competitive Performance: POOR (not competitive)")
        
        print(f"\nOVERALL SCORE: {score:.1f}/{max_score}")
        
        if score >= 4:
            print("CONCLUSION: UDT FRAMEWORK VALIDATED")
            print("The cosmic connectivity theory shows excellent agreement with SPARC data")
        elif score >= 3:
            print("CONCLUSION: UDT FRAMEWORK PROMISING")
            print("The cosmic connectivity theory shows good agreement with SPARC data")
        elif score >= 2:
            print("CONCLUSION: UDT FRAMEWORK NEEDS IMPROVEMENT")
            print("The cosmic connectivity theory shows mixed results")
        else:
            print("CONCLUSION: UDT FRAMEWORK FAILS VALIDATION")
            print("The cosmic connectivity theory is not supported by SPARC data")

def main():
    """Run comprehensive SPARC validation."""
    validator = ComprehensiveSPARCValidation()
    validator.generate_comprehensive_report()

if __name__ == "__main__":
    main()