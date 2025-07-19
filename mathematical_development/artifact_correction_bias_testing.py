#!/usr/bin/env python
"""
Artifact Correction Bias Testing Framework
==========================================

CRITICAL SCIENTIFIC INTEGRITY TEST:
This script validates that artifact correction methods do not introduce
pro-UDT bias. It ensures corrections are scientifically rigorous and
actually help determine which cosmological model is correct.

REQUIREMENTS:
1. Corrections applied to ΛCDM data must still favor ΛCDM
2. Corrections applied to UDT data must still favor UDT  
3. Systematic uncertainties must increase, not decrease
4. Control tests must demonstrate method robustness

Author: Charles Rotter
Scientific Purpose: Prevent methodological bias in cosmological model testing
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import sys
import os

# Add parent directory for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from udt.core.temporal_geometry import distance_from_redshift

# Embedded cosmology classes for testing
class UDTCosmology:
    """UDT cosmology for bias testing"""
    def __init__(self, R0):
        self.R0 = R0
        
    def luminosity_distance(self, z):
        """Luminosity distance in UDT framework"""
        return distance_from_redshift(z, self.R0)
        
    def distance_modulus(self, z):
        """Distance modulus in UDT framework"""
        d_L = self.luminosity_distance(z)
        return 5 * np.log10(d_L) + 25  # d_L in Mpc

class LCDMCosmology:
    """ΛCDM cosmology for bias testing"""
    def __init__(self, H0, Omega_m=0.3, Omega_Lambda=0.7):
        self.H0 = H0
        self.Omega_m = Omega_m
        self.Omega_Lambda = Omega_Lambda
        
    def luminosity_distance(self, z):
        """Luminosity distance in flat ΛCDM"""
        from scipy.integrate import quad
        
        def E(z_prime):
            return np.sqrt(self.Omega_m * (1 + z_prime)**3 + self.Omega_Lambda)
        
        # Comoving distance integral
        integral, _ = quad(lambda z_prime: 1/E(z_prime), 0, z)
        
        # Convert to distance in Mpc
        c_km_s = 299792.458
        d_c = c_km_s * integral / self.H0  # Mpc
        
        # Luminosity distance
        d_L = d_c * (1 + z)
        return d_L
        
    def distance_modulus(self, z):
        """Distance modulus in ΛCDM framework"""
        d_L = self.luminosity_distance(z)
        return 5 * np.log10(d_L) + 25  # d_L in Mpc

class ArtifactCorrectionBiasTester:
    """
    Comprehensive bias testing for artifact correction methods
    """
    
    def __init__(self):
        self.results = {}
        
    def generate_synthetic_contaminated_data(self, true_cosmology, contaminating_cosmology, 
                                           contamination_level=0.1, n_points=100):
        """
        Generate synthetic data with known contamination
        
        Parameters:
        -----------
        true_cosmology : Cosmology object
            The actual cosmology that generated the data
        contaminating_cosmology : Cosmology object  
            The cosmology used in data processing (introduces bias)
        contamination_level : float
            Magnitude of contamination (0.1 = 10% systematic shift)
        n_points : int
            Number of data points to generate
            
        Returns:
        --------
        contaminated_data : dict
            Synthetic dataset with known contamination
        """
        
        # Generate redshift range
        z_range = np.linspace(0.01, 0.8, n_points)
        
        # Generate true distance moduli
        true_mu = []
        for z in z_range:
            if hasattr(true_cosmology, 'distance_modulus'):
                mu_true = true_cosmology.distance_modulus(z)
            else:
                # Fallback calculation
                d_L = true_cosmology.luminosity_distance(z)
                mu_true = 5 * np.log10(d_L) + 25  # d_L in Mpc
            true_mu.append(mu_true)
        
        true_mu = np.array(true_mu)
        
        # Apply contamination bias (systematic shift based on contaminating cosmology)
        contaminated_mu = []
        for z in z_range:
            if hasattr(contaminating_cosmology, 'distance_modulus'):
                mu_contam = contaminating_cosmology.distance_modulus(z)
            else:
                d_L = contaminating_cosmology.luminosity_distance(z)
                mu_contam = 5 * np.log10(d_L) + 25
            contaminated_mu.append(mu_contam)
        
        contaminated_mu = np.array(contaminated_mu)
        
        # Apply contamination as systematic shift
        bias_correction = contamination_level * (contaminated_mu - true_mu)
        observed_mu = true_mu + bias_correction
        
        # Add realistic uncertainties
        uncertainties = 0.15 + 0.1 * z_range  # Realistic supernova uncertainties
        
        # Add random scatter
        observed_mu += np.random.normal(0, uncertainties)
        
        return {
            'redshift': z_range,
            'distance_modulus': observed_mu,
            'uncertainty': uncertainties,
            'true_mu': true_mu,
            'contamination_bias': bias_correction,
            'true_cosmology': true_cosmology.__class__.__name__,
            'contaminating_cosmology': contaminating_cosmology.__class__.__name__,
            'contamination_level': contamination_level
        }
    
    def apply_artifact_correction(self, data, correction_method='reverse_engineering'):
        """
        Apply artifact correction method to contaminated data
        
        CRITICAL: This correction must not introduce pro-UDT bias!
        """
        
        if correction_method == 'reverse_engineering':
            # Simulate reverse-engineering approach
            # Remove systematic bias without knowledge of true cosmology
            
            # Estimate contamination bias from redshift-dependent patterns
            z = data['redshift']
            mu_obs = data['distance_modulus']
            uncertainty = data['uncertainty']
            
            # Model-independent contamination detection
            # Look for systematic trends vs redshift that indicate processing bias
            
            # Fit polynomial to residuals to identify systematic bias
            from numpy.polynomial import Polynomial
            p = Polynomial.fit(z, mu_obs, deg=2)
            systematic_trend = p(z)
            
            # Remove systematic trend (this is the "artifact correction")
            # Conservative approach: only remove clear systematic bias
            detrended_mu = mu_obs - (systematic_trend - np.mean(systematic_trend))
            
            # Inflate uncertainties to account for correction uncertainty
            inflated_uncertainty = np.sqrt(uncertainty**2 + (0.05)**2)  # Add 0.05 mag systematic
            
            corrected_data = data.copy()
            corrected_data['distance_modulus'] = detrended_mu
            corrected_data['uncertainty'] = inflated_uncertainty
            corrected_data['correction_applied'] = systematic_trend - np.mean(systematic_trend)
            
            return corrected_data
        
        else:
            raise ValueError(f"Unknown correction method: {correction_method}")
    
    def fit_cosmology_to_data(self, data, cosmology_class):
        """
        Fit cosmological model to data and return fit quality
        """
        
        def chi_squared(params):
            if cosmology_class == UDTCosmology:
                R0 = params[0]
                cosmology = UDTCosmology(R0)
            elif cosmology_class == LCDMCosmology:
                H0 = params[0]
                cosmology = LCDMCosmology(H0)
            else:
                raise ValueError("Unknown cosmology class")
            
            chi2 = 0
            for i, z in enumerate(data['redshift']):
                if hasattr(cosmology, 'distance_modulus'):
                    mu_pred = cosmology.distance_modulus(z)
                else:
                    d_L = cosmology.luminosity_distance(z)
                    mu_pred = 5 * np.log10(d_L) + 25
                
                mu_obs = data['distance_modulus'][i]
                uncertainty = data['uncertainty'][i]
                
                chi2 += ((mu_obs - mu_pred) / uncertainty)**2
            
            return chi2
        
        # Initial parameter guesses
        if cosmology_class == UDTCosmology:
            initial_guess = [3500.0]  # R0 in Mpc
            bounds = [(1000, 10000)]
        elif cosmology_class == LCDMCosmology:
            initial_guess = [70.0]  # H0 in km/s/Mpc
            bounds = [(50, 100)]
        
        # Minimize chi-squared
        result = minimize(chi_squared, initial_guess, bounds=bounds)
        
        chi2_min = result.fun
        best_params = result.x
        n_data = len(data['redshift'])
        n_params = len(initial_guess)
        dof = n_data - n_params
        
        return {
            'chi2': chi2_min,
            'dof': dof,
            'chi2_per_dof': chi2_min / dof,
            'parameters': best_params,
            'fit_success': result.success
        }
    
    def control_test_lcdm_recovery(self):
        """
        CRITICAL TEST: Corrections applied to ΛCDM data must still favor ΛCDM
        
        This is the key test to ensure we don't introduce pro-UDT bias
        """
        
        print("=" * 70)
        print("CONTROL TEST: LCDM DATA RECOVERY")
        print("=" * 70)
        print("PURPOSE: Verify artifact correction doesn't bias toward UDT")
        print()
        
        # Generate LCDM data with UDT contamination
        lcdm_truth = LCDMCosmology(H0=70)
        udt_contaminant = UDTCosmology(R0=3500)
        
        contaminated_data = self.generate_synthetic_contaminated_data(
            true_cosmology=lcdm_truth,
            contaminating_cosmology=udt_contaminant,
            contamination_level=0.15,  # 15% contamination
            n_points=150
        )
        
        print(f"Generated LCDM data contaminated by UDT processing")
        print(f"True cosmology: {contaminated_data['true_cosmology']}")
        print(f"Contaminating cosmology: {contaminated_data['contaminating_cosmology']}")
        print(f"Contamination level: {contaminated_data['contamination_level']:.1%}")
        print()
        
        # Test 1: Fit contaminated data directly
        print("1. FITTING CONTAMINATED DATA (before correction):")
        lcdm_fit_raw = self.fit_cosmology_to_data(contaminated_data, LCDMCosmology)
        udt_fit_raw = self.fit_cosmology_to_data(contaminated_data, UDTCosmology)
        
        print(f"   LCDM fit: chi2/dof = {lcdm_fit_raw['chi2_per_dof']:.2f}")
        print(f"   UDT fit:  chi2/dof = {udt_fit_raw['chi2_per_dof']:.2f}")
        
        if udt_fit_raw['chi2_per_dof'] < lcdm_fit_raw['chi2_per_dof']:
            print("   -> Contamination makes UDT appear better (as expected)")
        else:
            print("   -> WARNING: Contamination didn't favor UDT as expected")
        print()
        
        # Test 2: Apply artifact correction
        corrected_data = self.apply_artifact_correction(contaminated_data)
        
        print("2. FITTING CORRECTED DATA (after artifact correction):")
        lcdm_fit_corrected = self.fit_cosmology_to_data(corrected_data, LCDMCosmology)
        udt_fit_corrected = self.fit_cosmology_to_data(corrected_data, UDTCosmology)
        
        print(f"   LCDM fit: chi2/dof = {lcdm_fit_corrected['chi2_per_dof']:.2f}")
        print(f"   UDT fit:  chi2/dof = {udt_fit_corrected['chi2_per_dof']:.2f}")
        
        # CRITICAL TEST: LCDM must fit better after correction
        if lcdm_fit_corrected['chi2_per_dof'] < udt_fit_corrected['chi2_per_dof']:
            print("   [PASS] Correction properly favors true LCDM cosmology")
            bias_test_pass = True
        else:
            print("   [FAIL] Correction introduces pro-UDT bias!")
            bias_test_pass = False
        print()
        
        # Test 3: Uncertainty inflation check
        uncertainty_ratio = np.mean(corrected_data['uncertainty']) / np.mean(contaminated_data['uncertainty'])
        print("3. UNCERTAINTY INFLATION CHECK:")
        print(f"   Average uncertainty before correction: {np.mean(contaminated_data['uncertainty']):.3f}")
        print(f"   Average uncertainty after correction:  {np.mean(corrected_data['uncertainty']):.3f}")
        print(f"   Uncertainty inflation factor: {uncertainty_ratio:.2f}")
        
        if uncertainty_ratio > 1.0:
            print("   [PASS] Uncertainties properly inflated")
            uncertainty_test_pass = True
        else:
            print("   [FAIL] Uncertainties decreased (suspicious!)")
            uncertainty_test_pass = False
        print()
        
        # Overall assessment
        print("CONTROL TEST RESULTS:")
        print(f"   Bias test: {'PASS' if bias_test_pass else 'FAIL'}")
        print(f"   Uncertainty test: {'PASS' if uncertainty_test_pass else 'FAIL'}")
        
        if bias_test_pass and uncertainty_test_pass:
            print("   [OVERALL PASS] Artifact correction methodology is scientifically valid")
        else:
            print("   [OVERALL FAIL] Artifact correction introduces bias - METHOD IS INVALID")
        
        return {
            'bias_test_pass': bias_test_pass,
            'uncertainty_test_pass': uncertainty_test_pass,
            'raw_fits': {'lcdm': lcdm_fit_raw, 'udt': udt_fit_raw},
            'corrected_fits': {'lcdm': lcdm_fit_corrected, 'udt': udt_fit_corrected},
            'uncertainty_inflation': uncertainty_ratio
        }
    
    def control_test_udt_recovery(self):
        """
        SYMMETRY TEST: Corrections applied to UDT data must still favor UDT
        
        This ensures the correction method is symmetric and not systematically biased
        """
        
        print("\n" + "=" * 70)
        print("SYMMETRY TEST: UDT DATA RECOVERY")
        print("=" * 70)
        print("PURPOSE: Verify correction method is symmetric (not systematically biased)")
        print()
        
        # Generate UDT data with ΛCDM contamination
        udt_truth = UDTCosmology(R0=3500)
        lcdm_contaminant = LCDMCosmology(H0=70)
        
        contaminated_data = self.generate_synthetic_contaminated_data(
            true_cosmology=udt_truth,
            contaminating_cosmology=lcdm_contaminant,
            contamination_level=0.15,
            n_points=150
        )
        
        print(f"Generated UDT data contaminated by LCDM processing")
        print(f"True cosmology: {contaminated_data['true_cosmology']}")
        print(f"Contaminating cosmology: {contaminated_data['contaminating_cosmology']}")
        print()
        
        # Apply correction and test
        corrected_data = self.apply_artifact_correction(contaminated_data)
        
        udt_fit_corrected = self.fit_cosmology_to_data(corrected_data, UDTCosmology)
        lcdm_fit_corrected = self.fit_cosmology_to_data(corrected_data, LCDMCosmology)
        
        print("FITTING CORRECTED UDT DATA:")
        print(f"   UDT fit:  chi2/dof = {udt_fit_corrected['chi2_per_dof']:.2f}")
        print(f"   LCDM fit: chi2/dof = {lcdm_fit_corrected['chi2_per_dof']:.2f}")
        
        if udt_fit_corrected['chi2_per_dof'] < lcdm_fit_corrected['chi2_per_dof']:
            print("   [PASS] Correction properly favors true UDT cosmology")
            symmetry_pass = True
        else:
            print("   [FAIL] Correction method is not symmetric!")
            symmetry_pass = False
        
        return {
            'symmetry_test_pass': symmetry_pass,
            'udt_fit': udt_fit_corrected,
            'lcdm_fit': lcdm_fit_corrected
        }
    
    def comprehensive_bias_test(self):
        """
        Run complete bias testing suite
        """
        
        print("ARTIFACT CORRECTION BIAS TESTING FRAMEWORK")
        print("==========================================")
        print("Scientific Purpose: Validate correction methodology doesn't introduce bias")
        print("Critical Requirement: Methods must be model-agnostic and conservative")
        print()
        
        # Run control tests
        lcdm_test = self.control_test_lcdm_recovery()
        udt_test = self.control_test_udt_recovery()
        
        # Overall assessment
        print("\n" + "=" * 70)
        print("COMPREHENSIVE BIAS TEST RESULTS")
        print("=" * 70)
        
        all_tests_pass = (lcdm_test['bias_test_pass'] and 
                         lcdm_test['uncertainty_test_pass'] and 
                         udt_test['symmetry_test_pass'])
        
        if all_tests_pass:
            print("[SCIENTIFIC VALIDATION PASS] Artifact correction methodology is VALID")
            print("   -> Corrections do not introduce pro-UDT bias")
            print("   -> Method is symmetric and conservative")
            print("   -> Systematic uncertainties are properly inflated")
            print("   -> Safe to use for cosmological model comparison")
        else:
            print("[SCIENTIFIC VALIDATION FAIL] Artifact correction methodology is INVALID")
            print("   -> Method introduces systematic bias")
            print("   -> DO NOT USE for scientific analysis")
            print("   -> Requires methodological revision")
        
        print("\nDOCUMENTATION REQUIREMENTS:")
        print("- Include these test results in all scientific publications")
        print("- Provide complete test code for peer review")
        print("- Update tests with each methodological change")
        print("- Require independent verification by collaborators")
        
        return {
            'overall_valid': all_tests_pass,
            'lcdm_recovery_test': lcdm_test,
            'udt_recovery_test': udt_test
        }

def main():
    """
    Run bias testing for artifact correction methods
    """
    
    tester = ArtifactCorrectionBiasTester()
    results = tester.comprehensive_bias_test()
    
    # Save results for documentation
    import json
    
    # Convert numpy arrays to lists for JSON serialization
    def convert_numpy(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        return obj
    
    # Clean results for JSON
    def clean_for_json(d):
        if isinstance(d, dict):
            return {k: clean_for_json(v) for k, v in d.items()}
        elif isinstance(d, list):
            return [clean_for_json(v) for v in d]
        else:
            return convert_numpy(d)
    
    clean_results = clean_for_json(results)
    
    results_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'results')
    os.makedirs(results_dir, exist_ok=True)
    
    with open(os.path.join(results_dir, 'artifact_correction_bias_test_results.json'), 'w') as f:
        json.dump(clean_results, f, indent=2)
    
    print(f"\nResults saved to: {os.path.join(results_dir, 'artifact_correction_bias_test_results.json')}")
    
    return results

if __name__ == "__main__":
    main()