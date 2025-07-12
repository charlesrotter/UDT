"""
Real SPARC Data Analysis for Information Curvature Theory
Test β = 2.5 hypothesis against actual galaxy rotation curves
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

class RealSPARCAnalysis:
    """
    Analyze real SPARC rotation curves for β parameter
    """
    
    def __init__(self):
        self.G = 4.3e-6  # (km/s)² kpc / M☉
        
    def load_rotation_curve(self, galaxy_name):
        """
        Load rotation curve data for specific galaxy
        """
        file_path = f"data/sparc_curves/{galaxy_name}_rotcur.dat"
        
        if not os.path.exists(file_path):
            print(f"❌ File not found: {file_path}")
            return None
        
        try:
            # Load SPARC rotation curve format
            # Columns typically: Radius, Vobs, errV, Vgas, Vdisk, Vbul
            data = pd.read_csv(file_path, sep='\s+', comment='#',
                             names=['R_kpc', 'V_obs', 'V_err', 'V_gas', 'V_disk', 'V_bul'])
            
            # Quality control
            valid_data = data[
                (data['R_kpc'] > 0) &
                (data['V_obs'] > 0) &
                (data['V_err'] > 0) &
                (data['V_err'] < 0.5 * data['V_obs'])  # Reasonable errors
            ].copy()
            
            print(f"✅ {galaxy_name}: {len(valid_data)} data points")
            return valid_data
            
        except Exception as e:
            print(f"❌ Error loading {galaxy_name}: {e}")
            return None
    
    def information_curvature_model(self, r, GM, kappa_c2, R0, beta):
        """
        Information Curvature rotation curve model
        v²(r) = GM/r + κc² (r/R₀)^(β-1)
        """
        r = np.maximum(r, 1e-6)
        R0 = max(R0, 1e-6)
        
        v_newton_sq = GM / r
        v_info_sq = kappa_c2 * (r / R0)**(beta - 1)
        
        return np.sqrt(v_newton_sq + v_info_sq)
    
    def fit_real_galaxy(self, galaxy_name):
        """
        Fit Information Curvature model to real galaxy data
        """
        print(f"\n🔬 FITTING {galaxy_name}")
        print("-" * 25)
        
        # Load data
        data = self.load_rotation_curve(galaxy_name)
        if data is None:
            return None
        
        r_kpc = data['R_kpc'].values
        v_obs = data['V_obs'].values
        v_err = data['V_err'].values
        
        # Initial parameter guesses
        GM_guess = np.mean(v_obs**2 * r_kpc)
        kappa_c2_guess = np.mean(v_obs**2)
        R0_guess = np.median(r_kpc)
        beta_guess = 2.5  # Our hypothesis
        
        initial_guess = [GM_guess, kappa_c2_guess, R0_guess, beta_guess]
        
        # Realistic parameter bounds
        bounds = (
            [1e8, 1e2, 0.1, 1.0],    # Lower bounds
            [1e14, 1e7, 100, 5.0]    # Upper bounds
        )
        
        try:
            # Weighted least squares fit
            popt, pcov = curve_fit(
                self.information_curvature_model,
                r_kpc, v_obs, sigma=v_err,
                p0=initial_guess, bounds=bounds,
                maxfev=10000
            )
            
            # Extract results
            GM, kappa_c2, R0, beta = popt
            param_errors = np.sqrt(np.diag(pcov))
            GM_err, kappa_err, R0_err, beta_err = param_errors
            
            # Calculate fit quality
            v_model = self.information_curvature_model(r_kpc, *popt)
            residuals = (v_obs - v_model) / v_err
            chi2 = np.sum(residuals**2)
            dof = len(r_kpc) - 4
            chi2_reduced = chi2 / dof if dof > 0 else np.inf
            
            print(f"β = {beta:.3f} ± {beta_err:.3f}")
            print(f"χ²/ν = {chi2_reduced:.2f}")
            print(f"R₀ = {R0:.1f} ± {R0_err:.1f} kpc")
            
            # Assess fit quality
            if chi2_reduced < 2.0 and beta_err < 0.5:
                quality = "excellent"
                print("Quality: EXCELLENT ✅")
            elif chi2_reduced < 4.0 and beta_err < 1.0:
                quality = "good"
                print("Quality: GOOD ✓")
            else:
                quality = "poor"
                print("Quality: POOR ⚠️")
            
            return {
                'galaxy': galaxy_name,
                'success': True,
                'beta': beta,
                'beta_error': beta_err,
                'GM': GM,
                'kappa_c2': kappa_c2,
                'R0': R0,
                'chi2_reduced': chi2_reduced,
                'quality': quality,
                'n_points': len(r_kpc),
                'data': {
                    'r_kpc': r_kpc,
                    'v_obs': v_obs,
                    'v_err': v_err,
                    'v_model': v_model
                }
            }
            
        except Exception as e:
            print(f"❌ Fit failed: {e}")
            return {
                'galaxy': galaxy_name,
                'success': False,
                'error': str(e)
            }
    
    def analyze_real_sample(self, galaxy_list):
        """
        Analyze complete sample of real galaxies
        """
        print("🌌 REAL SPARC ANALYSIS")
        print("Information Curvature Theory - Actual Data Test")
        print("=" * 50)
        
        results = []
        successful_fits = []
        
        for galaxy in galaxy_list:
            result = self.fit_real_galaxy(galaxy)
            results.append(result)
            
            if result and result['success']:
                successful_fits.append(result)
        
        if len(successful_fits) < 3:
            print(f"\n❌ INSUFFICIENT SUCCESSFUL FITS")
            print(f"Only {len(successful_fits)} galaxies fitted successfully")
            print("Need at least 3 for statistical analysis")
            return results
        
        # Statistical analysis of β measurements
        beta_values = [r['beta'] for r in successful_fits]
        beta_errors = [r['beta_error'] for r in successful_fits]
        galaxy_names = [r['galaxy'] for r in successful_fits]
        
        print(f"\n📊 REAL DATA STATISTICAL ANALYSIS")
        print("=" * 35)
        print(f"Successful fits: {len(successful_fits)} galaxies")
        
        # Calculate statistics
        beta_mean = np.mean(beta_values)
        beta_std = np.std(beta_values, ddof=1)
        beta_sem = beta_std / np.sqrt(len(beta_values))
        
        print(f"Measured β: {beta_mean:.3f} ± {beta_sem:.3f}")
        print(f"Standard deviation: {beta_std:.3f}")
        print(f"Range: [{np.min(beta_values):.3f}, {np.max(beta_values):.3f}]")
        
        # Test against β = 2.5 hypothesis
        hypothesis_beta = 2.5
        t_statistic = (beta_mean - hypothesis_beta) / beta_sem
        deviation_sigma = abs(t_statistic)
        
        print(f"\nHypothesis test (β = 2.5):")
        print(f"t-statistic: {t_statistic:.2f}")
        print(f"Deviation: {deviation_sigma:.1f}σ")
        
        # Assessment
        if deviation_sigma < 1.0:
            status = "STRONG SUPPORT"
            print("🎊 STRONG SUPPORT for β = 2.5!")
        elif deviation_sigma < 2.0:
            status = "CONSISTENT"
            print("✅ CONSISTENT with β = 2.5")
        elif deviation_sigma < 3.0:
            status = "MARGINAL"
            print("⚠️ MARGINAL consistency with β = 2.5")
        else:
            status = "INCONSISTENT"
            print("❌ INCONSISTENT with β = 2.5")
            print(f"Real β appears to be ~{beta_mean:.2f}")
        
        print(f"\n📋 INDIVIDUAL GALAXY RESULTS:")
        for i, result in enumerate(successful_fits, 1):
            print(f"{i:2d}. {result['galaxy']}: β = {result['beta']:.3f} ± {result['beta_error']:.3f}")
        
        print(f"\n🎯 CONCLUSION:")
        print(f"Status: {status}")
        if status == "INCONSISTENT":
            print(f"🔄 THEORY REFINEMENT NEEDED")
            print(f"Consider updating hypothesis to β = {beta_mean:.2f}")
        else:
            print(f"✅ Information Curvature Theory shows promise!")
        
        return {
            'results': results,
            'successful_fits': successful_fits,
            'beta_mean': beta_mean,
            'beta_sem': beta_sem,
            'beta_std': beta_std,
            'status': status,
            'deviation_sigma': deviation_sigma
        }

def main():
    """
    Main analysis function
    """
    # List of galaxies to analyze (update based on successful downloads)
    galaxies = [
        'NGC2403', 'NGC3198', 'UGC02885', 'DDO154', 'NGC7793',
        'IC2574', 'NGC925', 'DDO170', 'NGC1560', 'NGC5055'
    ]
    
    analyzer = RealSPARCAnalysis()
    results = analyzer.analyze_real_sample(galaxies)
    
    return results

if __name__ == "__main__":
    analysis_results = main()
