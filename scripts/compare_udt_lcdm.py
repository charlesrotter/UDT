#!/usr/bin/env python3
"""
Direct UDT vs LCDM Comparison
=============================

Direct comparison of UDT and LCDM models on identical raw supernova data.
Uses same data and fitting procedures to determine which model fits better.

Author: Charles Rotter
Date: 2025-01-17
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.integrate import quad
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from udt.utils.data_loader import load_csp_database, load_pantheon_data


def udt_distance_modulus(z, R0, M_B):
    """UDT distance modulus: mu = m - M = 5*log10(d_L) + 25"""
    d_L_Mpc = z * R0  # UDT distance relation
    mu = 5 * np.log10(d_L_Mpc) + 25
    return mu + M_B


def lcdm_distance_modulus(z, H0, Omega_m, M_B):
    """LCDM distance modulus with flat geometry (Omega_L = 1 - Omega_m)"""
    c = 299792.458  # km/s
    Omega_lambda = 1 - Omega_m
    
    def E(z_prime):
        return 1 / np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_lambda)
    
    # Integrate E(z) from 0 to z
    d_L_Mpc = np.zeros_like(z)
    for i, z_val in enumerate(z):
        if z_val > 0:
            integral, _ = quad(E, 0, z_val)
            d_L_Mpc[i] = (c / H0) * (1 + z_val) * integral
        else:
            d_L_Mpc[i] = c * z_val / H0  # Linear approximation for z->0
    
    mu = 5 * np.log10(d_L_Mpc) + 25
    return mu + M_B


def fit_udt_model(z_obs, m_obs, m_err):
    """Fit UDT model to supernova data"""
    def chi_squared(params):
        R0, M_B = params
        m_pred = udt_distance_modulus(z_obs, R0, M_B)
        residuals = (m_obs - m_pred) / m_err
        return np.sum(residuals**2)
    
    # Initial guess
    initial_guess = [3000, -19.0]
    bounds = [(100, 10000), (-22, -17)]
    
    result = minimize(chi_squared, initial_guess, bounds=bounds, method='L-BFGS-B')
    
    if result.success:
        R0_fit, M_B_fit = result.x
        chi2_min = result.fun
        m_pred = udt_distance_modulus(z_obs, R0_fit, M_B_fit)
        residuals = m_obs - m_pred
        rms = np.sqrt(np.mean(residuals**2))
        
        return {
            'success': True,
            'R0': R0_fit,
            'M_B': M_B_fit,
            'chi2': chi2_min,
            'rms': rms,
            'residuals': residuals,
            'm_predicted': m_pred
        }
    else:
        return {'success': False}


def fit_lcdm_model(z_obs, m_obs, m_err):
    """Fit LCDM model to supernova data"""
    def chi_squared(params):
        H0, Omega_m, M_B = params
        m_pred = lcdm_distance_modulus(z_obs, H0, Omega_m, M_B)
        residuals = (m_obs - m_pred) / m_err
        return np.sum(residuals**2)
    
    # Initial guess (Planck 2018 values)
    initial_guess = [70, 0.3, -19.0]
    bounds = [(50, 100), (0.1, 0.9), (-22, -17)]
    
    result = minimize(chi_squared, initial_guess, bounds=bounds, method='L-BFGS-B')
    
    if result.success:
        H0_fit, Omega_m_fit, M_B_fit = result.x
        chi2_min = result.fun
        m_pred = lcdm_distance_modulus(z_obs, H0_fit, Omega_m_fit, M_B_fit)
        residuals = m_obs - m_pred
        rms = np.sqrt(np.mean(residuals**2))
        
        return {
            'success': True,
            'H0': H0_fit,
            'Omega_m': Omega_m_fit,
            'M_B': M_B_fit,
            'chi2': chi2_min,
            'rms': rms,
            'residuals': residuals,
            'm_predicted': m_pred
        }
    else:
        return {'success': False}


def compare_models(dataset_name, z_obs, m_obs, m_err):
    """Compare UDT and LCDM models on identical data"""
    print(f"\n{'='*60}")
    print(f"UDT vs LCDM COMPARISON - {dataset_name}")
    print(f"{'='*60}")
    print(f"Sample size: {len(z_obs)} supernovae")
    print(f"Redshift range: {z_obs.min():.4f} - {z_obs.max():.4f}")
    print(f"Data type: RAW observational magnitudes")
    print()
    
    # Fit UDT model
    print("Fitting UDT model...")
    udt_result = fit_udt_model(z_obs, m_obs, m_err)
    
    # Fit LCDM model
    print("Fitting LCDM model...")
    lcdm_result = fit_lcdm_model(z_obs, m_obs, m_err)
    
    if udt_result['success'] and lcdm_result['success']:
        # Calculate degrees of freedom
        dof_udt = len(z_obs) - 2  # R0, M_B
        dof_lcdm = len(z_obs) - 3  # H0, Omega_m, M_B
        
        print(f"\nUDT RESULTS:")
        print(f"  R0 = {udt_result['R0']:.0f} Mpc")
        print(f"  M_B = {udt_result['M_B']:.3f}")
        print(f"  chi^2 = {udt_result['chi2']:.1f}")
        print(f"  chi^2/dof = {udt_result['chi2']/dof_udt:.2f}")
        print(f"  RMS = {udt_result['rms']:.3f} mag")
        
        print(f"\nLCDM RESULTS:")
        print(f"  H0 = {lcdm_result['H0']:.1f} km/s/Mpc")
        print(f"  Omega_m = {lcdm_result['Omega_m']:.3f}")
        print(f"  M_B = {lcdm_result['M_B']:.3f}")
        print(f"  chi^2 = {lcdm_result['chi2']:.1f}")
        print(f"  chi^2/dof = {lcdm_result['chi2']/dof_lcdm:.2f}")
        print(f"  RMS = {lcdm_result['rms']:.3f} mag")
        
        # Model comparison
        delta_chi2 = lcdm_result['chi2'] - udt_result['chi2']
        rms_improvement = (lcdm_result['rms'] - udt_result['rms']) / lcdm_result['rms'] * 100
        
        print(f"\nMODEL COMPARISON:")
        print(f"  Delta chi^2 = {delta_chi2:+.1f} (LCDM - UDT)")
        print(f"  RMS improvement = {rms_improvement:+.1f}% (UDT better if positive)")
        
        # Statistical significance
        if delta_chi2 > 10:
            verdict = "UDT STRONGLY PREFERRED"
        elif delta_chi2 > 5:
            verdict = "UDT PREFERRED"
        elif delta_chi2 > 0:
            verdict = "UDT SLIGHTLY PREFERRED"
        elif delta_chi2 > -5:
            verdict = "COMPARABLE FIT QUALITY"
        elif delta_chi2 > -10:
            verdict = "LCDM SLIGHTLY PREFERRED"
        else:
            verdict = "LCDM STRONGLY PREFERRED"
        
        print(f"  VERDICT: {verdict}")
        
        # Additional analysis
        print(f"\nADDITIONAL ANALYSIS:")
        print(f"  UDT effective H0 = {299792.458 / udt_result['R0']:.1f} km/s/Mpc")
        print(f"  LCDM best-fit H0 = {lcdm_result['H0']:.1f} km/s/Mpc")
        print(f"  Hubble constant agreement: {abs(299792.458/udt_result['R0'] - lcdm_result['H0']):.1f} km/s/Mpc difference")
        
        return {
            'udt': udt_result,
            'lcdm': lcdm_result,
            'delta_chi2': delta_chi2,
            'rms_improvement': rms_improvement,
            'verdict': verdict
        }
    else:
        print("ERROR: One or both model fits failed")
        return None


def main():
    print("="*80)
    print("DIRECT UDT vs LCDM COMPARISON")
    print("="*80)
    print("Testing identical raw data with both models")
    print("Uses proper cosmological integration for LCDM")
    print("="*80)
    
    results = {}
    
    # Test CSP data
    print("\nLoading CSP DR3 data...")
    csp_dir = 'data/CSP_Photometry_DR3/DR3'
    try:
        csp_df = load_csp_database(csp_dir)
        if len(csp_df) > 0:
            # Use raw data
            csp_df = csp_df[csp_df['redshift'] <= 0.1]  # z cut
            z_csp = csp_df['redshift'].values
            m_csp = csp_df['B_peak_raw'].values
            m_err_csp = csp_df['B_error_raw'].values
            
            results['csp'] = compare_models('CSP DR3', z_csp, m_csp, m_err_csp)
    except Exception as e:
        print(f"Error loading CSP data: {e}")
    
    # Test Pantheon+ data
    print("\nLoading Pantheon+ data...")
    pantheon_file = 'data/Pantheon_SH0ES.dat'
    try:
        pantheon_df = load_pantheon_data(pantheon_file)
        if len(pantheon_df) > 0:
            # Use raw data
            pantheon_df = pantheon_df[pantheon_df['zCMB'] <= 0.1]  # z cut
            z_pantheon = pantheon_df['zCMB'].values
            m_pantheon = pantheon_df['mB'].values  # Raw SALT2 magnitude
            m_err_pantheon = pantheon_df['mBERR'].values
            
            results['pantheon'] = compare_models('Pantheon+', z_pantheon, m_pantheon, m_err_pantheon)
    except Exception as e:
        print(f"Error loading Pantheon+ data: {e}")
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY: UDT vs LCDM PERFORMANCE")
    print("="*80)
    
    for dataset_name, result in results.items():
        if result:
            print(f"\n{dataset_name.upper()}:")
            print(f"  Delta chi^2 = {result['delta_chi2']:+.1f}")
            print(f"  RMS improvement = {result['rms_improvement']:+.1f}%")
            print(f"  Verdict: {result['verdict']}")
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print("This comparison uses identical raw data for both models")
    print("Positive Delta chi^2 means UDT fits better than LCDM")
    print("Positive RMS improvement means UDT has lower scatter")
    print("="*80)


if __name__ == "__main__":
    main()