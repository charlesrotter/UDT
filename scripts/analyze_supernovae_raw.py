#!/usr/bin/env python3
"""
Supernova Analysis using UDT Framework - Raw Data Only
=====================================================

Analyzes Type Ia supernovae using the Universal Distance Dilation Theory.
This version specifically uses RAW data to avoid contamination:
- CSP: Raw B-band peak magnitudes from individual light curves
- Pantheon+: Raw mB (SALT2 peak magnitude) instead of corrected m_b_corr

Author: Charles Rotter
Date: 2025-01-17
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Add parent directory to path to import udt package
sys.path.insert(0, str(Path(__file__).parent.parent))

from udt.core.cosmology import (fit_supernova_hubble_diagram, 
                               pure_temporal_magnitude,
                               calculate_hubble_parameter)
from udt.utils.data_loader import load_csp_database, load_pantheon_data
from udt.utils.plotting import plot_hubble_diagram, setup_plot_style


def analyze_supernova_sample(sn_data, dataset_name="Supernovae", data_type="unknown"):
    """Analyze a supernova dataset with UDT model."""
    print(f"\nAnalyzing {dataset_name} ({data_type})...")
    print(f"  Number of SNe: {len(sn_data)}")
    print(f"  Redshift range: {sn_data['z'].min():.4f} - {sn_data['z'].max():.4f}")
    
    # Prepare data
    z_obs = sn_data['z'].values
    m_obs = sn_data['m'].values
    m_err = sn_data['m_err'].values if 'm_err' in sn_data else np.ones_like(m_obs) * 0.15
    
    # Fit UDT model
    fit_result = fit_supernova_hubble_diagram(z_obs, m_obs, m_err, 
                                            R0_bounds=(100, 10000),
                                            M_B_bounds=(-22, -17))
    
    if fit_result['success']:
        print(f"  [OK] Fit successful")
        print(f"  R0 = {fit_result['R0']:.0f} Mpc")
        print(f"  M_B = {fit_result['M_B']:.2f}")
        print(f"  RMS = {fit_result['rms']:.3f} mag")
        
        chi2_dof = fit_result['chi2'] / (len(z_obs) - 2)
        print(f"  chi2/dof = {chi2_dof:.2f}")
        
        # Flag suspicious chi2/dof values
        if chi2_dof < 0.8:
            print(f"  [WARNING] chi2/dof = {chi2_dof:.2f} < 0.8 (possible data contamination)")
        elif chi2_dof > 5.0:
            print(f"  [WARNING] chi2/dof = {chi2_dof:.2f} > 5.0 (poor fit or underestimated errors)")
        else:
            print(f"  [OK] chi2/dof = {chi2_dof:.2f} (reasonable fit)")
        
        # Calculate effective Hubble parameter
        H_eff = calculate_hubble_parameter(fit_result['R0'])
        print(f"  H_eff = {H_eff:.1f} km/s/Mpc")
        
        # Add chi2/dof to result
        fit_result['chi2_dof'] = chi2_dof
        
    else:
        print(f"  [ERROR] Fit failed")
    
    return fit_result


def load_pantheon_raw_data(pantheon_file, z_max=0.1):
    """Load Pantheon+ data using RAW mB instead of corrected m_b_corr."""
    print(f"Loading Pantheon+ RAW data from: {pantheon_file}")
    
    # Load data
    pantheon_df = load_pantheon_data(pantheon_file)
    print(f"  Loaded {len(pantheon_df)} supernovae from Pantheon+")
    
    # Data contamination check
    available_columns = list(pantheon_df.columns)
    print(f"  Available magnitude columns: {[col for col in available_columns if 'm' in col.lower()]}")
    
    # Use RAW mB (SALT2 peak magnitude) instead of corrected m_b_corr
    if 'mB' in pantheon_df.columns:
        print("  [OK] Using RAW mB (SALT2 peak magnitude) - CLEAN DATA")
        pantheon_df['z'] = pantheon_df['zCMB']
        pantheon_df['m'] = pantheon_df['mB']  # RAW SALT2 magnitude
        pantheon_df['m_err'] = pantheon_df['mBERR']
        data_type = "RAW mB"
    else:
        print("  [WARNING] mB not found, falling back to m_b_corr - CONTAMINATED DATA")
        pantheon_df['z'] = pantheon_df['zCMB']
        pantheon_df['m'] = pantheon_df['m_b_corr']
        pantheon_df['m_err'] = pantheon_df['m_b_corr_err_DIAG']
        data_type = "CONTAMINATED m_b_corr"
    
    # Apply redshift cut
    pantheon_df = pantheon_df[pantheon_df['z'] <= z_max]
    print(f"  After z <= {z_max} cut: {len(pantheon_df)} supernovae")
    print(f"  Redshift range: {pantheon_df['z'].min():.4f} - {pantheon_df['z'].max():.4f}")
    
    return pantheon_df, data_type


def load_csp_raw_data(csp_dir, z_max=0.1):
    """Load CSP data using RAW B-band peak magnitudes."""
    print(f"Loading CSP RAW data from: {csp_dir}")
    
    # Load data
    csp_df = load_csp_database(csp_dir)
    print(f"  Loaded {len(csp_df)} supernovae from CSP DR3")
    
    # Use RAW B-band peak magnitudes
    print("  [OK] Using RAW B_peak_raw from individual light curves - CLEAN DATA")
    csp_df['z'] = csp_df['redshift']
    csp_df['m'] = csp_df['B_peak_raw']  # RAW B-band peak
    csp_df['m_err'] = csp_df['B_error_raw']
    data_type = "RAW B_peak_raw"
    
    # Apply redshift cut
    csp_df = csp_df[csp_df['z'] <= z_max]
    print(f"  After z <= {z_max} cut: {len(csp_df)} supernovae")
    print(f"  Redshift range: {csp_df['z'].min():.4f} - {csp_df['z'].max():.4f}")
    
    return csp_df, data_type


def compare_datasets(results):
    """Compare performance between datasets."""
    print("\n" + "=" * 80)
    print("DATASET COMPARISON - RAW DATA ONLY")
    print("=" * 80)
    
    if 'csp' in results and 'pantheon' in results:
        csp_result = results['csp']
        pantheon_result = results['pantheon']
        
        if csp_result['success'] and pantheon_result['success']:
            print(f"\nCSP DR3 (RAW B_peak_raw):")
            print(f"  R0 = {csp_result['R0']:.0f} Mpc")
            print(f"  RMS = {csp_result['rms']:.3f} mag")
            print(f"  chi2/dof = {csp_result['chi2_dof']:.2f}")
            
            print(f"\nPantheon+ (RAW mB):")
            print(f"  R0 = {pantheon_result['R0']:.0f} Mpc")
            print(f"  RMS = {pantheon_result['rms']:.3f} mag")
            print(f"  chi2/dof = {pantheon_result['chi2_dof']:.2f}")
            
            # Calculate differences
            delta_R0 = csp_result['R0'] - pantheon_result['R0']
            delta_RMS = csp_result['rms'] - pantheon_result['rms']
            delta_chi2 = csp_result['chi2_dof'] - pantheon_result['chi2_dof']
            
            print(f"\nCROSSCHECK ANALYSIS:")
            print(f"  Delta R0 = {delta_R0:+.0f} Mpc ({delta_R0/pantheon_result['R0']*100:+.1f}%)")
            print(f"  Delta RMS = {delta_RMS:+.3f} mag ({delta_RMS/pantheon_result['rms']*100:+.1f}%)")
            print(f"  Delta chi2/dof = {delta_chi2:+.2f}")
            
            # Assessment
            print(f"\nDATA QUALITY ASSESSMENT:")
            
            # Check for contamination indicators
            contamination_flags = []
            if pantheon_result['chi2_dof'] < 0.8:
                contamination_flags.append("Pantheon+ chi2/dof < 0.8 (suspiciously good)")
            if pantheon_result['rms'] < 0.15:
                contamination_flags.append("Pantheon+ RMS < 0.15 (unrealistically low)")
            if abs(delta_R0) / pantheon_result['R0'] > 0.5:
                contamination_flags.append("Large R0 difference (>50%)")
            
            if contamination_flags:
                print(f"  [WARNING] CONTAMINATION WARNINGS:")
                for flag in contamination_flags:
                    print(f"    - {flag}")
            else:
                print(f"  [OK] Both datasets show consistent results")
                
            # R0 consistency check
            if abs(delta_R0) / pantheon_result['R0'] < 0.2:
                print(f"  [OK] R0 values agree within 20% - good crosscheck")
            else:
                print(f"  [WARNING] R0 values differ by >20% - investigate dataset differences")


def main():
    parser = argparse.ArgumentParser(description='Analyze supernovae with UDT framework using RAW data only')
    parser.add_argument('--dataset', type=str, default='both',
                       choices=['csp', 'pantheon', 'both'],
                       help='Which dataset to analyze')
    parser.add_argument('--data-dir', type=str,
                       default='data',
                       help='Base data directory')
    parser.add_argument('--output-dir', type=str,
                       default='results/supernova_raw_analysis',
                       help='Directory for output plots and results')
    parser.add_argument('--z-max', type=float, default=0.1,
                       help='Maximum redshift to include')
    parser.add_argument('--plot', action='store_true',
                       help='Generate Hubble diagram plots')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("=" * 80)
    print("SUPERNOVA ANALYSIS - UDT FRAMEWORK (RAW DATA ONLY)")
    print("=" * 80)
    print(f"Theory: tau(r) = R0/(R0 + r)")
    print(f"Distance: d_L = z × R0")
    print(f"No expansion - pure temporal geometry")
    print(f"DATA CONTAMINATION PREVENTION: Using only RAW observational data")
    print("=" * 80)
    
    results = {}
    
    # Load and analyze CSP data
    if args.dataset in ['csp', 'both']:
        print(f"\n" + "=" * 50)
        print("CSP DR3 ANALYSIS")
        print("=" * 50)
        
        csp_dir = os.path.join(args.data_dir, 'CSP_Photometry_DR3', 'DR3')
        try:
            csp_df, csp_data_type = load_csp_raw_data(csp_dir, args.z_max)
            if len(csp_df) > 0:
                # Analyze
                csp_result = analyze_supernova_sample(csp_df, "CSP DR3", csp_data_type)
                results['csp'] = csp_result
                
                # Plot if requested
                if args.plot and csp_result['success']:
                    plot_path = os.path.join(args.output_dir, 'csp_raw_hubble_diagram.png')
                    plot_hubble_diagram(
                        csp_df['z'].values,
                        csp_df['m'].values,
                        csp_df['m_err'].values,
                        csp_result['m_predicted'],
                        title="CSP DR3 - Raw B-band Peak Magnitudes",
                        save_path=plot_path
                    )
                    print(f"  Plot saved to: {plot_path}")
                    
        except Exception as e:
            print(f"  Error loading CSP data: {e}")
    
    # Load and analyze Pantheon+ data
    if args.dataset in ['pantheon', 'both']:
        print(f"\n" + "=" * 50)
        print("PANTHEON+ ANALYSIS")
        print("=" * 50)
        
        pantheon_file = os.path.join(args.data_dir, 'Pantheon_SH0ES.dat')
        try:
            pantheon_df, pantheon_data_type = load_pantheon_raw_data(pantheon_file, args.z_max)
            if len(pantheon_df) > 0:
                # Analyze
                pantheon_result = analyze_supernova_sample(pantheon_df, "Pantheon+", pantheon_data_type)
                results['pantheon'] = pantheon_result
                
                # Plot if requested
                if args.plot and pantheon_result['success']:
                    plot_path = os.path.join(args.output_dir, 'pantheon_raw_hubble_diagram.png')
                    plot_hubble_diagram(
                        pantheon_df['z'].values,
                        pantheon_df['m'].values,
                        pantheon_df['m_err'].values,
                        pantheon_result['m_predicted'],
                        title="Pantheon+ - Raw SALT2 mB Magnitudes",
                        save_path=plot_path
                    )
                    print(f"  Plot saved to: {plot_path}")
                    
        except Exception as e:
            print(f"  Error loading Pantheon+ data: {e}")
    
    # Compare datasets
    if len(results) > 1:
        compare_datasets(results)
    
    # Save results
    if results:
        results_data = []
        for dataset_name, result in results.items():
            if result:
                results_data.append({
                    'dataset': dataset_name,
                    'data_type': 'RAW',
                    'R0_Mpc': result['R0'] if result['success'] else np.nan,
                    'M_B': result['M_B'] if result['success'] else np.nan,
                    'rms_mag': result['rms'] if result['success'] else np.nan,
                    'chi2_dof': result['chi2_dof'] if result['success'] else np.nan,
                    'success': result['success']
                })
        
        results_df = pd.DataFrame(results_data)
        results_file = os.path.join(args.output_dir, 'supernova_raw_results.csv')
        results_df.to_csv(results_file, index=False)
        print(f"\nResults saved to: {results_file}")
    
    print("\n" + "=" * 80)
    print("RAW DATA ANALYSIS COMPLETE")
    print("=" * 80)
    print("This analysis uses only raw observational data to avoid contamination:")
    print("- CSP: Raw B-band peak magnitudes from individual light curves")
    print("- Pantheon+: Raw SALT2 mB magnitudes (not corrected m_b_corr)")
    print("Chi2/dof values should be ~0.8-1.5 for good fits with clean data")
    print("=" * 80)


if __name__ == "__main__":
    main()