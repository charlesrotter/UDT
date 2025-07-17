#!/usr/bin/env python3
"""
Supernova Analysis using UDT Framework
=====================================

Analyzes Type Ia supernovae using the Universal Distance Dilation Theory.
This script replaces csp_udt_temporal.py with a cleaner structure.

Author: Charles Rotter
Date: 2025-01-17
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

# Add parent directory to path to import udt package
sys.path.insert(0, str(Path(__file__).parent.parent))

from udt.core.cosmology import (fit_supernova_hubble_diagram, 
                               pure_temporal_magnitude,
                               calculate_hubble_parameter,
                               compare_with_standard_cosmology)
from udt.utils.data_loader import load_csp_database, load_pantheon_data
from udt.utils.plotting import plot_hubble_diagram, setup_plot_style


def analyze_supernova_sample(sn_data, dataset_name="Supernovae"):
    """Analyze a supernova dataset with UDT model."""
    print(f"\nAnalyzing {dataset_name}...")
    print(f"  Number of SNe: {len(sn_data)}")
    print(f"  Redshift range: {sn_data['z'].min():.4f} - {sn_data['z'].max():.4f}")
    
    # Prepare data
    z_obs = sn_data['z'].values
    m_obs = sn_data['m'].values
    m_err = sn_data['m_err'].values if 'm_err' in sn_data else np.ones_like(m_obs) * 0.15
    
    # Fit UDT model
    # Use wider bounds to allow finding R0 ~3000 Mpc as in original analysis
    fit_result = fit_supernova_hubble_diagram(z_obs, m_obs, m_err, 
                                            R0_bounds=(100, 10000),
                                            M_B_bounds=(-22, -17))
    
    if fit_result['success']:
        print(f"  OK Fit successful")
        print(f"  R0 = {fit_result['R0']:.0f} Mpc")
        print(f"  M_B = {fit_result['M_B']:.2f}")
        print(f"  RMS = {fit_result['rms']:.3f} mag")
        print(f"  chi2/dof = {fit_result['chi2'] / (len(z_obs) - 2):.2f}")
        
        # Calculate effective Hubble parameter
        H_eff = calculate_hubble_parameter(fit_result['R0'])
        print(f"  H_eff = {H_eff:.1f} km/s/Mpc")
    else:
        print(f"  ERROR Fit failed")
    
    return fit_result


def main():
    parser = argparse.ArgumentParser(description='Analyze supernovae with UDT framework')
    parser.add_argument('--dataset', type=str, default='csp',
                       choices=['csp', 'pantheon', 'both'],
                       help='Which dataset to analyze')
    parser.add_argument('--data-dir', type=str,
                       default='data',
                       help='Base data directory')
    parser.add_argument('--output-dir', type=str,
                       default='results/supernova_analysis',
                       help='Directory for output plots and results')
    parser.add_argument('--z-max', type=float, default=0.1,
                       help='Maximum redshift to include')
    parser.add_argument('--plot', action='store_true',
                       help='Generate Hubble diagram plots')
    parser.add_argument('--compare', action='store_true',
                       help='Compare with standard cosmology')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("=" * 60)
    print("SUPERNOVA ANALYSIS - UDT FRAMEWORK")
    print("=" * 60)
    print(f"Theory: tau(r) = R0/(R0 + r)")
    print(f"Distance: d_L = z × R0")
    print(f"No expansion - pure temporal geometry")
    print("=" * 60)
    
    results = {}
    
    # Load and analyze CSP data
    if args.dataset in ['csp', 'both']:
        print(f"\nLoading CSP data...")
        csp_dir = os.path.join(args.data_dir, 'CSP_Photometry_DR3', 'DR3')
        try:
            csp_df = load_csp_database(csp_dir)
            if len(csp_df) > 0:
                # Prepare columns
                csp_df['z'] = csp_df['redshift']
                csp_df['m'] = csp_df['B_peak_raw']
                csp_df['m_err'] = csp_df['B_error_raw']
                
                # Apply redshift cut
                csp_df = csp_df[csp_df['z'] <= args.z_max]
                
                # Analyze
                csp_result = analyze_supernova_sample(csp_df, "CSP DR3")
                results['csp'] = csp_result
                
                # Plot if requested
                if args.plot and csp_result['success']:
                    plot_path = os.path.join(args.output_dir, 'csp_hubble_diagram.png')
                    plot_hubble_diagram(
                        csp_df['z'].values,
                        csp_df['m'].values,
                        csp_df['m_err'].values,
                        csp_result['m_predicted'],
                        title="CSP DR3 - UDT Hubble Diagram",
                        save_path=plot_path
                    )
                    print(f"  Plot saved to: {plot_path}")
                    
        except Exception as e:
            print(f"  Error loading CSP data: {e}")
    
    # Load and analyze Pantheon+ data
    if args.dataset in ['pantheon', 'both']:
        print(f"\nLoading Pantheon+ data...")
        pantheon_file = os.path.join(args.data_dir, 'Pantheon_SH0ES.dat')
        try:
            pantheon_df = load_pantheon_data(pantheon_file)
            if len(pantheon_df) > 0:
                # Prepare columns
                pantheon_df['z'] = pantheon_df['zCMB']
                pantheon_df['m'] = pantheon_df['m_b_corr']
                pantheon_df['m_err'] = pantheon_df['m_b_corr_err_DIAG']
                
                # Apply redshift cut
                pantheon_df = pantheon_df[pantheon_df['z'] <= args.z_max]
                
                # Analyze
                pantheon_result = analyze_supernova_sample(pantheon_df, "Pantheon+")
                results['pantheon'] = pantheon_result
                
                # Plot if requested
                if args.plot and pantheon_result['success']:
                    plot_path = os.path.join(args.output_dir, 'pantheon_hubble_diagram.png')
                    plot_hubble_diagram(
                        pantheon_df['z'].values,
                        pantheon_df['m'].values,
                        pantheon_df['m_err'].values,
                        pantheon_result['m_predicted'],
                        title="Pantheon+ - UDT Hubble Diagram",
                        save_path=plot_path
                    )
                    print(f"  Plot saved to: {plot_path}")
                    
        except Exception as e:
            print(f"  Error loading Pantheon+ data: {e}")
    
    # Summary
    print("\n" + "=" * 60)
    print("ANALYSIS SUMMARY")
    print("=" * 60)
    
    for dataset_name, result in results.items():
        if result and result['success']:
            print(f"\n{dataset_name.upper()}:")
            print(f"  R0 = {result['R0']:.0f} Mpc")
            print(f"  Scale ratio to galactic R0 (~38 kpc): {result['R0']*1000/38:.0f}:1")
            print(f"  RMS residual: {result['rms']:.3f} mag")
    
    # Comparison with standard cosmology
    if args.compare and any(r['success'] for r in results.values() if r):
        print("\n" + "=" * 60)
        print("COMPARISON WITH STANDARD COSMOLOGY")
        print("=" * 60)
        
        z_test = np.logspace(-3, np.log10(args.z_max), 100)
        
        for dataset_name, result in results.items():
            if result and result['success']:
                comparison = compare_with_standard_cosmology(
                    z_test, result['R0'], result['M_B']
                )
                
                print(f"\n{dataset_name.upper()} comparison:")
                print(f"  Mean Delta-m at z<0.05: {np.mean(comparison['delta_m'][z_test < 0.05]):.3f} mag")
                print(f"  Max |Delta-m| at z<{args.z_max}: {np.max(np.abs(comparison['delta_m'])):.3f} mag")
                
                if args.plot:
                    # Plot comparison
                    setup_plot_style()
                    import matplotlib.pyplot as plt
                    
                    fig, ax = plt.subplots(figsize=(10, 6))
                    ax.plot(z_test, comparison['delta_m'], 'b-', linewidth=2)
                    ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
                    ax.set_xlabel('Redshift z')
                    ax.set_ylabel('Delta-m (UDT - LCDM)')
                    ax.set_title(f'{dataset_name.upper()} - UDT vs Standard Cosmology')
                    ax.grid(True, alpha=0.3)
                    ax.set_xscale('log')
                    
                    comp_plot = os.path.join(args.output_dir, f'{dataset_name}_comparison.png')
                    plt.savefig(comp_plot, dpi=300, bbox_inches='tight')
                    print(f"  Comparison plot saved to: {comp_plot}")
    
    # Save results
    if results:
        results_data = []
        for dataset_name, result in results.items():
            if result:
                results_data.append({
                    'dataset': dataset_name,
                    'R0_Mpc': result['R0'] if result['success'] else np.nan,
                    'M_B': result['M_B'] if result['success'] else np.nan,
                    'rms_mag': result['rms'] if result['success'] else np.nan,
                    'chi2_dof': result['chi2'] / (len(result['residuals']) - 2) if result['success'] else np.nan,
                    'success': result['success']
                })
        
        results_df = pd.DataFrame(results_data)
        results_file = os.path.join(args.output_dir, 'supernova_udt_results.csv')
        results_df.to_csv(results_file, index=False)
        print(f"\nResults saved to: {results_file}")


if __name__ == "__main__":
    main()