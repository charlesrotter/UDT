#!/usr/bin/env python3
"""
SPARC Galaxy Analysis using UDT Framework
========================================

Analyzes SPARC galaxy rotation curves using the Universal Distance Dilation Theory.
This script replaces temporal_unification_breakthrough.py with a cleaner structure.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

# Add parent directory to path to import udt package
sys.path.insert(0, str(Path(__file__).parent.parent))

from udt.core.galactic_dynamics import fit_galaxy_rotation_curve, pure_temporal_velocity
from udt.utils.data_loader import load_sparc_database
from udt.utils.plotting import plot_rotation_curve, setup_plot_style


def analyze_galaxy(galaxy_data, output_dir=None):
    """Analyze a single galaxy with UDT model."""
    name = galaxy_data['name']
    radius = galaxy_data['radius']
    velocity = galaxy_data['velocity']
    velocity_error = galaxy_data['velocity_error']
    
    print(f"\nAnalyzing {name}...")
    print(f"  Data points: {len(radius)}")
    print(f"  Radius range: {radius.min():.1f} - {radius.max():.1f} kpc")
    
    # Fit the model
    fit_result = fit_galaxy_rotation_curve(radius, velocity, velocity_error)
    
    if fit_result['success']:
        print(f"  OK Fit successful")
        print(f"  R0_gal = {fit_result['R0_gal']:.1f} kpc")
        print(f"  V_scale = {fit_result['V_scale']:.1f} km/s")
        print(f"  RMS = {fit_result['rms']:.1f} km/s")
        print(f"  chi2/dof = {fit_result['chi2'] / (len(radius) - 2):.2f}")
        
        # Plot if output directory specified
        if output_dir:
            plot_path = os.path.join(output_dir, f"{name}_fit.png")
            plot_rotation_curve(radius, velocity, velocity_error,
                              v_predicted=fit_result['v_predicted'],
                              galaxy_name=name, save_path=plot_path)
    else:
        print(f"  ERROR Fit failed")
    
    # Add metadata to result
    fit_result['name'] = name
    fit_result['n_points'] = len(radius)
    fit_result['r_max'] = radius.max()
    
    return fit_result


def main():
    parser = argparse.ArgumentParser(description='Analyze SPARC galaxies with UDT framework')
    parser.add_argument('--data-dir', type=str, 
                       default='data/sparc_database',
                       help='Path to SPARC data directory')
    parser.add_argument('--output-dir', type=str,
                       default='results/sparc_analysis',
                       help='Directory for output plots and results')
    parser.add_argument('--max-galaxies', type=int, default=None,
                       help='Maximum number of galaxies to analyze')
    parser.add_argument('--plot', action='store_true',
                       help='Generate plots for each galaxy')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("=" * 60)
    print("SPARC GALAXY ANALYSIS - UDT FRAMEWORK")
    print("=" * 60)
    print(f"Theory: tau(r) = R0/(R0 + r)")
    print(f"Enhancement: 1/tau^2 = (1 + r/R0)^2")
    print("=" * 60)
    
    # Load galaxy data
    print(f"\nLoading data from: {args.data_dir}")
    try:
        galaxies = load_sparc_database(args.data_dir)
    except Exception as e:
        print(f"Error loading data: {e}")
        return
    
    # Limit number if requested
    if args.max_galaxies:
        galaxies = galaxies[:args.max_galaxies]
    
    print(f"\nAnalyzing {len(galaxies)} galaxies...")
    
    # Analyze each galaxy
    results = []
    successful_fits = 0
    
    for galaxy_data in galaxies:
        result = analyze_galaxy(
            galaxy_data, 
            output_dir=args.output_dir if args.plot else None
        )
        results.append(result)
        if result['success']:
            successful_fits += 1
    
    # Summary statistics
    print("\n" + "=" * 60)
    print("ANALYSIS SUMMARY")
    print("=" * 60)
    print(f"Total galaxies analyzed: {len(results)}")
    if len(results) > 0:
        print(f"Successful fits: {successful_fits} ({successful_fits/len(results)*100:.1f}%)")
    else:
        print("No galaxies were successfully analyzed")
    
    # Calculate statistics for successful fits
    successful_results = [r for r in results if r['success']]
    if successful_results:
        R0_values = [r['R0_gal'] for r in successful_results]
        rms_values = [r['rms'] for r in successful_results]
        
        print(f"\nR0_gal statistics:")
        print(f"  Mean: {np.mean(R0_values):.1f} kpc")
        print(f"  Median: {np.median(R0_values):.1f} kpc")
        print(f"  Std: {np.std(R0_values):.1f} kpc")
        
        print(f"\nRMS residual statistics:")
        print(f"  Mean: {np.mean(rms_values):.1f} km/s")
        print(f"  Median: {np.median(rms_values):.1f} km/s")
        print(f"  Std: {np.std(rms_values):.1f} km/s")
    
    # Save results
    results_df = pd.DataFrame(results)
    results_file = os.path.join(args.output_dir, 'sparc_udt_results.csv')
    results_df.to_csv(results_file, index=False)
    print(f"\nResults saved to: {results_file}")
    
    # Create summary plot if requested
    if args.plot and successful_results:
        setup_plot_style()
        import matplotlib.pyplot as plt
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # R0 distribution
        ax1.hist(R0_values, bins=20, alpha=0.7, color='blue', edgecolor='black')
        ax1.axvline(np.median(R0_values), color='red', linestyle='--', 
                   label=f'Median = {np.median(R0_values):.1f} kpc')
        ax1.set_xlabel('R₀_gal (kpc)')
        ax1.set_ylabel('Number of Galaxies')
        ax1.set_title('Distribution of Galactic Scale Parameter')
        ax1.legend()
        
        # RMS distribution
        ax2.hist(rms_values, bins=20, alpha=0.7, color='green', edgecolor='black')
        ax2.axvline(np.median(rms_values), color='red', linestyle='--',
                   label=f'Median = {np.median(rms_values):.1f} km/s')
        ax2.set_xlabel('RMS Residual (km/s)')
        ax2.set_ylabel('Number of Galaxies')
        ax2.set_title('Distribution of Fit Quality')
        ax2.legend()
        
        plt.tight_layout()
        summary_plot = os.path.join(args.output_dir, 'sparc_summary.png')
        plt.savefig(summary_plot, dpi=300, bbox_inches='tight')
        print(f"Summary plot saved to: {summary_plot}")


if __name__ == "__main__":
    main()