#!/usr/bin/env python
"""
Analyze Baryon Acoustic Oscillations (BAO) data with UDT predictions

This script compares UDT predictions against real BAO measurements from
multiple surveys (SDSS, BOSS, eBOSS, WiggleZ, DES, etc.)

BAO provides a standard ruler for cosmological distance measurements,
making it ideal for testing UDT's modified distance-redshift relations.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from datetime import datetime

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from udt.cosmology import UDTCosmology
from udt.constants import *

def load_bao_data(file_path):
    """Load BAO data from uncorBAO.txt file"""
    data = []
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or len(line.strip()) == 0:
                continue
            if 'File written by' in line or 'If you use' in line or 'Last update' in line:
                continue
            if 'zeff val error parameter' in line:
                continue
                
            parts = line.split()
            if len(parts) >= 6:
                try:
                    z = float(parts[0])
                    val = float(parts[1])
                    err = float(parts[2])
                    param = parts[3]
                    # Rest is experiment info
                    experiment = ' '.join(parts[6:]) if len(parts) > 6 else 'Unknown'
                    
                    data.append({
                        'z': z,
                        'value': val,
                        'error': err,
                        'parameter': param,
                        'experiment': experiment
                    })
                except ValueError:
                    continue
                    
    return data

def calculate_udt_bao_predictions(z_array, R0, cosmology):
    """Calculate UDT predictions for various BAO parameters"""
    predictions = {}
    
    # Calculate comoving distance
    D_c = cosmology.comoving_distance(z_array)
    
    # Angular diameter distance
    D_A = D_c / (1 + z_array)
    
    # Hubble parameter
    H_z = cosmology.hubble_parameter(z_array)
    
    # Sound horizon at drag epoch (fiducial value)
    r_d = 147.09  # Mpc (Planck 2018 value)
    
    # Calculate various BAO parameters
    predictions['DArd'] = D_A / r_d
    predictions['DHrd'] = c_LIGHT_SI / 1000 / H_z / r_d  # DH = c/H(z)
    predictions['DVrd'] = (z_array * D_A**2 * c_LIGHT_SI / 1000 / H_z)**(1/3) / r_d
    
    # For ratio parameters, we need reference values
    # These are survey-specific, so we'll handle them separately
    
    return predictions

def chi_squared_bao(R0, data, cosmology_class):
    """Calculate chi-squared for BAO data"""
    cosmology = cosmology_class(R0)
    chi2 = 0
    n_points = 0
    
    for point in data:
        z = point['z']
        obs_val = point['value']
        obs_err = point['error']
        param = point['parameter']
        
        # Calculate UDT prediction based on parameter type
        if param == 'DArd':
            D_A = cosmology.comoving_distance(z) / (1 + z)
            r_d = 147.09  # Mpc
            pred_val = D_A / r_d
        elif param == 'DHrd':
            H_z = cosmology.hubble_parameter(z)
            r_d = 147.09  # Mpc
            pred_val = c_LIGHT_SI / 1000 / H_z / r_d
        elif param == 'DVrd' or param == 'rdDV':
            D_A = cosmology.comoving_distance(z) / (1 + z)
            H_z = cosmology.hubble_parameter(z)
            D_V = (z * D_A**2 * c_LIGHT_SI / 1000 / H_z)**(1/3)
            r_d = 147.09  # Mpc
            if param == 'DVrd':
                pred_val = D_V / r_d
            else:  # rdDV
                pred_val = r_d / D_V
        elif param in ['DVratio', 'DAratio']:
            # These require reference redshifts - skip for now
            continue
        else:
            continue
            
        chi2 += ((obs_val - pred_val) / obs_err)**2
        n_points += 1
        
    return chi2, n_points

def analyze_bao_with_udt():
    """Main analysis function"""
    print("\n" + "="*70)
    print("BAO DATA ANALYSIS WITH UDT")
    print("="*70)
    
    # Load BAO data
    data_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'bao', 'uncorBAO.txt')
    if not os.path.exists(data_file):
        print(f"Error: BAO data file not found at {data_file}")
        return
        
    data = load_bao_data(data_file)
    print(f"\nLoaded {len(data)} BAO measurements")
    
    # Group by parameter type
    param_types = {}
    for point in data:
        param = point['parameter']
        if param not in param_types:
            param_types[param] = []
        param_types[param].append(point)
    
    print("\nParameter types found:")
    for param, points in param_types.items():
        print(f"  {param}: {len(points)} measurements")
    
    # Fit R0 using all compatible measurements
    print("\nFitting UDT parameter R0...")
    
    # Initial guess for R0
    R0_init = 14000  # Mpc
    
    # Minimize chi-squared
    result = minimize(lambda x: chi_squared_bao(x[0], data, UDTCosmology)[0], 
                     [R0_init], 
                     bounds=[(1000, 50000)],
                     method='L-BFGS-B')
    
    R0_best = result.x[0]
    chi2_best, n_points = chi_squared_bao(R0_best, data, UDTCosmology)
    chi2_dof = chi2_best / (n_points - 1) if n_points > 1 else np.inf
    
    print(f"\nBest-fit parameters:")
    print(f"  R0 = {R0_best:.1f} Mpc")
    print(f"  χ² = {chi2_best:.2f}")
    print(f"  n_points = {n_points}")
    print(f"  χ²/dof = {chi2_dof:.2f}")
    
    # Create plots for each parameter type
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    cosmology = UDTCosmology(R0_best)
    z_theory = np.linspace(0.01, 2.5, 100)
    
    # Plot different BAO parameters
    param_info = {
        'DArd': ('D_A/r_d', 0),
        'DHrd': ('D_H/r_d', 1),
        'DVrd': ('D_V/r_d', 2),
        'rdDV': ('r_d/D_V', 3)
    }
    
    for param, (label, idx) in param_info.items():
        if param not in param_types or idx >= len(axes):
            continue
            
        ax = axes[idx]
        
        # Plot data points
        points = param_types[param]
        z_data = [p['z'] for p in points]
        val_data = [p['value'] for p in points]
        err_data = [p['error'] for p in points]
        
        ax.errorbar(z_data, val_data, yerr=err_data, fmt='o', 
                   label='BAO Data', markersize=6, capsize=3)
        
        # Plot UDT prediction
        if param == 'DArd':
            D_A = cosmology.comoving_distance(z_theory) / (1 + z_theory)
            r_d = 147.09
            pred = D_A / r_d
        elif param == 'DHrd':
            H_z = cosmology.hubble_parameter(z_theory)
            r_d = 147.09
            pred = c_LIGHT_SI / 1000 / H_z / r_d
        elif param == 'DVrd':
            D_A = cosmology.comoving_distance(z_theory) / (1 + z_theory)
            H_z = cosmology.hubble_parameter(z_theory)
            D_V = (z_theory * D_A**2 * c_LIGHT_SI / 1000 / H_z)**(1/3)
            r_d = 147.09
            pred = D_V / r_d
        elif param == 'rdDV':
            D_A = cosmology.comoving_distance(z_theory) / (1 + z_theory)
            H_z = cosmology.hubble_parameter(z_theory)
            D_V = (z_theory * D_A**2 * c_LIGHT_SI / 1000 / H_z)**(1/3)
            r_d = 147.09
            pred = r_d / D_V
        else:
            continue
            
        ax.plot(z_theory, pred, 'r-', linewidth=2, 
                label=f'UDT (R₀={R0_best:.0f} Mpc)')
        
        ax.set_xlabel('Redshift z')
        ax.set_ylabel(label)
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.suptitle(f'UDT Fit to BAO Data (χ²/dof = {chi2_dof:.2f})', fontsize=14)
    plt.tight_layout()
    
    # Save results
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'results', 'bao_analysis')
    os.makedirs(output_dir, exist_ok=True)
    
    # Save plot
    plot_file = os.path.join(output_dir, f'bao_udt_fit_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png')
    plt.savefig(plot_file, dpi=150)
    print(f"\nPlot saved to: {plot_file}")
    
    # Save detailed results
    results_file = os.path.join(output_dir, f'bao_results_{datetime.now().strftime("%Y%m%d_%H%M%S")}.txt')
    with open(results_file, 'w') as f:
        f.write("BAO Analysis Results with UDT\n")
        f.write("="*50 + "\n\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Data file: {data_file}\n")
        f.write(f"Number of measurements: {len(data)}\n")
        f.write(f"Number used in fit: {n_points}\n\n")
        f.write("Best-fit parameters:\n")
        f.write(f"  R0 = {R0_best:.1f} Mpc\n")
        f.write(f"  χ² = {chi2_best:.2f}\n")
        f.write(f"  χ²/dof = {chi2_dof:.2f}\n\n")
        f.write("Data points by experiment:\n")
        
        # Group by experiment
        experiments = {}
        for point in data:
            exp = point['experiment']
            if exp not in experiments:
                experiments[exp] = []
            experiments[exp].append(point)
            
        for exp, points in experiments.items():
            f.write(f"\n{exp}:\n")
            for p in points:
                f.write(f"  z={p['z']:.3f}, {p['parameter']}={p['value']:.3f}±{p['error']:.3f}\n")
    
    print(f"Results saved to: {results_file}")
    
    plt.show()
    
    return R0_best, chi2_dof

if __name__ == "__main__":
    R0_best, chi2_dof = analyze_bao_with_udt()
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nSummary:")
    print(f"  Best-fit R0: {R0_best:.1f} Mpc")
    print(f"  χ²/dof: {chi2_dof:.2f}")
    
    if chi2_dof < 2.0:
        print("\n✓ Good fit to BAO data!")
    else:
        print("\n⚠ Moderate fit - may need additional parameters or corrections")