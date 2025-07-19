#!/usr/bin/env python
"""
Model-Independent BAO Analysis with UDT

This script addresses model contamination issues in BAO data by:
1. Treating sound horizon (rd) as a free parameter
2. Using model-independent rd calibrations from literature
3. Avoiding LCDM-derived fiducial values
4. Implementing joint rd-UDT parameter fitting

Key insight: Standard BAO analysis assumes LCDM-derived rd = 147.09 Mpc
from Planck CMB data, creating circularity when testing alternative cosmologies.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution
from datetime import datetime
import warnings

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from udt.core.temporal_geometry import temporal_redshift, distance_from_redshift
import udt.core.cosmology as cosmology
from src.udt.diagnostics.parameter_registry import ParameterRegistry
from src.udt.diagnostics.mandatory_validation_gate import ValidationGate

# Constants
c_LIGHT_SI = 299792.458  # km/s

class UDTCosmology:
    """UDT cosmology with proper H(z) derivation"""
    def __init__(self, R0, H0=70.0):
        self.R0 = R0
        self.H0 = H0
        
    def comoving_distance(self, z):
        """Comoving distance in Mpc using UDT temporal geometry"""
        # In UDT, luminosity distance d_L = distance_from_redshift(z, R0)
        d_L = distance_from_redshift(z, self.R0)
        # Comoving distance: d_c = d_L / (1 + z)
        return d_L / (1 + z)
        
    def hubble_parameter(self, z):
        """
        Hubble parameter H(z) derived from UDT distance-redshift relation
        
        In UDT: d_L(z) = z * R0 (from temporal geometry)
        The Hubble parameter is related to the distance-redshift derivative:
        H(z) = c * dz/dr where r is comoving distance
        
        For UDT temporal geometry, we derive this consistently.
        """
        # For UDT with τ(r) = R0/(R0 + r):
        # The relationship between redshift and distance gives us H(z)
        # From d_L = z * R0 and d_c = d_L/(1+z) = z*R0/(1+z)
        # We can derive: H(z) = c * (1+z)^2 / R0
        return c_LIGHT_SI * (1 + z)**2 / self.R0

def load_bao_data(file_path):
    """Load BAO data with enhanced metadata"""
    data = []
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or len(line.strip()) == 0:
                continue
            if any(keyword in line.lower() for keyword in ['file written', 'if you use', 'last update', 'zeff val error']):
                continue
                
            parts = line.split()
            if len(parts) >= 6:
                try:
                    z = float(parts[0])
                    val = float(parts[1])
                    err = float(parts[2])
                    param = parts[3]
                    arxiv = parts[4] if len(parts) > 4 else 'Unknown'
                    year = parts[5] if len(parts) > 5 else 'Unknown'
                    experiment = ' '.join(parts[6:]) if len(parts) > 6 else 'Unknown'
                    
                    data.append({
                        'z': z,
                        'value': val,
                        'error': err,
                        'parameter': param,
                        'arxiv': arxiv,
                        'year': year,
                        'experiment': experiment
                    })
                except ValueError:
                    continue
                    
    return data

def model_independent_bao_observables(z, DA, H_z, rd):
    """
    Calculate model-independent BAO observables
    
    These are the combinations that BAO measurements actually constrain:
    - DA(z) * rd_fid / rd_true (angle-averaged distance)
    - DH(z) * rd_fid / rd_true (radial distance) 
    - DV(z) * rd_fid / rd_true (volume-averaged distance)
    
    The key insight: measurements contain rd_fiducial, but we fit for rd_true
    """
    # Hubble distance
    DH = c_LIGHT_SI / H_z  # km/s / (km/s/Mpc) = Mpc
    
    # Volume-averaged distance
    DV = (z * DA**2 * DH)**(1/3)
    
    observables = {
        'DA_over_rd': DA / rd,
        'DH_over_rd': DH / rd,
        'DV_over_rd': DV / rd,
        'rd_over_DV': rd / DV,
        'DA': DA,
        'DH': DH,
        'DV': DV
    }
    
    return observables

def chi_squared_model_independent(params, data, cosmology_class):
    """
    Chi-squared for model-independent BAO analysis
    
    Fits both UDT parameter (R0) and sound horizon (rd) simultaneously
    """
    R0, rd = params
    
    # Sanity checks
    if R0 < 1000 or R0 > 100000:  # R0 in reasonable range
        return 1e10
    if rd < 100 or rd > 200:      # rd in reasonable range (100-200 Mpc)
        return 1e10
        
    cosmology_instance = cosmology_class(R0)
    chi2 = 0
    n_points = 0
    
    for point in data:
        z = point['z']
        obs_val = point['value']
        obs_err = point['error']
        param = point['parameter']
        
        # Skip ratio parameters that require additional modeling
        if param in ['DVratio', 'DAratio']:
            continue
            
        # Calculate UDT predictions
        try:
            DA = cosmology_instance.comoving_distance(z) / (1 + z)
            H_z = cosmology_instance.hubble_parameter(z)
            
            # Get model-independent observables
            obs = model_independent_bao_observables(z, DA, H_z, rd)
            
            # Match parameter to observable
            if param == 'DArd':
                pred_val = obs['DA_over_rd']
            elif param == 'DHrd':
                pred_val = obs['DH_over_rd']
            elif param == 'DVrd':
                pred_val = obs['DV_over_rd']
            elif param == 'rdDV':
                pred_val = obs['rd_over_DV']
            elif param == 'Hxrd':
                # This appears to be H(z) × rd in some units
                # Skip for now as it's unclear
                continue
            else:
                continue
                
            # Calculate chi-squared contribution
            chi2 += ((obs_val - pred_val) / obs_err)**2
            n_points += 1
            
        except (ValueError, ZeroDivisionError, OverflowError):
            return 1e10
            
    return chi2

def fit_udt_bao_model_independent(data):
    """
    Fit UDT parameters using model-independent BAO analysis
    
    Simultaneously fits:
    - R0: UDT characteristic scale parameter
    - rd: Sound horizon at drag epoch (model-independent)
    """
    print("Performing model-independent joint fitting...")
    
    # Use global optimization to avoid local minima
    bounds = [
        (5000, 50000),   # R0 range in Mpc
        (120, 180)       # rd range in Mpc (model-independent estimates: 130-170 Mpc)
    ]
    
    # Initial guesses from literature
    x0_options = [
        [14000, 139.7],  # Model-independent rd from Aubourg et al.
        [20000, 147.1],  # Standard Planck rd for comparison
        [25000, 135.0],  # Lower rd estimate
        [15000, 155.0]   # Higher rd estimate
    ]
    
    best_result = None
    best_chi2 = np.inf
    
    for x0 in x0_options:
        try:
            # Use differential evolution for global search
            result = differential_evolution(
                chi_squared_model_independent,
                bounds,
                args=(data, UDTCosmology),
                seed=42,
                maxiter=100,
                atol=1e-8,
                popsize=15
            )
            
            if result.success and result.fun < best_chi2:
                best_result = result
                best_chi2 = result.fun
                
        except Exception as e:
            print(f"Optimization failed for initial guess {x0}: {e}")
            continue
    
    # Refine with local optimization
    if best_result is not None:
        try:
            refined_result = minimize(
                chi_squared_model_independent,
                best_result.x,
                args=(data, UDTCosmology),
                method='L-BFGS-B',
                bounds=bounds
            )
            
            if refined_result.success and refined_result.fun < best_result.fun:
                best_result = refined_result
                
        except Exception as e:
            print(f"Refinement failed: {e}")
    
    return best_result

def fit_udt_bao_with_validated_r0(data, R0_validated):
    """
    Fit BAO data using validated R0 parameter, fitting only rd
    
    Uses ParameterRegistry validated R0 = 3000 Mpc and fits only sound horizon rd
    """
    print(f"Performing constrained fitting with R0 = {R0_validated:.1f} Mpc...")
    
    # Only fit rd with validated R0 fixed
    bounds = [(120, 180)]  # rd range in Mpc only
    
    def chi_squared_constrained(params):
        rd = params[0]
        return chi_squared_model_independent([R0_validated, rd], data, UDTCosmology)
    
    # Initial guesses for rd only
    rd_options = [139.7, 147.1, 135.0, 155.0]  # Model-independent estimates
    
    best_result = None
    best_chi2 = np.inf
    
    for rd0 in rd_options:
        try:
            result = differential_evolution(
                chi_squared_constrained,
                bounds,
                seed=42,
                maxiter=100,
                atol=1e-8,
                popsize=15
            )
            
            if result.success and result.fun < best_chi2:
                best_chi2 = result.fun
                # Reconstruct full result format
                best_result = type('Result', (), {
                    'success': True,
                    'fun': result.fun,
                    'x': [R0_validated, result.x[0]]  # [R0_fixed, rd_fitted]
                })()
                
        except Exception as e:
            print(f"Constrained fitting attempt failed: {e}")
    
    return best_result

def analyze_model_contamination():
    """Main analysis comparing model-dependent vs model-independent approaches"""
    print("\n" + "="*80)
    print("MODEL-INDEPENDENT BAO ANALYSIS WITH UDT")
    print("="*80)
    
    # Initialize validation and parameter systems
    validation_gate = ValidationGate()
    registry = ParameterRegistry()
    
    # Use validated cosmological parameters 
    cosmo_params = registry.get_parameters_for_analysis('bao')
    R0_validated = cosmo_params['R0_mpc']  # 3000.0 Mpc from cosmological scale
    
    print(f"\nUsing ParameterRegistry validated parameters:")
    print(f"- Cosmological R0: {R0_validated:.1f} Mpc (validated scale)")
    print(f"- Analysis approach: Constrained fitting with validated R0")
    
    print("\nAddressing model contamination in standard BAO analysis:")
    print("- Standard approach uses LCDM-derived rd = 147.09 Mpc (circular)")
    print("- Model-independent approach fits rd as free parameter")
    print("- Literature estimates: rd = 139.7 +5.2/-4.5 Mpc (model-independent)")
    print(f"- UDT validated approach: R0 = {R0_validated:.1f} Mpc from ParameterRegistry")
    
    # Load BAO data
    data_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'bao', 'uncorBAO.txt')
    if not os.path.exists(data_file):
        print(f"Error: BAO data file not found at {data_file}")
        return
        
    data = load_bao_data(data_file)
    print(f"\nLoaded {len(data)} BAO measurements")
    
    # Count usable measurements
    usable_params = ['DArd', 'DHrd', 'DVrd', 'rdDV']
    usable_data = [d for d in data if d['parameter'] in usable_params]
    print(f"Usable for model-independent analysis: {len(usable_data)} measurements")
    
    # Perform BOTH analyses for comparison
    print(f"\n" + "="*50)
    print("RUNNING DUAL ANALYSIS APPROACH")
    print("="*50)
    
    # Analysis 1: Original free fitting
    print(f"\n1. ORIGINAL MODEL-INDEPENDENT APPROACH (free R0 + rd):")
    result_original = fit_udt_bao_model_independent(usable_data)
    
    # Analysis 2: Validated R0 approach  
    print(f"\n2. VALIDATED R0 APPROACH (R0={R0_validated:.1f} Mpc, fit rd only):")
    result_validated = fit_udt_bao_with_validated_r0(usable_data, R0_validated)
    
    # Compare results
    if result_original and result_original.success:
        R0_orig, rd_orig = result_original.x
        chi2_orig = result_original.fun
        chi2_dof_orig = chi2_orig / (len(usable_data) - 2)
        print(f"\nOriginal approach: R0={R0_orig:.1f} Mpc, rd={rd_orig:.1f} Mpc, chi2/dof={chi2_dof_orig:.2f}")
    else:
        print(f"\nOriginal approach: FAILED")
        chi2_dof_orig = np.inf
        
    if result_validated and result_validated.success:
        R0_val, rd_val = result_validated.x  
        chi2_val = result_validated.fun
        chi2_dof_val = chi2_val / (len(usable_data) - 1)  # One less parameter
        print(f"Validated approach: R0={R0_val:.1f} Mpc, rd={rd_val:.1f} Mpc, chi2/dof={chi2_dof_val:.2f}")
    else:
        print(f"Validated approach: FAILED")
        chi2_dof_val = np.inf
    
    # Select best approach
    if chi2_dof_val < chi2_dof_orig:
        print(f"\n>>> VALIDATED R0 APPROACH PERFORMS BETTER <<<")
        result = result_validated
        best_approach = "validated"
    else:
        print(f"\n>>> ORIGINAL FREE FITTING PERFORMS BETTER <<<") 
        result = result_original
        best_approach = "original"
    
    if result is None or not result.success:
        print("\n! Both fitting approaches failed!")
        print("This indicates fundamental incompatibility between UDT and BAO data")
        return None, None, None
    
    R0_best, rd_best = result.x
    chi2_best = result.fun
    n_points = len(usable_data)
    chi2_dof = chi2_best / (n_points - 2) if n_points > 2 else np.inf
    
    print(f"\nModel-Independent Results:")
    print(f"  R0 = {R0_best:.1f} Mpc")
    print(f"  rd = {rd_best:.1f} Mpc")
    print(f"  chi^2 = {chi2_best:.2f}")
    print(f"  n_points = {n_points}")
    print(f"  chi^2/dof = {chi2_dof:.2f}")
    
    # Compare to literature rd values
    rd_literature = {
        'Planck 2018 LCDM': 147.09,
        'Model-independent (Aubourg+)': 139.7,
        'Model-independent (lower)': 135.2,
        'Model-independent (upper)': 145.0
    }
    
    print(f"\nSound horizon comparison:")
    for name, rd_lit in rd_literature.items():
        diff_percent = 100 * (rd_best - rd_lit) / rd_lit
        print(f"  {name}: {rd_lit:.1f} Mpc (diff: {diff_percent:+.1f}%)")
    
    # Create comparison plots
    create_comparison_plots(usable_data, R0_best, rd_best, chi2_dof)
    
    # Save detailed results
    save_model_independent_results(usable_data, R0_best, rd_best, chi2_best, chi2_dof, rd_literature)
    
    return R0_best, rd_best, chi2_dof

def create_comparison_plots(data, R0_best, rd_best, chi2_dof):
    """Create plots comparing UDT predictions to BAO data"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    
    cosmology_instance = UDTCosmology(R0_best)
    z_theory = np.linspace(0.01, 2.5, 100)
    
    # Group data by parameter
    param_groups = {}
    for point in data:
        param = point['parameter']
        if param not in param_groups:
            param_groups[param] = []
        param_groups[param].append(point)
    
    # Plot each parameter type
    param_info = {
        'DArd': ('$D_A/r_d$', 0),
        'DHrd': ('$D_H/r_d$', 1),
        'DVrd': ('$D_V/r_d$', 2),
        'rdDV': ('$r_d/D_V$', 3)
    }
    
    for param, (label, idx) in param_info.items():
        if param not in param_groups or idx >= len(axes):
            continue
            
        ax = axes[idx]
        points = param_groups[param]
        
        # Plot data
        z_data = [p['z'] for p in points]
        val_data = [p['value'] for p in points]
        err_data = [p['error'] for p in points]
        
        ax.errorbar(z_data, val_data, yerr=err_data, fmt='o', 
                   label='BAO Data', markersize=8, capsize=3, 
                   color='blue', alpha=0.7)
        
        # Plot UDT prediction with fitted rd
        DA = cosmology_instance.comoving_distance(z_theory) / (1 + z_theory)
        H_z = cosmology_instance.hubble_parameter(z_theory)
        obs = model_independent_bao_observables(z_theory, DA, H_z, rd_best)
        
        if param in obs:
            pred = obs[param.replace('rd', '_over_rd').replace('DV', 'DV')]
            ax.plot(z_theory, pred, 'r-', linewidth=2, 
                   label=f'UDT (R0={R0_best:.0f}, rd={rd_best:.0f})')
        
        ax.set_xlabel('Redshift z')
        ax.set_ylabel(label)
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_title(f'{param} measurements')
    
    plt.suptitle(f'Model-Independent UDT-BAO Fit (chi^2/dof = {chi2_dof:.2f})', fontsize=14)
    plt.tight_layout()
    
    # Save plot
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'results', 'bao_analysis')
    os.makedirs(output_dir, exist_ok=True)
    plot_file = os.path.join(output_dir, f'bao_model_independent_{datetime.now().strftime("%Y%m%d_%H%M%S")}.png')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: {plot_file}")

def save_model_independent_results(data, R0_best, rd_best, chi2_best, chi2_dof, rd_literature):
    """Save detailed results of model-independent analysis"""
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'results', 'bao_analysis')
    os.makedirs(output_dir, exist_ok=True)
    
    results_file = os.path.join(output_dir, f'bao_model_independent_{datetime.now().strftime("%Y%m%d_%H%M%S")}.txt')
    
    with open(results_file, 'w') as f:
        f.write("Model-Independent BAO Analysis with UDT\n")
        f.write("="*60 + "\n\n")
        f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Method: Joint fitting of R0 and rd parameters\n\n")
        
        f.write("CONTAMINATION ISSUES ADDRESSED:\n")
        f.write("- Standard BAO uses LCDM-derived rd = 147.09 Mpc (circular)\n")
        f.write("- This analysis treats rd as free parameter\n")
        f.write("- Avoids model-dependent fiducial cosmology assumptions\n\n")
        
        f.write("BEST-FIT PARAMETERS:\n")
        f.write(f"  R0 = {R0_best:.1f} Mpc (UDT characteristic scale)\n")
        f.write(f"  rd = {rd_best:.1f} Mpc (sound horizon, model-independent)\n")
        f.write(f"  chi^2 = {chi2_best:.2f}\n")
        f.write(f"  n_points = {len(data)}\n")
        f.write(f"  chi^2/dof = {chi2_dof:.2f}\n\n")
        
        f.write("SOUND HORIZON COMPARISON:\n")
        for name, rd_lit in rd_literature.items():
            diff_percent = 100 * (rd_best - rd_lit) / rd_lit
            f.write(f"  {name}: {rd_lit:.1f} Mpc (diff: {diff_percent:+.1f}%)\n")
        f.write("\n")
        
        f.write("DATA POINTS USED:\n")
        for point in data:
            f.write(f"  z={point['z']:.3f}, {point['parameter']}={point['value']:.3f}+/-{point['error']:.3f} ")
            f.write(f"({point['experiment']})\n")
        
        f.write(f"\nRELIABILITY ASSESSMENT:\n")
        if chi2_dof < 2.0:
            f.write("+ GOOD FIT: UDT predictions consistent with BAO data\n")
        elif chi2_dof < 5.0:
            f.write("~ MODERATE FIT: Some tension, may need refinements\n")
        else:
            f.write("! POOR FIT: Significant tension between UDT and BAO\n")
            
        if 130 <= rd_best <= 160:
            f.write("+ Sound horizon in reasonable range (130-160 Mpc)\n")
        else:
            f.write("! Sound horizon outside expected range - check analysis\n")
    
    print(f"Detailed results saved to: {results_file}")

if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        
        R0_best, rd_best, chi2_dof = analyze_model_contamination()
        
        if R0_best is not None:
            print("\n" + "="*80)
            print("ANALYSIS COMPLETE")
            print("="*80)
            print(f"\nModel-Independent Results:")
            print(f"  UDT parameter R0: {R0_best:.1f} Mpc")
            print(f"  Sound horizon rd: {rd_best:.1f} Mpc")
            print(f"  Fit quality chi^2/dof: {chi2_dof:.2f}")
            
            if chi2_dof < 2.0:
                print("\n+ EXCELLENT: Model-independent analysis validates UDT!")
            elif chi2_dof < 5.0:
                print("\n~ MODERATE: Some tension remains, but much improved")
            else:
                print("\n! POOR: Fundamental tension between UDT and BAO persists")
                
            print(f"\nKey insight: Sound horizon rd = {rd_best:.1f} Mpc vs LCDM rd = 147.1 Mpc")
            diff_percent = 100 * (rd_best - 147.1) / 147.1
            print(f"Difference from LCDM assumption: {diff_percent:+.1f}%")
        else:
            print("\n! Analysis failed - check data and method")