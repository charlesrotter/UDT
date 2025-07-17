#!/usr/bin/env python3
"""
Test R₀ Optimization Landscape for UDT Supernova Fit

Check if higher R₀ values could improve UDT performance.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
import os

def load_pantheon_raw():
    """Load raw Pantheon+ data"""
    pantheon_file = "data/Pantheon_SH0ES.dat"
    
    data = []
    with open(pantheon_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split()
            if len(parts) >= 12:
                try:
                    z_cmb = float(parts[4])      # z_CMB
                    m_b = float(parts[8])        # m_b_corr
                    m_b_err = float(parts[9])    # m_b_corr_err_DIAG
                    
                    if 0.01 < z_cmb < 2.3 and m_b_err < 1.0 and m_b > 0:
                        data.append([z_cmb, m_b, m_b_err])
                except (ValueError, IndexError):
                    continue
    
    data = np.array(data)
    return data[:, 0], data[:, 1], data[:, 2]  # z, m_b, m_b_err

def udt_distance(z, R0):
    """UDT luminosity distance"""
    return z * R0

def magnitude_from_distance(d_L, M_B):
    """Convert distance to magnitude"""
    mu = 5 * np.log10(d_L) + 25
    return M_B + mu

def chi_squared_for_R0(R0, z_data, m_b_data, m_b_err):
    """Calculate chi-squared for given R₀"""
    if R0 <= 0:
        return 1e10
    
    # For each R₀, find optimal M_B
    def chi2_for_MB(M_B):
        if M_B < -25 or M_B > -15:
            return 1e10
        
        d_L_udt = udt_distance(z_data, R0)
        m_b_udt = magnitude_from_distance(d_L_udt, M_B)
        return np.sum((m_b_data - m_b_udt)**2 / m_b_err**2)
    
    # Optimize M_B for this R₀
    result = minimize_scalar(chi2_for_MB, bounds=(-25, -15), method='bounded')
    return result.fun

def main():
    print("R_0 OPTIMIZATION LANDSCAPE ANALYSIS")
    print("=" * 50)
    print()
    
    # Load data
    print("Loading Pantheon+ raw data...")
    z_data, m_b_data, m_b_err = load_pantheon_raw()
    print(f"Loaded {len(z_data)} supernovae")
    print()
    
    # Test R_0 range
    print("Testing R_0 optimization landscape...")
    R0_values = np.logspace(2, 5, 50)  # 100 to 100,000 Mpc
    chi2_values = []
    
    for i, R0 in enumerate(R0_values):
        chi2 = chi_squared_for_R0(R0, z_data, m_b_data, m_b_err)
        chi2_values.append(chi2)
        
        if i % 10 == 0:
            print(f"  R_0 = {R0:.0f} Mpc, chi^2 = {chi2:.1f}")
    
    chi2_values = np.array(chi2_values)
    
    # Find minimum
    min_idx = np.argmin(chi2_values)
    best_R0 = R0_values[min_idx]
    best_chi2 = chi2_values[min_idx]
    
    print()
    print(f"OPTIMIZATION RESULTS:")
    print(f"  Best R_0 = {best_R0:.1f} Mpc")
    print(f"  Best chi^2 = {best_chi2:.1f}")
    print(f"  chi^2/DOF = {best_chi2/(len(z_data)-2):.2f}")
    print()
    
    # Check if higher R_0 improves fit
    high_R0_values = R0_values[R0_values > best_R0]
    high_chi2_values = chi2_values[R0_values > best_R0]
    
    if len(high_R0_values) > 0:
        print("HIGHER R_0 ANALYSIS:")
        for R0, chi2 in zip(high_R0_values[:5], high_chi2_values[:5]):
            improvement = (best_chi2 - chi2) / best_chi2 * 100
            print(f"  R_0 = {R0:.0f} Mpc: chi^2 = {chi2:.1f} ({improvement:+.1f}%)")
        print()
    
    # Create plot
    plt.figure(figsize=(10, 6))
    plt.semilogx(R0_values, chi2_values, 'b-', linewidth=2)
    plt.axvline(best_R0, color='red', linestyle='--', label=f'Best R_0 = {best_R0:.0f} Mpc')
    plt.axvline(3000, color='green', linestyle='--', label='Expected R_0 = 3000 Mpc')
    plt.xlabel('R_0 (Mpc)')
    plt.ylabel('chi^2')
    plt.title('UDT Supernova Fit: chi^2 vs R_0')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('results/udt_r0_optimization.png', dpi=300, bbox_inches='tight')
    print("Plot saved: results/udt_r0_optimization.png")
    
    # Theoretical implications
    print("THEORETICAL IMPLICATIONS:")
    
    if best_R0 > 10000:
        print("  High R_0 optimum suggests UDT needs modification")
        print("  Linear d_L = z * R_0 may be too simple")
    elif 2000 < best_R0 < 5000:
        print("  R_0 optimum is reasonable for cosmological scales")
        print("  But UDT still fails due to systematic trends")
    else:
        print("  Unexpected R_0 range - check analysis")
    
    print()
    print("MODIFICATION SUGGESTIONS:")
    print("1. Non-linear temporal geometry: tau(r) != R_0/(R_0 + r)")
    print("2. Redshift-dependent R_0: R_0(z)")
    print("3. Additional cosmological parameters")
    print("4. Hybrid UDT-expansion model")

if __name__ == "__main__":
    main()