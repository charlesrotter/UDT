#!/usr/bin/env python3
"""
Test UDT Cosmic Horizon Properties

Check if UDT leads to infinite mass at cosmic horizon and other singularities.

UDT metric: ds^2 = -c^2 * tau^2(r) * dt^2 + h(r) * dr^2 + r^2 * dOmega^2
where tau(r) = R_0 / (R_0 + r)
"""

import numpy as np
import matplotlib.pyplot as plt

def analyze_udt_horizon():
    """Analyze UDT cosmic horizon properties"""
    print("UDT COSMIC HORIZON ANALYSIS")
    print("=" * 50)
    print()
    
    # UDT parameters
    R0_cosmo = 3000  # Mpc (cosmological scale)
    c = 299792.458   # km/s (speed of light)
    
    print(f"UDT Parameters:")
    print(f"  R_0 = {R0_cosmo} Mpc")
    print(f"  Temporal function: tau(r) = R_0 / (R_0 + r)")
    print()
    
    # Distance range
    r_values = np.logspace(-2, 6, 1000)  # 0.01 to 1,000,000 Mpc
    
    # UDT temporal function
    tau_r = R0_cosmo / (R0_cosmo + r_values)
    
    # Effective light speed
    c_eff = c * tau_r
    
    # Distance modulus for UDT
    # In UDT: d_L = z * R_0, and for large distances z ≈ r/R_0
    z_approx = r_values / R0_cosmo
    d_L_udt = z_approx * R0_cosmo  # This equals r_values
    
    print("HORIZON ANALYSIS:")
    print()
    
    # Check for cosmic horizon
    print("1. COSMIC HORIZON LOCATION:")
    
    # In UDT, is there a distance where tau(r) → 0?
    # tau(r) = R_0/(R_0 + r) → 0 as r → ∞
    print("   tau(r) = R_0/(R_0 + r) -> 0 as r -> infinity")
    print("   No finite cosmic horizon - tau approaches zero asymptotically")
    print()
    
    # Check specific distances
    test_distances = [R0_cosmo, 10*R0_cosmo, 100*R0_cosmo, 1000*R0_cosmo]
    
    print("2. TAU(R) VALUES AT LARGE DISTANCES:")
    for r in test_distances:
        tau_val = R0_cosmo / (R0_cosmo + r)
        c_eff_val = c * tau_val
        print(f"   r = {r/R0_cosmo:.0f} × R_0: tau = {tau_val:.6f}, c_eff = {c_eff_val:.3f} km/s")
    print()
    
    # Mass analysis
    print("3. MASS IMPLICATIONS:")
    print()
    
    # In UDT, what happens to mass at large distances?
    # The key question: does the temporal geometry affect mass measurements?
    
    print("   UDT temporal geometry affects:")
    print("   - Light travel time: dt_obs = dt_emit / tau(r)")
    print("   - Frequency: f_obs = f_emit * tau(r)")
    print("   - Energy: E_obs = E_emit * tau(r)")
    print()
    
    # Energy/mass redshift factor
    redshift_factor = tau_r  # Energy scaling
    
    print("4. ENERGY/MASS REDSHIFT:")
    print("   Energy observed = Energy emitted × tau(r)")
    print("   As r -> infinity: tau(r) -> 0, so E_obs -> 0")
    print("   Objects become infinitely redshifted, not infinitely massive")
    print()
    
    # Check if there's a mass divergence
    print("5. MASS DIVERGENCE CHECK:")
    
    # In standard cosmology, mass density scales as rho ∝ (1+z)^3
    # In UDT, what's the equivalent?
    
    # For UDT: z ≈ r/R_0 for large r, so tau ≈ R_0/r
    # This means z ≈ 1/tau, so tau ≈ 1/(1+z) for large z
    
    # If UDT preserves mass-energy: rho(r) proportional to tau^(-3)?
    rho_scaling = tau_r**(-3)
    
    print("   If mass density scales as rho proportional to tau^(-3):")
    for i, r in enumerate(test_distances):
        tau_val = R0_cosmo / (R0_cosmo + r)
        rho_factor = tau_val**(-3)
        print(f"   r = {r/R0_cosmo:.0f} x R_0: rho_factor = {rho_factor:.2e}")
    
    print()
    print("   YES: Mass density diverges as tau -> 0 (r -> infinity)")
    print("   This suggests INFINITE MASS at cosmic distances")
    print()
    
    # Plot analysis
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: Temporal function
    ax1 = axes[0, 0]
    ax1.loglog(r_values/R0_cosmo, tau_r)
    ax1.set_xlabel('r/R_0')
    ax1.set_ylabel('tau(r)')
    ax1.set_title('UDT Temporal Function')
    ax1.grid(True, alpha=0.3)
    ax1.axhline(0.1, color='red', linestyle='--', alpha=0.5, label='tau = 0.1')
    ax1.axhline(0.01, color='red', linestyle='--', alpha=0.5, label='tau = 0.01')
    ax1.legend()
    
    # Plot 2: Effective light speed
    ax2 = axes[0, 1]
    ax2.loglog(r_values/R0_cosmo, c_eff)
    ax2.set_xlabel('r/R_0')
    ax2.set_ylabel('c_eff (km/s)')
    ax2.set_title('Effective Light Speed')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Mass density scaling
    ax3 = axes[1, 0]
    ax3.loglog(r_values/R0_cosmo, rho_scaling)
    ax3.set_xlabel('r/R_0')
    ax3.set_ylabel('rho/rho_0')
    ax3.set_title('Mass Density Scaling (proportional to tau^-3)')
    ax3.grid(True, alpha=0.3)
    ax3.axhline(1e6, color='red', linestyle='--', alpha=0.5, label='10^6 x rho_0')
    ax3.legend()
    
    # Plot 4: Redshift
    ax4 = axes[1, 1]
    z_udt = r_values / R0_cosmo
    ax4.loglog(r_values/R0_cosmo, z_udt)
    ax4.set_xlabel('r/R_0')
    ax4.set_ylabel('Redshift z')
    ax4.set_title('UDT Redshift vs Distance')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('results/udt_cosmic_horizon_analysis.png', dpi=300, bbox_inches='tight')
    print("Plot saved: results/udt_cosmic_horizon_analysis.png")
    print()
    
    # Theoretical conclusions
    print("THEORETICAL CONCLUSIONS:")
    print()
    print("1. NO FINITE COSMIC HORIZON:")
    print("   - tau(r) -> 0 asymptotically as r -> infinity")
    print("   - No sharp boundary, but effective horizon where tau approximately 0")
    print()
    print("2. INFINITE MASS PROBLEM:")
    print("   - If mass density proportional to tau^(-3), then rho -> infinity as r -> infinity")
    print("   - This violates energy conservation")
    print("   - Suggests UDT needs mass-energy modification")
    print()
    print("3. OBSERVATIONAL CONSEQUENCES:")
    print("   - Objects at large r become infinitely redshifted")
    print("   - But also infinitely massive (paradox)")
    print("   - May explain 'missing mass' as displaced to large r")
    print()
    print("4. POSSIBLE RESOLUTIONS:")
    print("   - Modify mass-energy conservation in UDT")
    print("   - Add cutoff at finite distance")
    print("   - Non-linear temporal geometry")
    print("   - Quantum effects at large scales")

def main():
    analyze_udt_horizon()

if __name__ == "__main__":
    main()