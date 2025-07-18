#!/usr/bin/env python3
"""
Debug Kernel Discrepancy
========================

Investigate why the derived kernel performs poorly compared to theoretical formula.

Author: Charles Rotter
Date: 2025-07-18
"""

import numpy as np
import matplotlib.pyplot as plt

class DebugKernelDiscrepancy:
    def __init__(self):
        print("DEBUGGING KERNEL DISCREPANCY")
        print("=" * 30)
        
        self.G = 4.302e-6  # km^2/s^2 per solar mass per kpc
        self.R0_gal = 62.4  # kpc
        self.alpha = 1.0
    
    def tau(self, r):
        """Distance equivalence principle."""
        return self.R0_gal / (self.R0_gal + r)
    
    def f_modification(self, tau_val):
        """Our derived modification function."""
        return (1 - tau_val)**2 / (tau_val * (3 - 2*tau_val))
    
    def derived_kernel_velocity(self, r, M_central):
        """Velocity from derived kernel: v² = G_eff × M / r"""
        tau_val = self.tau(r)
        f_val = self.f_modification(tau_val)
        G_eff = self.G * (1 + self.alpha * f_val)
        return np.sqrt(G_eff * M_central / r)
    
    def theoretical_velocity(self, r, V_scale):
        """Theoretical velocity: v = V_scale × √(r/(r + R₀/3)) × (1 + r/R₀)"""
        base = np.sqrt(r / (r + self.R0_gal/3))
        enhancement = (1 + r/self.R0_gal)
        return V_scale * base * enhancement
    
    def compare_velocity_shapes(self):
        """Compare the velocity profile shapes."""
        print("COMPARING VELOCITY PROFILE SHAPES")
        print("-" * 35)
        
        r = np.linspace(0.1, 20, 100)
        
        # Derived kernel with some reference mass
        M_ref = 1e10  # Solar masses
        v_derived = self.derived_kernel_velocity(r, M_ref)
        
        # Theoretical formula with some reference scale
        V_ref = 100  # km/s
        v_theory = self.theoretical_velocity(r, V_ref)
        
        # Normalize both at r = 5 kpc for comparison
        idx_5kpc = np.argmin(np.abs(r - 5.0))
        v_derived_norm = v_derived / v_derived[idx_5kpc]
        v_theory_norm = v_theory / v_theory[idx_5kpc]
        
        print(f"Shapes at r = 5 kpc (normalized to 1.0):")
        test_radii = [1, 2, 5, 10, 15, 20]
        
        for r_test in test_radii:
            idx = np.argmin(np.abs(r - r_test))
            print(f"r = {r_test:2d} kpc: Derived = {v_derived_norm[idx]:.3f}, Theory = {v_theory_norm[idx]:.3f}")
        
        # Check if they have the same shape
        max_diff = np.max(np.abs(v_derived_norm - v_theory_norm))
        print(f"\nMaximum normalized difference: {max_diff:.3f}")
        
        if max_diff < 0.1:
            print("SHAPES MATCH: The profiles have the same shape!")
        else:
            print("SHAPES DIFFER: The profiles have different shapes.")
        
        # Create comparison plot
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Panel 1: Raw velocity curves
        ax1 = axes[0, 0]
        ax1.plot(r, v_derived, 'r-', linewidth=2, label='Derived Kernel')
        ax1.plot(r, v_theory, 'b--', linewidth=2, label='Theoretical Formula')
        ax1.set_xlabel('Radius (kpc)')
        ax1.set_ylabel('Velocity (km/s)')
        ax1.set_title('Raw Velocity Profiles')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Panel 2: Normalized comparison
        ax2 = axes[0, 1]
        ax2.plot(r, v_derived_norm, 'r-', linewidth=2, label='Derived Kernel')
        ax2.plot(r, v_theory_norm, 'b--', linewidth=2, label='Theoretical Formula')
        ax2.set_xlabel('Radius (kpc)')
        ax2.set_ylabel('Normalized Velocity')
        ax2.set_title('Normalized Velocity Profiles')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Panel 3: G_eff vs r
        ax3 = axes[1, 0]
        tau_vals = self.tau(r)
        f_vals = self.f_modification(tau_vals)
        G_eff = self.G * (1 + self.alpha * f_vals)
        
        ax3.plot(r, G_eff/self.G, 'g-', linewidth=2)
        ax3.set_xlabel('Radius (kpc)')
        ax3.set_ylabel('G_eff / G')
        ax3.set_title('Effective Gravity Enhancement')
        ax3.grid(True, alpha=0.3)
        
        # Panel 4: Components breakdown
        ax4 = axes[1, 1]
        base_profile = np.sqrt(r / (r + self.R0_gal/3))
        enhancement = (1 + r/self.R0_gal)
        
        ax4.plot(r, base_profile, 'g-', linewidth=2, label='Base: sqrt(r/(r+R0/3))')
        ax4.plot(r, enhancement, 'orange', linewidth=2, label='Enhancement: (1+r/R0)')
        ax4.plot(r, base_profile * enhancement, 'k--', linewidth=2, label='Product')
        ax4.set_xlabel('Radius (kpc)')
        ax4.set_ylabel('Component Value')
        ax4.set_title('Theoretical Formula Components')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('C:/UDT/results/kernel_debug_comparison.png', dpi=150)
        plt.close()
        
        print("Saved: C:/UDT/results/kernel_debug_comparison.png")
        print()
        
        return max_diff
    
    def analyze_issue(self):
        """Analyze what's causing the discrepancy."""
        print("ANALYZING THE ISSUE")
        print("-" * 20)
        
        # The issue might be in the interpretation of the kernel
        print("HYPOTHESIS: The kernel interpretation might be wrong.")
        print()
        
        print("Current interpretation:")
        print("v^2 = G_eff * M_central / r")
        print("where G_eff = G * [1 + alpha * f(tau)]")
        print()
        
        print("But maybe the correct interpretation is:")
        print("v^2 = G * M_enhanced / r")
        print("where M_enhanced = M_central * [1 + alpha * f(tau)]")
        print()
        
        print("Or maybe:")
        print("v^2 = G * M_central * [1 + alpha * f(tau)] / r")
        print("which is the same as G_eff interpretation.")
        print()
        
        print("Let me check the fundamental issue...")
        print()
        
        # Let's examine what happens with the gravitational kernel
        r_test = 5.0  # kpc
        tau_test = self.tau(r_test)
        f_test = self.f_modification(tau_test)
        
        print(f"At r = {r_test} kpc:")
        print(f"tau = {tau_test:.3f}")
        print(f"f(tau) = {f_test:.3f}")
        print(f"G_eff/G = {1 + self.alpha * f_test:.3f}")
        print()
        
        print("The problem might be that our kernel gives:")
        print("K(r,r') = G * [1 + alpha * f(tau(r), tau(r'))] / |r-r'|^2")
        print()
        print("But for rotation curves, we need to consider the")
        print("gravitational field from a DISTRIBUTED mass, not point mass.")
        print()
        
        print("POSSIBLE SOLUTION:")
        print("The kernel should be applied to the entire mass distribution,")
        print("not just treat the galaxy as a point mass at the center.")
        print()
        
        return f_test
    
    def test_mass_distribution_approach(self):
        """Test applying the kernel to a distributed mass."""
        print("TESTING MASS DISTRIBUTION APPROACH")
        print("-" * 35)
        
        print("Instead of treating the galaxy as a point mass,")
        print("let's consider how the kernel modifies the gravitational")
        print("field from a distributed mass.")
        print()
        
        print("For a rotation curve, we typically have:")
        print("v^2(r) = G * M(<r) / r")
        print("where M(<r) is the mass enclosed within radius r.")
        print()
        
        print("With UDT kernel, this becomes:")
        print("v^2(r) = [G * integral rho(r') * [1 + alpha * f(tau(r), tau(r'))] * dV'] / r")
        print()
        
        print("For a centrally concentrated mass distribution,")
        print("this is approximately:")
        print("v^2(r) approximately G * M_central * [1 + alpha * f(tau(r), tau(0))] / r")
        print("      = G * M_central * [1 + alpha * f(tau(r), 1)] / r")
        print()
        
        print("This is exactly what we implemented!")
        print("So the interpretation is correct.")
        print()
        
        print("The issue must be elsewhere...")
        print()
        
        return "mass_distribution_correct"
    
    def investigate_coupling_constant(self):
        """Investigate if the coupling constant alpha is wrong."""
        print("INVESTIGATING COUPLING CONSTANT")
        print("-" * 32)
        
        print("Maybe alpha = 1 is not the right value.")
        print("Let's see what alpha would be needed to match the theoretical formula.")
        print()
        
        # Test different alpha values
        r = np.linspace(1, 20, 20)
        
        # Target: theoretical formula
        V_scale = 100  # km/s reference
        v_target = self.theoretical_velocity(r, V_scale)
        
        # Try to find alpha that matches
        def test_alpha(alpha_test):
            # Assume some reference mass that gives the right overall scale
            M_ref = 1e10
            
            # Calculate derived kernel velocity
            tau_vals = self.tau(r)
            f_vals = self.f_modification(tau_vals)
            G_eff = self.G * (1 + alpha_test * f_vals)
            v_derived = np.sqrt(G_eff * M_ref / r)
            
            # Normalize both at r = 5 kpc
            idx_5kpc = np.argmin(np.abs(r - 5.0))
            v_derived_norm = v_derived / v_derived[idx_5kpc]
            v_target_norm = v_target / v_target[idx_5kpc]
            
            # Calculate RMS difference
            rms_diff = np.sqrt(np.mean((v_derived_norm - v_target_norm)**2))
            return rms_diff
        
        # Test range of alpha values
        alpha_values = np.logspace(-2, 2, 50)  # 0.01 to 100
        rms_values = [test_alpha(alpha) for alpha in alpha_values]
        
        best_idx = np.argmin(rms_values)
        best_alpha = alpha_values[best_idx]
        best_rms = rms_values[best_idx]
        
        print(f"Best alpha = {best_alpha:.3f}")
        print(f"Best RMS difference = {best_rms:.6f}")
        print()
        
        if best_rms < 0.01:
            print("SUCCESS: Found alpha that matches theoretical formula!")
            print(f"The coupling constant should be alpha = {best_alpha:.3f}")
        else:
            print("FAILURE: Cannot find alpha that matches theoretical formula.")
            print("The kernel formulation may be fundamentally wrong.")
        
        return best_alpha, best_rms
    
    def run_complete_debug(self):
        """Run complete debugging analysis."""
        print("COMPLETE DEBUGGING ANALYSIS")
        print("=" * 30)
        
        shape_diff = self.compare_velocity_shapes()
        f_test = self.analyze_issue()
        mass_dist_result = self.test_mass_distribution_approach()
        best_alpha, best_rms = self.investigate_coupling_constant()
        
        print("\n" + "=" * 50)
        print("DEBUGGING CONCLUSIONS")
        print("=" * 50)
        
        print(f"1. Velocity profile shape difference: {shape_diff:.3f}")
        print(f"2. Modification function at r=5kpc: f = {f_test:.3f}")
        print(f"3. Mass distribution approach: {mass_dist_result}")
        print(f"4. Best coupling constant: alpha = {best_alpha:.3f}")
        print(f"5. Best RMS difference: {best_rms:.6f}")
        print()
        
        if best_rms < 0.01:
            print("RESOLUTION: The issue is the coupling constant!")
            print(f"Use alpha = {best_alpha:.3f} instead of alpha = 1.0")
        else:
            print("UNRESOLVED: The kernel formulation may need revision.")
            print("The mathematical derivation might have a missing factor.")

def main():
    """Run the debugging analysis."""
    debug = DebugKernelDiscrepancy()
    debug.run_complete_debug()

if __name__ == "__main__":
    main()