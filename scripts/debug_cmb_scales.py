#!/usr/bin/env python3
"""
Debug CMB Scale Issues
======================

Systematic investigation of the scale mismatch in UDT CMB predictions.
Current issue: predicted l1 = 6 vs observed l1 = 220 (factor of ~37 error).

This script diagnoses:
1. Conformal time calculation problems
2. UDT-to-standard cosmology mapping issues  
3. Sound horizon integration accuracy
4. Angular scale formula validation

Author: Charles Rotter
Date: 2025-01-17
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import os
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from src.udt.diagnostics.parameter_registry import ParameterRegistry

class CMBScaleDebugger:
    """Debug CMB scale calculations in UDT framework."""
    
    def __init__(self):
        """Initialize with both UDT and standard cosmology parameters."""
        # Load validated parameters from registry
        self.registry = ParameterRegistry()
        
        # UDT parameters - CRITICAL: Use CMB scale, not supernova scale!
        cmb_params = self.registry.get_parameters_for_analysis('cmb')
        self.R0_cmb = cmb_params['R0_mpc']  # 13041.1 Mpc (CMB acoustic peak scale)
        
        # Standard cosmology (for comparison)
        self.c_light = 299792.458  # km/s
        self.H0_standard = 70.0    # km/s/Mpc
        self.Omega_m = 0.31        # Matter density
        self.Omega_b = 0.049       # Baryon density
        self.Omega_r = 8.24e-5     # Radiation density
        self.Omega_lambda = 0.69   # Dark energy
        
        # Known observational results (for validation)
        self.z_recombination = 1100.0
        self.eta_rec_standard = 288.0  # Mpc/c (standard cosmology)
        self.r_s_standard = 147.3      # Mpc (Planck 2018)
        self.D_A_standard = 13975.0    # Mpc (Planck 2018)
        self.l1_observed = 220.0       # First acoustic peak
        
    def debug_conformal_time_calculation(self):
        """Debug the conformal time calculation at recombination."""
        print("CONFORMAL TIME DEBUGGING")
        print("=" * 50)
        
        # Method 1: Current UDT approach
        z_rec = 1100.0
        eta_rec_udt_v1 = self.R0_cmb / (1 + z_rec)  # Mpc
        eta_rec_conf_v1 = eta_rec_udt_v1 / self.c_light  # Mpc/c
        
        print(f"Method 1 (current UDT):")
        print(f"  eta_rec = R0/(1+z) = {self.R0_cmb}/(1+{z_rec}) = {eta_rec_udt_v1:.2f} Mpc")
        print(f"  eta_rec_conf = {eta_rec_conf_v1:.6f} Mpc/c")
        print(f"  This is {eta_rec_conf_v1/self.eta_rec_standard*100:.2f}% of standard value")
        print()
        
        # Method 2: Scale UDT to match standard at recombination
        eta_rec_conf_v2 = self.eta_rec_standard  # Use standard value
        eta_rec_Mpc_v2 = eta_rec_conf_v2 * self.c_light
        
        print(f"Method 2 (calibrated to standard):")
        print(f"  Use standard eta_rec = {self.eta_rec_standard:.1f} Mpc/c")
        print(f"  Corresponds to eta_rec = {eta_rec_Mpc_v2:.1f} Mpc")
        print(f"  Implied UDT scale: R0_eff = {eta_rec_Mpc_v2 * (1 + z_rec):.1f} Mpc")
        print()
        
        # Method 3: Alternative UDT interpretation
        # Maybe UDT z-distance relation needs modification
        eta_rec_conf_v3 = 0.1  # Try reasonable intermediate value
        
        print(f"Method 3 (intermediate scale):")
        print(f"  Try eta_rec = {eta_rec_conf_v3:.3f} Mpc/c")
        print(f"  This is {eta_rec_conf_v3/self.eta_rec_standard*100:.2f}% of standard")
        print()
        
        return {
            'v1_current': eta_rec_conf_v1,
            'v2_standard': eta_rec_conf_v2, 
            'v3_intermediate': eta_rec_conf_v3
        }
    
    def debug_sound_horizon_integration(self, eta_rec_conf):
        """Debug sound horizon calculation."""
        print("SOUND HORIZON DEBUGGING")
        print("=" * 50)
        
        # Convert to Mpc for distance calculations
        eta_rec_Mpc = eta_rec_conf * self.c_light
        
        print(f"Using eta_rec = {eta_rec_conf:.6f} Mpc/c = {eta_rec_Mpc:.2f} Mpc")
        
        # Debug integration limits
        eta_early_options = [
            eta_rec_conf / 1000.0,   # Current method
            eta_rec_conf / 100.0,    # Less early
            eta_rec_conf / 10.0,     # Much less early
            0.001                    # Fixed early time
        ]
        
        print("\nIntegration limit sensitivity:")
        for i, eta_early in enumerate(eta_early_options):
            try:
                # Simple constant sound speed for testing
                c_s_const = 0.577  # c/sqrt(3) for radiation-dominated era
                
                # Direct integration
                r_s_simple = c_s_const * (eta_rec_conf - eta_early) * self.c_light
                
                print(f"  Option {i+1}: eta_early = {eta_early:.6f}, r_s = {r_s_simple:.1f} Mpc")
                
                # Ratio to standard
                ratio = r_s_simple / self.r_s_standard
                print(f"           Ratio to standard: {ratio:.3f}")
                
            except Exception as e:
                print(f"  Option {i+1}: Failed - {e}")
        
        print()
        
        # Check if sound speed evolution matters
        print("Sound speed evolution check:")
        
        # Method A: Constant sound speed
        c_s_const = 0.577
        r_s_constant = c_s_const * eta_rec_conf * self.c_light
        print(f"  Constant c_s = {c_s_const:.3f}: r_s = {r_s_constant:.1f} Mpc")
        
        # Method B: Evolving sound speed (simplified)
        def simple_sound_speed(eta):
            # Rough approximation: c_s decreases as baryons become important
            eta_Mpc = eta * self.c_light
            z_eff = max(1, self.R0_cmb / eta_Mpc - 1)
            R_gamma = 3 * self.Omega_b / (4 * self.Omega_r * (1 + z_eff))
            c_s_squared = 1.0 / (3.0 * (1.0 + R_gamma))
            return np.sqrt(c_s_squared)
        
        try:
            eta_early = eta_rec_conf / 100.0
            r_s_evolving, _ = quad(simple_sound_speed, eta_early, eta_rec_conf)
            r_s_evolving *= self.c_light
            print(f"  Evolving c_s: r_s = {r_s_evolving:.1f} Mpc")
        except:
            print(f"  Evolving c_s: Integration failed")
        
        return r_s_constant
    
    def debug_angular_scale_formula(self, eta_rec_conf, r_s):
        """Debug angular scale distance and peak position formula."""
        print("ANGULAR SCALE DEBUGGING")  
        print("=" * 50)
        
        eta_rec_Mpc = eta_rec_conf * self.c_light
        
        # Method 1: Current UDT approach
        tau_rec = self.R0_cmb / (self.R0_cmb + eta_rec_Mpc)
        D_A_udt = eta_rec_Mpc * tau_rec
        l1_udt = np.pi * D_A_udt / r_s
        
        print(f"Method 1 (UDT with temporal geometry):")
        print(f"  eta_rec = {eta_rec_Mpc:.2f} Mpc")
        print(f"  tau(eta_rec) = {tau_rec:.6f}")
        print(f"  D_A = eta_rec * tau = {D_A_udt:.2f} Mpc")
        print(f"  l1 = pi * D_A / r_s = {l1_udt:.1f}")
        print(f"  Error vs observed: {abs(l1_udt - self.l1_observed)/self.l1_observed*100:.1f}%")
        print()
        
        # Method 2: Standard flat cosmology (no UDT modification)
        D_A_flat = eta_rec_Mpc  # No temporal geometry
        l1_flat = np.pi * D_A_flat / r_s
        
        print(f"Method 2 (Flat space, no UDT):")
        print(f"  D_A = eta_rec = {D_A_flat:.2f} Mpc") 
        print(f"  l1 = pi * D_A / r_s = {l1_flat:.1f}")
        print(f"  Error vs observed: {abs(l1_flat - self.l1_observed)/self.l1_observed*100:.1f}%")
        print()
        
        # Method 3: Scale to match observations
        D_A_required = self.l1_observed * r_s / np.pi
        scale_factor = D_A_required / eta_rec_Mpc
        
        print(f"Method 3 (Scale to match l1 = 220):")
        print(f"  Required D_A = l1 * r_s / pi = {D_A_required:.1f} Mpc")
        print(f"  Required scale factor = {scale_factor:.2f}")
        print(f"  This suggests eta_rec should be {D_A_required:.1f} Mpc")
        print(f"  Or R0_eff should be {D_A_required * (1 + self.z_recombination):.1f} Mpc")
        print()
        
        return {
            'D_A_udt': D_A_udt,
            'D_A_flat': D_A_flat,
            'D_A_required': D_A_required,
            'l1_udt': l1_udt,
            'l1_flat': l1_flat,
            'scale_factor_needed': scale_factor
        }
    
    def compare_with_standard_cosmology(self):
        """Compare UDT approach with standard cosmology calculation."""
        print("STANDARD COSMOLOGY COMPARISON")
        print("=" * 50)
        
        # Standard cosmology calculation
        print(f"Standard cosmology (Planck 2018):")
        print(f"  Conformal time at recombination: {self.eta_rec_standard:.1f} Mpc/c")
        print(f"  Sound horizon: {self.r_s_standard:.1f} Mpc")
        print(f"  Angular diameter distance: {self.D_A_standard:.1f} Mpc")
        print(f"  First acoustic peak: l1 = {self.l1_observed:.1f}")
        print(f"  Acoustic scale: theta_s = r_s/D_A = {self.r_s_standard/self.D_A_standard:.6f} rad")
        print()
        
        # What UDT needs to match this
        print("UDT requirements to match standard:")
        
        # If we want l1 = 220, what do we need?
        r_s_target = self.r_s_standard  # Assume sound horizon is roughly correct
        D_A_target = self.l1_observed * r_s_target / np.pi
        eta_rec_target = D_A_target  # Assuming flat space approximation
        
        print(f"  Target angular distance: {D_A_target:.1f} Mpc")
        print(f"  Target conformal time: {eta_rec_target:.1f} Mpc")
        print(f"  Target conformal time: {eta_rec_target/self.c_light:.3f} Mpc/c")
        
        # What R0 would give this?
        R0_required = eta_rec_target * (1 + self.z_recombination)
        print(f"  Required UDT scale: R0 = {R0_required:.1f} Mpc")
        print(f"  Current R0: {self.R0_cmb:.1f} Mpc")
        print(f"  Scale factor needed: {R0_required/self.R0_cmb:.2f}")
        
        return {
            'R0_required': R0_required,
            'scale_factor': R0_required/self.R0_cmb
        }
    
    def create_scale_diagnostic_plots(self, output_dir="results/cmb_scale_debug"):
        """Create diagnostic plots for scale issues."""
        os.makedirs(output_dir, exist_ok=True)
        
        fig = plt.figure(figsize=(16, 12))
        
        # 1. Conformal time vs R0
        plt.subplot(2, 3, 1)
        R0_range = np.linspace(100, 10000, 100)
        eta_conf_range = R0_range / (1 + self.z_recombination) / self.c_light
        
        plt.plot(R0_range, eta_conf_range, 'b-', linewidth=2)
        plt.axhline(y=self.eta_rec_standard, color='red', linestyle='--', 
                   label=f'Standard = {self.eta_rec_standard:.0f} Mpc/c')
        plt.axvline(x=self.R0_cmb, color='green', linestyle='--', 
                   label=f'Current R0 = {self.R0_cmb:.0f} Mpc')
        plt.xlabel('R0 (Mpc)')
        plt.ylabel('Conformal Time (Mpc/c)')
        plt.title('Conformal Time at Recombination')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 2. Sound horizon scaling
        plt.subplot(2, 3, 2)
        r_s_range = 0.577 * eta_conf_range * self.c_light  # Constant c_s approximation
        
        plt.plot(R0_range, r_s_range, 'g-', linewidth=2)
        plt.axhline(y=self.r_s_standard, color='red', linestyle='--',
                   label=f'Standard = {self.r_s_standard:.0f} Mpc')
        plt.axvline(x=self.R0_cmb, color='green', linestyle='--')
        plt.xlabel('R0 (Mpc)')
        plt.ylabel('Sound Horizon (Mpc)')
        plt.title('Sound Horizon vs R0')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 3. Angular diameter distance
        plt.subplot(2, 3, 3)
        eta_Mpc_range = eta_conf_range * self.c_light
        tau_range = R0_range / (R0_range + eta_Mpc_range)
        D_A_range = eta_Mpc_range * tau_range
        
        plt.plot(R0_range, D_A_range, 'm-', linewidth=2, label='UDT')
        plt.plot(R0_range, eta_Mpc_range, 'c--', linewidth=2, label='Flat space')
        plt.axhline(y=self.l1_observed * self.r_s_standard / np.pi, 
                   color='red', linestyle='--', label='Required for l1=220')
        plt.axvline(x=self.R0_cmb, color='green', linestyle='--')
        plt.xlabel('R0 (Mpc)')
        plt.ylabel('Angular Distance (Mpc)')
        plt.title('Angular Diameter Distance vs R0')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 4. First acoustic peak position
        plt.subplot(2, 3, 4)
        l1_range = np.pi * D_A_range / r_s_range
        l1_flat_range = np.pi * eta_Mpc_range / r_s_range
        
        plt.plot(R0_range, l1_range, 'm-', linewidth=2, label='UDT')
        plt.plot(R0_range, l1_flat_range, 'c--', linewidth=2, label='Flat space')
        plt.axhline(y=self.l1_observed, color='red', linestyle='--',
                   label=f'Observed = {self.l1_observed:.0f}')
        plt.axvline(x=self.R0_cmb, color='green', linestyle='--')
        plt.xlabel('R0 (Mpc)')
        plt.ylabel('First Peak l1')
        plt.title('First Acoustic Peak vs R0')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.ylim(0, 500)
        
        # 5. Scale factor analysis
        plt.subplot(2, 3, 5)
        scale_factors = R0_range / self.R0_cmb
        l1_scaled = l1_range[0] * scale_factors  # Linear scaling approximation
        
        plt.plot(scale_factors, l1_scaled, 'orange', linewidth=2)
        plt.axhline(y=self.l1_observed, color='red', linestyle='--')
        plt.axvline(x=1.0, color='green', linestyle='--', label='Current')
        required_scale = self.l1_observed / l1_range[0]
        plt.axvline(x=required_scale, color='blue', linestyle='--', 
                   label=f'Required = {required_scale:.1f}')
        plt.xlabel('Scale Factor (R0/R0_current)')
        plt.ylabel('First Peak l1')
        plt.title('Peak Position vs Scale Factor')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 6. Error analysis
        plt.subplot(2, 3, 6)
        l1_error = abs(l1_range - self.l1_observed) / self.l1_observed * 100
        
        plt.plot(R0_range, l1_error, 'red', linewidth=2)
        plt.axvline(x=self.R0_cmb, color='green', linestyle='--')
        plt.axhline(y=5, color='orange', linestyle='--', label='5% error')
        plt.axhline(y=10, color='blue', linestyle='--', label='10% error')
        plt.xlabel('R0 (Mpc)')
        plt.ylabel('Error in l1 (%)')
        plt.title('Peak Position Error vs R0')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.ylim(0, 200)
        
        plt.tight_layout()
        plot_file = os.path.join(output_dir, 'cmb_scale_debugging.png')
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"Scale debugging plots saved to: {plot_file}")
        return plot_file
    
    def run_full_diagnosis(self):
        """Run complete diagnosis of CMB scale issues."""
        print("CMB SCALE ISSUE DIAGNOSIS")
        print("=" * 60)
        print()
        
        # 1. Debug conformal time
        eta_options = self.debug_conformal_time_calculation()
        print()
        
        # 2. Debug sound horizon for each conformal time option
        print("Testing sound horizon with different conformal times:")
        r_s_current = self.debug_sound_horizon_integration(eta_options['v1_current'])
        print()
        
        # 3. Debug angular scales
        scale_results = self.debug_angular_scale_formula(eta_options['v1_current'], r_s_current)
        print()
        
        # 4. Compare with standard cosmology
        standard_comparison = self.compare_with_standard_cosmology()
        print()
        
        # 5. Create diagnostic plots
        print("CREATING DIAGNOSTIC PLOTS")
        print("-" * 30)
        plot_file = self.create_scale_diagnostic_plots()
        print()
        
        # 6. Summary and recommendations
        print("DIAGNOSIS SUMMARY")
        print("=" * 50)
        print(f"Primary issue: Conformal time too small by factor of ~{self.eta_rec_standard/eta_options['v1_current']:.0f}")
        print(f"Current UDT gives: eta_rec = {eta_options['v1_current']:.6f} Mpc/c")
        print(f"Standard cosmology: eta_rec = {self.eta_rec_standard:.1f} Mpc/c")
        print()
        print(f"This leads to l1 = {scale_results['l1_udt']:.1f} vs observed l1 = {self.l1_observed:.1f}")
        print(f"Error: {abs(scale_results['l1_udt'] - self.l1_observed)/self.l1_observed*100:.1f}%")
        print()
        print("RECOMMENDED FIXES:")
        print("1. Recalibrate UDT R0 to match standard conformal time")
        print(f"   Current R0 = {self.R0_cmb:.1f} Mpc")
        print(f"   Required R0 = {standard_comparison['R0_required']:.1f} Mpc")
        print(f"   Scale factor = {standard_comparison['scale_factor']:.2f}")
        print()
        print("2. Alternative: Modify UDT z-distance relation")
        print("   Current: eta = R0/(1+z)")
        print("   Consider: eta = f(z) * R0/(1+z) with calibration factor")
        print()
        print("3. Check if UDT temporal geometry should affect CMB differently")
        print("   Current: D_A = eta * tau(eta)")
        print("   Consider: D_A = eta (flat space) for CMB")
        
        return {
            'eta_current': eta_options['v1_current'],
            'eta_required': self.eta_rec_standard,
            'scale_factor_needed': standard_comparison['scale_factor'],
            'l1_current': scale_results['l1_udt'],
            'l1_target': self.l1_observed
        }


def main():
    """Main debugging routine."""
    debugger = CMBScaleDebugger()
    results = debugger.run_full_diagnosis()
    return results


if __name__ == "__main__":
    main()