#!/usr/bin/env python3
"""
Calibrate CMB R0 Parameter
===========================

Precise calibration of R0_cmb to match observed CMB acoustic peak positions.
Uses iterative optimization to find R0_cmb that gives l1 = 220.

The challenge is that both sound horizon and angular diameter distance 
depend on R0_cmb, so we need to solve:
l1 = pi * D_A(R0_cmb) / r_s(R0_cmb) = 220

Author: Charles Rotter
Date: 2025-01-17
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import minimize_scalar, fsolve
import os

class CMBCalibrator:
    """Calibrate R0_cmb parameter for correct CMB acoustic peak positions."""
    
    def __init__(self):
        """Initialize with physical constants and target values."""
        
        # Physical constants
        self.c_light = 299792.458       # km/s
        self.z_recombination = 1100.0   # Recombination redshift
        
        # Standard cosmology parameters
        self.Omega_b = 0.049            # Baryon density
        self.Omega_r = 8.24e-5          # Radiation density
        
        # Target observational values
        self.l1_target = 220.0          # First acoustic peak (Planck)
        self.eta_rec_standard = 288.0   # Standard conformal time (Mpc/c)
        self.r_s_standard = 147.3       # Standard sound horizon (Mpc)
        self.D_A_standard = 13975.0     # Standard angular diameter distance (Mpc)
        
        # Current best guess
        self.R0_cmb_guess = 3000.0 * 3786  # From debug analysis
        
        print("CMB R0 CALIBRATOR")
        print("=" * 40)
        print(f"Target first peak: l1 = {self.l1_target}")
        print(f"Standard values:")
        print(f"  eta_rec = {self.eta_rec_standard:.1f} Mpc/c")
        print(f"  r_s = {self.r_s_standard:.1f} Mpc")
        print(f"  D_A = {self.D_A_standard:.1f} Mpc")
        print(f"Initial R0_cmb guess: {self.R0_cmb_guess:.0f} Mpc")
        print()
    
    def conformal_time_at_recombination(self, R0_cmb):
        """Calculate conformal time at recombination for given R0_cmb."""
        eta_rec_Mpc = R0_cmb / (1 + self.z_recombination)
        eta_rec_conf = eta_rec_Mpc / self.c_light
        return eta_rec_conf
    
    def sound_horizon(self, R0_cmb):
        """Calculate sound horizon for given R0_cmb."""
        eta_rec = self.conformal_time_at_recombination(R0_cmb)
        
        # Simplified sound speed integration
        # Use constant sound speed for initial calibration
        c_s_avg = 0.577  # c/sqrt(3) approximation
        
        # Integration from early times to recombination
        eta_early = eta_rec / 1000.0
        r_s = c_s_avg * (eta_rec - eta_early) * self.c_light
        
        return r_s
    
    def angular_diameter_distance(self, R0_cmb):
        """Calculate angular diameter distance for given R0_cmb."""
        eta_rec = self.conformal_time_at_recombination(R0_cmb)
        eta_rec_Mpc = eta_rec * self.c_light
        
        # UDT temporal geometry
        tau_rec = R0_cmb / (R0_cmb + eta_rec_Mpc)
        D_A = eta_rec_Mpc * tau_rec
        
        return D_A
    
    def first_acoustic_peak(self, R0_cmb):
        """Calculate first acoustic peak position for given R0_cmb."""
        r_s = self.sound_horizon(R0_cmb)
        D_A = self.angular_diameter_distance(R0_cmb)
        
        l1 = np.pi * D_A / r_s
        return l1
    
    def objective_function(self, R0_cmb):
        """Objective function for optimization: minimize |l1 - target|."""
        l1_pred = self.first_acoustic_peak(R0_cmb)
        error = abs(l1_pred - self.l1_target)
        return error
    
    def calibrate_R0_cmb(self):
        """Find optimal R0_cmb that gives l1 = 220."""
        print("CALIBRATING R0_CMB")
        print("=" * 30)
        
        # Test initial guess
        l1_initial = self.first_acoustic_peak(self.R0_cmb_guess)
        print(f"Initial guess: R0_cmb = {self.R0_cmb_guess:.0f} Mpc")
        print(f"  gives l1 = {l1_initial:.1f} (target = {self.l1_target})")
        print(f"  error = {abs(l1_initial - self.l1_target):.1f}")
        print()
        
        # Optimization bounds (reasonable range)
        R0_min = 1e5   # 100,000 Mpc
        R0_max = 1e9   # 1 billion Mpc
        
        print(f"Optimizing R0_cmb in range [{R0_min:.0e}, {R0_max:.0e}] Mpc")
        
        # Find optimal R0_cmb
        result = minimize_scalar(self.objective_function, 
                               bounds=(R0_min, R0_max), 
                               method='bounded')
        
        R0_optimal = result.x
        l1_optimal = self.first_acoustic_peak(R0_optimal)
        error_optimal = result.fun
        
        print(f"\nOptimization result:")
        print(f"  R0_cmb_optimal = {R0_optimal:.0f} Mpc ({R0_optimal/1e6:.1f} Gpc)")
        print(f"  l1_predicted = {l1_optimal:.1f}")
        print(f"  error = {error_optimal:.2f}")
        print(f"  success = {result.success}")
        
        # Validate with detailed calculations
        print(f"\nDetailed validation:")
        eta_rec = self.conformal_time_at_recombination(R0_optimal)
        r_s = self.sound_horizon(R0_optimal)
        D_A = self.angular_diameter_distance(R0_optimal)
        
        print(f"  eta_rec = {eta_rec:.3f} Mpc/c (standard = {self.eta_rec_standard:.1f})")
        print(f"  r_s = {r_s:.1f} Mpc (standard = {self.r_s_standard:.1f})")
        print(f"  D_A = {D_A:.1f} Mpc (standard = {self.D_A_standard:.1f})")
        print(f"  acoustic scale = {r_s/D_A:.6f} rad")
        
        # Check scale ratios
        eta_ratio = eta_rec / self.eta_rec_standard
        r_s_ratio = r_s / self.r_s_standard
        D_A_ratio = D_A / self.D_A_standard
        
        print(f"\nRatios to standard cosmology:")
        print(f"  eta_rec ratio = {eta_ratio:.3f}")
        print(f"  r_s ratio = {r_s_ratio:.3f}")
        print(f"  D_A ratio = {D_A_ratio:.3f}")
        
        return {
            'R0_cmb_optimal': R0_optimal,
            'l1_predicted': l1_optimal,
            'error': error_optimal,
            'eta_rec': eta_rec,
            'r_s': r_s,
            'D_A': D_A,
            'success': result.success
        }
    
    def alternative_calibration_methods(self):
        """Try alternative calibration approaches."""
        print("\nALTERNATIVE CALIBRATION METHODS")
        print("=" * 40)
        
        # Method 1: Match conformal time directly
        print("Method 1: Match conformal time to standard")
        R0_eta_match = self.eta_rec_standard * self.c_light * (1 + self.z_recombination)
        l1_eta_match = self.first_acoustic_peak(R0_eta_match)
        
        print(f"  R0_cmb = {R0_eta_match:.0f} Mpc")
        print(f"  l1 = {l1_eta_match:.1f}")
        print(f"  error = {abs(l1_eta_match - self.l1_target):.1f}")
        print()
        
        # Method 2: Match sound horizon directly
        print("Method 2: Iterative sound horizon matching")
        
        def sound_horizon_error(R0_cmb):
            r_s_pred = self.sound_horizon(R0_cmb)
            return abs(r_s_pred - self.r_s_standard)
        
        result_rs = minimize_scalar(sound_horizon_error, 
                                  bounds=(1e5, 1e9), 
                                  method='bounded')
        
        R0_rs_match = result_rs.x
        l1_rs_match = self.first_acoustic_peak(R0_rs_match)
        
        print(f"  R0_cmb = {R0_rs_match:.0f} Mpc")
        print(f"  l1 = {l1_rs_match:.1f}")
        print(f"  error = {abs(l1_rs_match - self.l1_target):.1f}")
        print()
        
        # Method 3: Scale factor approach
        print("Method 3: Scale factor from debug analysis")
        
        # From debug analysis, we know we need scale factor ~ 40
        scale_factor_needed = 220.0 / 5.9  # Rough estimate
        R0_scale = self.R0_cmb_guess * scale_factor_needed
        l1_scale = self.first_acoustic_peak(R0_scale)
        
        print(f"  Scale factor = {scale_factor_needed:.1f}")
        print(f"  R0_cmb = {R0_scale:.0f} Mpc")
        print(f"  l1 = {l1_scale:.1f}")
        print(f"  error = {abs(l1_scale - self.l1_target):.1f}")
        
        return {
            'eta_match': {'R0': R0_eta_match, 'l1': l1_eta_match},
            'rs_match': {'R0': R0_rs_match, 'l1': l1_rs_match},
            'scale_match': {'R0': R0_scale, 'l1': l1_scale}
        }
    
    def create_calibration_plots(self, calibration_result, output_dir="results/cmb_calibration"):
        """Create diagnostic plots for calibration process."""
        os.makedirs(output_dir, exist_ok=True)
        
        fig = plt.figure(figsize=(16, 12))
        
        # 1. R0_cmb vs l1 curve
        plt.subplot(2, 3, 1)
        R0_range = np.logspace(5, 9, 100)
        l1_range = [self.first_acoustic_peak(R0) for R0 in R0_range]
        
        plt.semilogx(R0_range, l1_range, 'b-', linewidth=2)
        plt.axhline(y=self.l1_target, color='red', linestyle='--', 
                   label=f'Target l1 = {self.l1_target}')
        plt.axvline(x=calibration_result['R0_cmb_optimal'], 
                   color='green', linestyle='--', 
                   label=f'Optimal R0 = {calibration_result["R0_cmb_optimal"]:.0e}')
        
        plt.xlabel('R0_cmb (Mpc)')
        plt.ylabel('First Peak l1')
        plt.title('R0_cmb Calibration Curve')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.ylim(0, 500)
        
        # 2. Sound horizon vs R0_cmb
        plt.subplot(2, 3, 2)
        r_s_range = [self.sound_horizon(R0) for R0 in R0_range]
        
        plt.semilogx(R0_range, r_s_range, 'g-', linewidth=2)
        plt.axhline(y=self.r_s_standard, color='red', linestyle='--',
                   label=f'Standard r_s = {self.r_s_standard:.1f} Mpc')
        plt.axvline(x=calibration_result['R0_cmb_optimal'], 
                   color='green', linestyle='--')
        
        plt.xlabel('R0_cmb (Mpc)')
        plt.ylabel('Sound Horizon (Mpc)')
        plt.title('Sound Horizon vs R0_cmb')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 3. Angular diameter distance vs R0_cmb
        plt.subplot(2, 3, 3)
        D_A_range = [self.angular_diameter_distance(R0) for R0 in R0_range]
        
        plt.semilogx(R0_range, D_A_range, 'm-', linewidth=2)
        plt.axhline(y=self.D_A_standard, color='red', linestyle='--',
                   label=f'Standard D_A = {self.D_A_standard:.0f} Mpc')
        plt.axvline(x=calibration_result['R0_cmb_optimal'], 
                   color='green', linestyle='--')
        
        plt.xlabel('R0_cmb (Mpc)')
        plt.ylabel('Angular Diameter Distance (Mpc)')
        plt.title('Angular Distance vs R0_cmb')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 4. Conformal time vs R0_cmb
        plt.subplot(2, 3, 4)
        eta_range = [self.conformal_time_at_recombination(R0) for R0 in R0_range]
        
        plt.semilogx(R0_range, eta_range, 'orange', linewidth=2)
        plt.axhline(y=self.eta_rec_standard, color='red', linestyle='--',
                   label=f'Standard eta_rec = {self.eta_rec_standard:.1f} Mpc/c')
        plt.axvline(x=calibration_result['R0_cmb_optimal'], 
                   color='green', linestyle='--')
        
        plt.xlabel('R0_cmb (Mpc)')
        plt.ylabel('Conformal Time (Mpc/c)')
        plt.title('Conformal Time vs R0_cmb')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 5. Optimization convergence (zoom)
        plt.subplot(2, 3, 5)
        R0_optimal = calibration_result['R0_cmb_optimal']
        R0_zoom = np.linspace(R0_optimal * 0.8, R0_optimal * 1.2, 100)
        l1_zoom = [self.first_acoustic_peak(R0) for R0 in R0_zoom]
        error_zoom = [abs(l1 - self.l1_target) for l1 in l1_zoom]
        
        plt.plot(R0_zoom, error_zoom, 'r-', linewidth=2)
        plt.axvline(x=R0_optimal, color='green', linestyle='--',
                   label=f'Minimum at {R0_optimal:.0f}')
        
        plt.xlabel('R0_cmb (Mpc)')
        plt.ylabel('Error |l1 - target|')
        plt.title('Optimization Convergence')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 6. Scale comparison
        plt.subplot(2, 3, 6)
        scales = ['Standard', 'UDT Optimized']
        eta_values = [self.eta_rec_standard, calibration_result['eta_rec']]
        r_s_values = [self.r_s_standard, calibration_result['r_s']]
        D_A_values = [self.D_A_standard/1000, calibration_result['D_A']/1000]  # In Gpc
        
        x = np.arange(len(scales))
        width = 0.25
        
        plt.bar(x - width, eta_values, width, label='eta_rec (Mpc/c)', alpha=0.7)
        plt.bar(x, r_s_values, width, label='r_s (Mpc)', alpha=0.7)
        plt.bar(x + width, D_A_values, width, label='D_A (Gpc)', alpha=0.7)
        
        plt.xlabel('Model')
        plt.ylabel('Value')
        plt.title('Scale Comparison')
        plt.xticks(x, scales)
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plot_file = os.path.join(output_dir, 'cmb_r0_calibration.png')
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"Calibration plots saved to: {plot_file}")
        return plot_file


def main():
    """Main calibration routine."""
    calibrator = CMBCalibrator()
    
    # Primary calibration
    result = calibrator.calibrate_R0_cmb()
    
    # Alternative methods
    alternatives = calibrator.alternative_calibration_methods()
    
    # Create diagnostic plots
    plot_file = calibrator.create_calibration_plots(result)
    
    # Summary
    print("\nCALIBRATION SUMMARY")
    print("=" * 40)
    print(f"Optimal R0_cmb: {result['R0_cmb_optimal']:.0f} Mpc ({result['R0_cmb_optimal']/1e6:.1f} Gpc)")
    print(f"Predicted l1: {result['l1_predicted']:.1f}")
    print(f"Target l1: {calibrator.l1_target}")
    print(f"Error: {result['error']:.2f}")
    print(f"Success: {result['success']}")
    
    if result['error'] < 1.0:
        print("+ Calibration successful!")
    elif result['error'] < 5.0:
        print("+ Calibration good (within 5)")
    else:
        print("! Calibration needs improvement")
    
    return result


if __name__ == "__main__":
    main()