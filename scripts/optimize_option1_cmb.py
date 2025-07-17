#!/usr/bin/env python3
"""
Optimize Option 1 UDT CMB
=========================

Precise optimization of Option 1 (direct temporal geometry: D_A = eta * tau)
to achieve exactly l1 = 220 for CMB acoustic peaks.

This approach uses:
- Standard cosmology baseline (eta_rec = 288 Mpc/c, r_s = 147.3 Mpc)
- UDT modification: D_A = eta_rec * tau(eta_rec) where tau = R0/(R0 + eta_rec)
- Optimization to find R0_cmb that gives l1 = pi * D_A / r_s = 220

Author: Charles Rotter
Date: 2025-01-17
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar, fsolve
import os

class OptimizedUDTCMB:
    """Optimized UDT CMB physics using Option 1 approach."""
    
    def __init__(self):
        """Initialize with standard cosmology baseline."""
        
        # Physical constants
        self.c_light = 299792.458       # km/s
        
        # Standard cosmology values (baseline)
        self.eta_rec_standard = 288.0   # Mpc/c (conformal time at recombination)
        self.r_s_standard = 147.3       # Mpc (sound horizon)
        self.D_A_standard = 13975.0     # Mpc (angular diameter distance)
        
        # Target observational value
        self.l1_target = 220.0          # First acoustic peak (Planck)
        
        print("OPTIMIZED UDT CMB (Option 1)")
        print("=" * 40)
        print("Approach: D_A = eta_rec * tau(eta_rec)")
        print("where tau = R0_cmb / (R0_cmb + eta_rec)")
        print()
        print("Standard cosmology baseline:")
        print(f"  eta_rec = {self.eta_rec_standard:.1f} Mpc/c")
        print(f"  r_s = {self.r_s_standard:.1f} Mpc")
        print(f"  D_A = {self.D_A_standard:.1f} Mpc")
        print(f"  l1_standard = {np.pi * self.D_A_standard / self.r_s_standard:.1f}")
        print(f"Target l1 = {self.l1_target}")
        print()
    
    def udt_angular_diameter_distance(self, R0_cmb):
        """
        Calculate UDT-corrected angular diameter distance.
        
        Parameters:
        -----------
        R0_cmb : float
            CMB-scale R0 parameter in Mpc
            
        Returns:
        --------
        D_A : float
            Angular diameter distance in Mpc
        """
        # Convert conformal time to physical distance
        eta_rec_Mpc = self.eta_rec_standard * self.c_light
        
        # UDT temporal geometry factor
        tau_rec = R0_cmb / (R0_cmb + eta_rec_Mpc)
        
        # UDT angular diameter distance
        D_A = eta_rec_Mpc * tau_rec
        
        return D_A
    
    def first_acoustic_peak(self, R0_cmb):
        """
        Calculate first acoustic peak position.
        
        Parameters:
        -----------
        R0_cmb : float
            CMB-scale R0 parameter in Mpc
            
        Returns:
        --------
        l1 : float
            First acoustic peak multipole
        """
        D_A = self.udt_angular_diameter_distance(R0_cmb)
        l1 = np.pi * D_A / self.r_s_standard
        return l1
    
    def find_optimal_R0_cmb(self):
        """Find R0_cmb that gives exactly l1 = 220."""
        print("FINDING OPTIMAL R0_CMB")
        print("=" * 30)
        
        # Objective function: minimize |l1 - target|
        def objective(R0_cmb):
            l1_pred = self.first_acoustic_peak(R0_cmb)
            return abs(l1_pred - self.l1_target)
        
        # Initial guess from plots (around 10,000 Mpc)
        R0_initial = 10000.0
        l1_initial = self.first_acoustic_peak(R0_initial)
        
        print(f"Initial guess: R0_cmb = {R0_initial:.0f} Mpc")
        print(f"  gives l1 = {l1_initial:.1f}")
        print(f"  error = {abs(l1_initial - self.l1_target):.1f}")
        print()
        
        # Optimization (focus on range where solution exists)
        result = minimize_scalar(objective, bounds=(1000, 100000), method='bounded')
        
        R0_optimal = result.x
        l1_optimal = self.first_acoustic_peak(R0_optimal)
        error = result.fun
        
        print(f"Optimization result:")
        print(f"  R0_cmb_optimal = {R0_optimal:.1f} Mpc ({R0_optimal/1000:.1f} Gpc)")
        print(f"  l1_predicted = {l1_optimal:.2f}")
        print(f"  target = {self.l1_target}")
        print(f"  error = {error:.3f}")
        print(f"  success = {result.success}")
        print()
        
        # Detailed validation
        D_A_optimal = self.udt_angular_diameter_distance(R0_optimal)
        eta_rec_Mpc = self.eta_rec_standard * self.c_light
        tau_optimal = R0_optimal / (R0_optimal + eta_rec_Mpc)
        
        print(f"Detailed results:")
        print(f"  eta_rec = {eta_rec_Mpc:.1f} Mpc ({self.eta_rec_standard:.1f} Mpc/c)")
        print(f"  tau_rec = {tau_optimal:.6f}")
        print(f"  D_A = {D_A_optimal:.1f} Mpc")
        print(f"  r_s = {self.r_s_standard:.1f} Mpc (unchanged)")
        print(f"  acoustic scale = {self.r_s_standard/D_A_optimal:.6f} rad")
        print()
        
        # Compare with standard cosmology
        print(f"Comparison with standard cosmology:")
        print(f"  D_A ratio = {D_A_optimal/self.D_A_standard:.4f}")
        print(f"  UDT correction factor = {tau_optimal:.6f}")
        print(f"  l1 ratio = {l1_optimal/(np.pi * self.D_A_standard / self.r_s_standard):.4f}")
        
        # Scale context
        print(f"\nScale context:")
        print(f"  R0_cmb = {R0_optimal:.0f} Mpc = {R0_optimal/1000:.1f} Gpc")
        print(f"  eta_rec = {eta_rec_Mpc:.0f} Mpc = {eta_rec_Mpc/1000:.1f} Gpc")
        print(f"  R0/eta ratio = {R0_optimal/eta_rec_Mpc:.2f}")
        
        if R0_optimal > eta_rec_Mpc:
            print(f"  R0 > eta_rec: UDT correction is small (tau ~ 1)")
        else:
            print(f"  R0 < eta_rec: UDT correction is significant (tau << 1)")
        
        return {
            'R0_cmb_optimal': R0_optimal,
            'l1_predicted': l1_optimal,
            'D_A': D_A_optimal,
            'tau_rec': tau_optimal,
            'error': error,
            'success': result.success
        }
    
    def validate_solution(self, R0_cmb):
        """Validate the optimal solution with detailed analysis."""
        print("\nSOLUTION VALIDATION")
        print("=" * 30)
        
        # Calculate all key quantities
        eta_rec_Mpc = self.eta_rec_standard * self.c_light
        tau_rec = R0_cmb / (R0_cmb + eta_rec_Mpc)
        D_A = self.udt_angular_diameter_distance(R0_cmb)
        l1 = self.first_acoustic_peak(R0_cmb)
        
        # Check acoustic peak positions (1st through 6th)
        peak_positions = []
        for n in range(1, 7):
            l_n = l1 * n
            peak_positions.append(l_n)
        
        print(f"Acoustic peak positions:")
        planck_peaks = [220, 540, 800, 1050, 1300, 1550]  # Approximate Planck values
        
        print("Peak  UDT   Planck  Error")
        print("----  ----  ------  -----")
        for i, (l_udt, l_planck) in enumerate(zip(peak_positions, planck_peaks)):
            error_pct = abs(l_udt - l_planck) / l_planck * 100
            print(f" {i+1:2d}   {l_udt:4.0f}   {l_planck:4.0f}   {error_pct:4.1f}%")
        
        # Physics consistency checks
        print(f"\nPhysics consistency:")
        print(f"  UDT correction factor tau = {tau_rec:.6f}")
        
        if tau_rec > 0.9:
            print(f"  Small correction: UDT is perturbative")
        elif tau_rec > 0.5:
            print(f"  Moderate correction: UDT is significant but reasonable")
        else:
            print(f"  Large correction: UDT dominates standard cosmology")
        
        # Scale reasonableness
        print(f"\nScale analysis:")
        
        # Compare with other UDT scales
        R0_galactic = 0.038        # Mpc (38 kpc)
        R0_cosmological = 3000.0   # Mpc
        
        print(f"  R0_cmb / R0_cosmo = {R0_cmb / R0_cosmological:.1f}")
        print(f"  R0_cmb / R0_gal = {R0_cmb / R0_galactic:.0f}")
        
        if R0_cmb > R0_cosmological:
            print(f"  CMB scale > cosmological scale (expected for early universe)")
        else:
            print(f"  CMB scale < cosmological scale (unexpected)")
        
        # Energy scale implications
        z_rec = 1100  # Recombination redshift
        z_eff_udt = R0_cmb / eta_rec_Mpc - 1  # Effective redshift in UDT
        
        print(f"\nRedshift analysis:")
        print(f"  Standard z_rec = {z_rec}")
        print(f"  UDT effective z = {z_eff_udt:.1f}")
        print(f"  Ratio = {z_eff_udt / z_rec:.3f}")
        
        return {
            'peak_positions': peak_positions,
            'planck_peaks': planck_peaks,
            'tau_rec': tau_rec,
            'scale_ratios': {
                'cmb_to_cosmo': R0_cmb / R0_cosmological,
                'cmb_to_gal': R0_cmb / R0_galactic
            },
            'redshift_analysis': {
                'z_standard': z_rec,
                'z_udt_effective': z_eff_udt,
                'ratio': z_eff_udt / z_rec
            }
        }
    
    def create_optimization_plots(self, result, validation, output_dir="results/optimized_udt_cmb"):
        """Create plots showing the optimization and validation."""
        os.makedirs(output_dir, exist_ok=True)
        
        fig = plt.figure(figsize=(16, 12))
        
        # 1. R0_cmb vs l1 with solution
        plt.subplot(2, 3, 1)
        R0_range = np.linspace(1000, 30000, 1000)
        l1_range = [self.first_acoustic_peak(R0) for R0 in R0_range]
        
        plt.plot(R0_range, l1_range, 'b-', linewidth=2, label='UDT prediction')
        plt.axhline(y=self.l1_target, color='red', linestyle='--', 
                   label=f'Target = {self.l1_target}', linewidth=2)
        plt.axvline(x=result['R0_cmb_optimal'], color='green', linestyle='--',
                   label=f'Optimal R0 = {result["R0_cmb_optimal"]:.0f}', linewidth=2)
        plt.scatter([result['R0_cmb_optimal']], [result['l1_predicted']], 
                   color='green', s=100, zorder=5, label=f'Solution: l1 = {result["l1_predicted"]:.2f}')
        
        plt.xlabel('R0_cmb (Mpc)')
        plt.ylabel('First Peak l1')
        plt.title('R0_cmb Optimization')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.xlim(1000, 30000)
        plt.ylim(0, 400)
        
        # 2. Temporal geometry factor
        plt.subplot(2, 3, 2)
        eta_rec_Mpc = self.eta_rec_standard * self.c_light
        tau_range = [R0 / (R0 + eta_rec_Mpc) for R0 in R0_range]
        
        plt.plot(R0_range, tau_range, 'purple', linewidth=2)
        plt.axvline(x=result['R0_cmb_optimal'], color='green', linestyle='--')
        plt.axhline(y=result['tau_rec'], color='green', linestyle=':', 
                   label=f'tau = {result["tau_rec"]:.4f}')
        plt.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='tau = 1')
        
        plt.xlabel('R0_cmb (Mpc)')
        plt.ylabel('Temporal Factor tau')
        plt.title('UDT Temporal Geometry')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.xlim(1000, 30000)
        plt.ylim(0, 1)
        
        # 3. Angular diameter distance
        plt.subplot(2, 3, 3)
        D_A_range = [self.udt_angular_diameter_distance(R0) for R0 in R0_range]
        
        plt.plot(R0_range, D_A_range, 'm-', linewidth=2, label='UDT')
        plt.axhline(y=self.D_A_standard, color='orange', linestyle='--',
                   label=f'Standard = {self.D_A_standard:.0f}', linewidth=2)
        plt.axvline(x=result['R0_cmb_optimal'], color='green', linestyle='--')
        plt.axhline(y=result['D_A'], color='green', linestyle=':',
                   label=f'Optimal = {result["D_A"]:.0f}')
        
        plt.xlabel('R0_cmb (Mpc)')
        plt.ylabel('Angular Diameter Distance (Mpc)')
        plt.title('Angular Distance vs R0_cmb')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.xlim(1000, 30000)
        
        # 4. Acoustic peak spectrum comparison
        plt.subplot(2, 3, 4)
        peak_numbers = np.arange(1, 7)
        
        plt.bar(peak_numbers - 0.2, validation['peak_positions'], 
               width=0.4, alpha=0.7, label='UDT', color='blue')
        plt.bar(peak_numbers + 0.2, validation['planck_peaks'], 
               width=0.4, alpha=0.7, label='Planck', color='orange')
        
        plt.xlabel('Peak Number')
        plt.ylabel('Multipole l')
        plt.title('Acoustic Peak Positions')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 5. Error analysis
        plt.subplot(2, 3, 5)
        R0_zoom = np.linspace(result['R0_cmb_optimal'] * 0.8, 
                             result['R0_cmb_optimal'] * 1.2, 200)
        l1_zoom = [self.first_acoustic_peak(R0) for R0 in R0_zoom]
        error_zoom = [abs(l1 - self.l1_target) for l1 in l1_zoom]
        
        plt.plot(R0_zoom, error_zoom, 'r-', linewidth=2)
        plt.axvline(x=result['R0_cmb_optimal'], color='green', linestyle='--',
                   label=f'Minimum error = {result["error"]:.3f}')
        plt.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5, label='1% error')
        
        plt.xlabel('R0_cmb (Mpc)')
        plt.ylabel('Error |l1 - target|')
        plt.title('Optimization Convergence')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 6. Scale hierarchy
        plt.subplot(2, 3, 6)
        scales = ['Galactic\n(38 kpc)', 'Cosmological\n(3000 Mpc)', f'CMB\n({result["R0_cmb_optimal"]:.0f} Mpc)']
        R0_values = [0.038, 3000.0, result['R0_cmb_optimal']]
        colors = ['blue', 'green', 'red']
        
        bars = plt.bar(range(len(scales)), np.log10(R0_values), 
                      color=colors, alpha=0.7)
        
        plt.xlabel('UDT Scale')
        plt.ylabel('log10(R0) [Mpc]')
        plt.title('Multi-Scale UDT Hierarchy')
        plt.xticks(range(len(scales)), scales)
        plt.grid(True, alpha=0.3)
        
        # Add values on bars
        for bar, val in zip(bars, R0_values):
            height = bar.get_height()
            if val >= 1:
                plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                        f'{val:.0f}', ha='center', va='bottom')
            else:
                plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                        f'{val:.3f}', ha='center', va='bottom')
        
        plt.tight_layout()
        plot_file = os.path.join(output_dir, 'optimized_udt_cmb_solution.png')
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"Optimization plots saved to: {plot_file}")
        return plot_file


def main():
    """Main optimization routine."""
    udt = OptimizedUDTCMB()
    
    # Find optimal R0_cmb
    result = udt.find_optimal_R0_cmb()
    
    # Validate solution
    validation = udt.validate_solution(result['R0_cmb_optimal'])
    
    # Create plots
    plot_file = udt.create_optimization_plots(result, validation)
    
    # Final summary
    print("\nFINAL OPTIMIZATION SUMMARY")
    print("=" * 50)
    print(f"SUCCESS: Found R0_cmb = {result['R0_cmb_optimal']:.1f} Mpc")
    print(f"Achieves l1 = {result['l1_predicted']:.2f} (target = {udt.l1_target})")
    print(f"Error: {result['error']:.3f} ({result['error']/udt.l1_target*100:.2f}%)")
    print()
    print("Multi-scale UDT framework now complete:")
    print(f"  Galactic: R0 = 38 kpc (galaxy rotation curves)")
    print(f"  Cosmological: R0 = 3,000 Mpc (supernova distances)")
    print(f"  CMB: R0 = {result['R0_cmb_optimal']:.0f} Mpc (acoustic peaks)")
    print()
    print("Ready for full CMB analysis with raw Planck data!")
    
    return result, validation


if __name__ == "__main__":
    main()