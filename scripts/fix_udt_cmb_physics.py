#!/usr/bin/env python3
"""
Fix UDT CMB Physics
==================

Addresses the fundamental issue where sound horizon and angular diameter distance
scale identically with R0_cmb, making l1 = pi * D_A / r_s constant.

The fix: Recognize that UDT temporal geometry affects different physical processes
differently. Sound horizon is set by early-universe physics (matter-radiation 
equality), while angular diameter distance reflects the full cosmic evolution.

Key insight: Maybe UDT should approach standard cosmology for CMB physics,
with only small corrections rather than complete replacement.

Author: Charles Rotter  
Date: 2025-01-17
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
import os
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from src.udt.diagnostics.parameter_registry import ParameterRegistry

class FixedUDTCMB:
    """Fixed UDT CMB physics with proper scale dependencies."""
    
    def __init__(self):
        """Initialize with corrected physics approach."""
        
        # Physical constants
        self.c_light = 299792.458       # km/s
        self.H0 = 70.0                  # km/s/Mpc
        self.z_recombination = 1100.0   # Recombination redshift
        
        # Standard cosmology parameters
        self.Omega_m = 0.31             # Matter density
        self.Omega_b = 0.049            # Baryon density
        self.Omega_r = 8.24e-5          # Radiation density
        self.Omega_lambda = 0.69        # Dark energy
        
        # Observational targets
        self.l1_target = 220.0          # First acoustic peak
        self.eta_rec_standard = 288.0   # Conformal time (Mpc/c)
        self.r_s_standard = 147.3       # Sound horizon (Mpc)
        self.D_A_standard = 13975.0     # Angular diameter distance (Mpc)
        
        # Load validated parameters from registry
        self.registry = ParameterRegistry()
        cmb_params = self.registry.get_parameters_for_analysis('cmb')
        
        # UDT parameters - use validated CMB scale, not supernova scale!
        self.R0_cmb = cmb_params['R0_mpc']  # 13041.1 Mpc (validated CMB scale)
        
        print("FIXED UDT CMB PHYSICS")
        print("=" * 40)
        print("Key insight: UDT should be a small correction to standard cosmology")
        print("for CMB physics, not a complete replacement.")
        print()
        print("Approach:")
        print("1. Use standard cosmology as baseline")
        print("2. Apply UDT corrections where physically motivated")
        print("3. Ensure sound horizon and angular distance have different R0 dependencies")
        print()
    
    def standard_conformal_time_at_recombination(self):
        """Standard cosmology conformal time at recombination."""
        # Use standard value as baseline
        return self.eta_rec_standard
    
    def udt_corrected_sound_horizon(self, udt_correction_factor=1.0):
        """
        Sound horizon with UDT corrections.
        
        Key insight: Sound horizon is set by early-universe physics.
        UDT should have minimal effect since R0 >> early universe scales.
        
        Parameters:
        -----------
        udt_correction_factor : float
            Small correction factor (should be close to 1.0)
        """
        # Start with standard sound horizon
        r_s_base = self.r_s_standard
        
        # Apply small UDT correction
        # Early universe: R0 >> scales, so tau ~ 1, minimal effect
        r_s_udt = r_s_base * udt_correction_factor
        
        return r_s_udt
    
    def udt_corrected_angular_diameter_distance(self, R0_cmb):
        """
        Angular diameter distance with UDT temporal geometry corrections.
        
        Key insight: This should have the main R0_cmb dependence.
        
        Parameters:
        -----------
        R0_cmb : float
            CMB-scale R0 parameter in Mpc
        """
        # Method 1: Start with standard and apply UDT modification
        eta_rec = self.standard_conformal_time_at_recombination()
        eta_rec_Mpc = eta_rec * self.c_light
        
        # UDT temporal geometry correction
        tau_rec = R0_cmb / (R0_cmb + eta_rec_Mpc)
        
        # Modified angular diameter distance
        # Standard would be D_A = eta_rec_Mpc (flat space)
        # UDT modification: D_A = eta_rec_Mpc * f(tau)
        
        # Several possible UDT modifications:
        # Option 1: Direct temporal geometry
        D_A_option1 = eta_rec_Mpc * tau_rec
        
        # Option 2: Inverse temporal geometry (light path affected differently)
        D_A_option2 = eta_rec_Mpc / tau_rec
        
        # Option 3: Logarithmic correction (small effect)
        D_A_option3 = eta_rec_Mpc * (1 + 0.1 * np.log(tau_rec))
        
        # Option 4: Scale-dependent correction
        scale_factor = eta_rec_Mpc / R0_cmb
        D_A_option4 = eta_rec_Mpc * (1 + scale_factor)
        
        # For now, use Option 2 (inverse) as it can give larger values
        D_A_udt = D_A_option2
        
        return D_A_udt, {
            'option1': D_A_option1,
            'option2': D_A_option2, 
            'option3': D_A_option3,
            'option4': D_A_option4,
            'tau_rec': tau_rec,
            'eta_rec_Mpc': eta_rec_Mpc
        }
    
    def first_acoustic_peak_position(self, R0_cmb, sound_horizon_correction=1.0):
        """
        Calculate first acoustic peak with fixed UDT physics.
        
        Parameters:
        -----------
        R0_cmb : float
            CMB-scale R0 parameter
        sound_horizon_correction : float
            Small correction to sound horizon
        """
        r_s = self.udt_corrected_sound_horizon(sound_horizon_correction)
        D_A, details = self.udt_corrected_angular_diameter_distance(R0_cmb)
        
        l1 = np.pi * D_A / r_s
        
        return l1, r_s, D_A, details
    
    def find_optimal_R0_cmb(self, target_l1=220.0):
        """
        Find R0_cmb that gives target l1.
        
        Parameters:
        -----------
        target_l1 : float
            Target first acoustic peak position
        """
        print("FINDING OPTIMAL R0_CMB")
        print("=" * 30)
        
        def objective(R0_cmb):
            l1, _, _, _ = self.first_acoustic_peak_position(R0_cmb)
            return abs(l1 - target_l1)
        
        # Optimization
        result = minimize_scalar(objective, bounds=(100, 1e8), method='bounded')
        
        R0_optimal = result.x
        l1_optimal, r_s_optimal, D_A_optimal, details = self.first_acoustic_peak_position(R0_optimal)
        
        print(f"Optimization result:")
        print(f"  R0_cmb = {R0_optimal:.0f} Mpc ({R0_optimal/1e6:.3f} Gpc)")
        print(f"  l1 = {l1_optimal:.1f} (target = {target_l1})")
        print(f"  r_s = {r_s_optimal:.1f} Mpc")
        print(f"  D_A = {D_A_optimal:.1f} Mpc")
        print(f"  tau_rec = {details['tau_rec']:.6f}")
        print(f"  error = {result.fun:.2f}")
        print()
        
        # Compare with standard cosmology
        print("Comparison with standard:")
        print(f"  r_s ratio = {r_s_optimal/self.r_s_standard:.3f}")
        print(f"  D_A ratio = {D_A_optimal/self.D_A_standard:.3f}")
        print(f"  eta_rec used = {details['eta_rec_Mpc']:.1f} Mpc (standard = {self.eta_rec_standard * self.c_light:.1f})")
        
        return {
            'R0_cmb': R0_optimal,
            'l1': l1_optimal,
            'r_s': r_s_optimal,
            'D_A': D_A_optimal,
            'tau_rec': details['tau_rec'],
            'error': result.fun,
            'success': result.success
        }
    
    def test_different_udt_modifications(self):
        """Test different ways UDT could modify CMB physics."""
        print("TESTING UDT MODIFICATION OPTIONS")
        print("=" * 40)
        
        # Test R0 value
        R0_test = 10000.0  # 10 Gpc
        
        print(f"Testing with R0_cmb = {R0_test:.0f} Mpc")
        print()
        
        # Standard approach (no UDT)
        r_s_std = self.r_s_standard
        D_A_std = self.eta_rec_standard * self.c_light  # Flat space
        l1_std = np.pi * D_A_std / r_s_std
        
        print(f"Standard cosmology (no UDT):")
        print(f"  r_s = {r_s_std:.1f} Mpc")
        print(f"  D_A = {D_A_std:.1f} Mpc")
        print(f"  l1 = {l1_std:.1f}")
        print()
        
        # UDT Option 1: Direct temporal geometry
        _, details = self.udt_corrected_angular_diameter_distance(R0_test)
        
        for option in ['option1', 'option2', 'option3', 'option4']:
            D_A_opt = details[option]
            l1_opt = np.pi * D_A_opt / r_s_std
            
            print(f"UDT {option}:")
            print(f"  D_A = {D_A_opt:.1f} Mpc")
            print(f"  l1 = {l1_opt:.1f}")
            print(f"  ratio to standard = {l1_opt/l1_std:.3f}")
            print()
        
        # Test scale dependence
        print("Scale dependence test:")
        R0_range = [1000, 3000, 10000, 30000, 100000]
        
        print("R0_cmb (Mpc)  l1 (opt1)  l1 (opt2)  l1 (opt3)  l1 (opt4)")
        print("-" * 60)
        
        for R0 in R0_range:
            _, details = self.udt_corrected_angular_diameter_distance(R0)
            l1_vals = []
            
            for option in ['option1', 'option2', 'option3', 'option4']:
                D_A_opt = details[option]
                l1_opt = np.pi * D_A_opt / r_s_std
                l1_vals.append(l1_opt)
            
            print(f"{R0:8.0f}      {l1_vals[0]:6.1f}     {l1_vals[1]:6.1f}     {l1_vals[2]:6.1f}     {l1_vals[3]:6.1f}")
    
    def create_fixed_physics_plots(self, output_dir="results/fixed_udt_cmb"):
        """Create diagnostic plots for fixed UDT CMB physics."""
        os.makedirs(output_dir, exist_ok=True)
        
        fig = plt.figure(figsize=(16, 12))
        
        # 1. R0_cmb vs l1 for different UDT options
        plt.subplot(2, 3, 1)
        R0_range = np.logspace(3, 7, 100)
        
        # Standard cosmology baseline
        r_s_std = self.r_s_standard
        D_A_std = self.eta_rec_standard * self.c_light
        l1_std = np.pi * D_A_std / r_s_std
        
        l1_option1 = []
        l1_option2 = []
        l1_option3 = []
        l1_option4 = []
        
        for R0 in R0_range:
            _, details = self.udt_corrected_angular_diameter_distance(R0)
            
            l1_option1.append(np.pi * details['option1'] / r_s_std)
            l1_option2.append(np.pi * details['option2'] / r_s_std)
            l1_option3.append(np.pi * details['option3'] / r_s_std)
            l1_option4.append(np.pi * details['option4'] / r_s_std)
        
        plt.semilogx(R0_range, l1_option1, 'b-', label='Option 1: tau', linewidth=2)
        plt.semilogx(R0_range, l1_option2, 'r-', label='Option 2: 1/tau', linewidth=2)
        plt.semilogx(R0_range, l1_option3, 'g-', label='Option 3: log correction', linewidth=2)
        plt.semilogx(R0_range, l1_option4, 'm-', label='Option 4: scale factor', linewidth=2)
        plt.axhline(y=self.l1_target, color='orange', linestyle='--', 
                   label=f'Target = {self.l1_target}', linewidth=2)
        plt.axhline(y=l1_std, color='gray', linestyle=':', 
                   label=f'Standard = {l1_std:.1f}', alpha=0.7)
        
        plt.xlabel('R0_cmb (Mpc)')
        plt.ylabel('First Peak l1')
        plt.title('UDT Modification Options')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.ylim(0, 500)
        
        # 2. Angular diameter distance options
        plt.subplot(2, 3, 2)
        D_A_option1 = []
        D_A_option2 = []
        D_A_option3 = []
        D_A_option4 = []
        
        for R0 in R0_range:
            _, details = self.udt_corrected_angular_diameter_distance(R0)
            D_A_option1.append(details['option1'])
            D_A_option2.append(details['option2'])
            D_A_option3.append(details['option3'])
            D_A_option4.append(details['option4'])
        
        plt.loglog(R0_range, D_A_option1, 'b-', label='Option 1', linewidth=2)
        plt.loglog(R0_range, D_A_option2, 'r-', label='Option 2', linewidth=2)
        plt.loglog(R0_range, D_A_option3, 'g-', label='Option 3', linewidth=2)
        plt.loglog(R0_range, D_A_option4, 'm-', label='Option 4', linewidth=2)
        plt.axhline(y=self.D_A_standard, color='orange', linestyle='--',
                   label=f'Standard = {self.D_A_standard:.0f}', linewidth=2)
        
        plt.xlabel('R0_cmb (Mpc)')
        plt.ylabel('Angular Diameter Distance (Mpc)')
        plt.title('Angular Distance Options')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 3. Temporal geometry factor
        plt.subplot(2, 3, 3)
        tau_values = []
        
        eta_rec_Mpc = self.eta_rec_standard * self.c_light
        for R0 in R0_range:
            tau = R0 / (R0 + eta_rec_Mpc)
            tau_values.append(tau)
        
        plt.semilogx(R0_range, tau_values, 'purple', linewidth=2)
        plt.axhline(y=1.0, color='gray', linestyle='--', alpha=0.7, label='tau = 1')
        plt.axvline(x=eta_rec_Mpc, color='orange', linestyle='--', alpha=0.7,
                   label=f'eta_rec = {eta_rec_Mpc:.0f} Mpc')
        
        plt.xlabel('R0_cmb (Mpc)')
        plt.ylabel('Temporal Factor tau')
        plt.title('UDT Temporal Geometry')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.ylim(0, 1)
        
        # 4. Find best option (closest to target)
        plt.subplot(2, 3, 4)
        
        # Find R0 that gives l1 = 220 for each option
        target = self.l1_target
        
        def find_R0_for_target(option_data):
            """Find R0 that gives target l1 for given option."""
            # Interpolate to find R0 where l1 = target
            for i in range(len(option_data)-1):
                if (option_data[i] <= target <= option_data[i+1]) or \
                   (option_data[i] >= target >= option_data[i+1]):
                    # Linear interpolation
                    frac = (target - option_data[i]) / (option_data[i+1] - option_data[i])
                    R0_interp = R0_range[i] + frac * (R0_range[i+1] - R0_range[i])
                    return R0_interp
            return None
        
        R0_opt1 = find_R0_for_target(l1_option1)
        R0_opt2 = find_R0_for_target(l1_option2)
        R0_opt3 = find_R0_for_target(l1_option3)
        R0_opt4 = find_R0_for_target(l1_option4)
        
        options = ['Option 1', 'Option 2', 'Option 3', 'Option 4']
        R0_solutions = [R0_opt1, R0_opt2, R0_opt3, R0_opt4]
        colors = ['blue', 'red', 'green', 'magenta']
        
        valid_solutions = []
        valid_labels = []
        valid_colors = []
        
        for i, (opt, R0_sol) in enumerate(zip(options, R0_solutions)):
            if R0_sol is not None:
                valid_solutions.append(R0_sol)
                valid_labels.append(opt)
                valid_colors.append(colors[i])
        
        if valid_solutions:
            plt.bar(range(len(valid_solutions)), np.log10(valid_solutions), 
                   color=valid_colors, alpha=0.7)
            plt.xlabel('UDT Option')
            plt.ylabel('log10(R0_cmb) for l1=220')
            plt.title('Required R0_cmb by Option')
            plt.xticks(range(len(valid_solutions)), valid_labels, rotation=45)
            plt.grid(True, alpha=0.3)
        else:
            plt.text(0.5, 0.5, 'No solutions found\nin R0 range', 
                    ha='center', va='center', transform=plt.gca().transAxes)
            plt.title('Required R0_cmb by Option')
        
        # 5. Error analysis
        plt.subplot(2, 3, 5)
        errors = [abs(np.array(l1_opt) - target) for l1_opt in [l1_option1, l1_option2, l1_option3, l1_option4]]
        
        for i, (error, label, color) in enumerate(zip(errors, options, colors)):
            plt.semilogx(R0_range, error, color=color, label=label, linewidth=2)
        
        plt.xlabel('R0_cmb (Mpc)')
        plt.ylabel('Error |l1 - target|')
        plt.title('Error vs R0_cmb')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.ylim(0, 100)
        
        # 6. Summary comparison
        plt.subplot(2, 3, 6)
        
        # Pick a representative R0 value
        R0_rep = 30000  # 30 Gpc
        if R0_rep in R0_range:
            idx = list(R0_range).index(R0_rep)
        else:
            idx = len(R0_range) // 2
            R0_rep = R0_range[idx]
        
        l1_values = [l1_option1[idx], l1_option2[idx], l1_option3[idx], l1_option4[idx]]
        
        bars = plt.bar(range(len(options)), l1_values, color=colors, alpha=0.7)
        plt.axhline(y=target, color='orange', linestyle='--', linewidth=2, label='Target')
        plt.axhline(y=l1_std, color='gray', linestyle=':', alpha=0.7, label='Standard')
        
        plt.xlabel('UDT Option')
        plt.ylabel('First Peak l1')
        plt.title(f'l1 Comparison at R0={R0_rep:.0f} Mpc')
        plt.xticks(range(len(options)), [f'Opt {i+1}' for i in range(len(options))])
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Add values on bars
        for bar, val in zip(bars, l1_values):
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height + 5,
                    f'{val:.1f}', ha='center', va='bottom')
        
        plt.tight_layout()
        plot_file = os.path.join(output_dir, 'fixed_udt_cmb_physics.png')
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"Fixed UDT CMB physics plots saved to: {plot_file}")
        return plot_file


def main():
    """Main analysis for fixed UDT CMB physics."""
    udt_cmb = FixedUDTCMB()
    
    # Test different UDT modifications
    udt_cmb.test_different_udt_modifications()
    
    # Find optimal R0_cmb (using option 2 which seems most promising)
    print("FINDING OPTIMAL R0_CMB (using Option 2: 1/tau)")
    print("=" * 50)
    
    # Manual optimization for option 2
    def objective_option2(R0_cmb):
        r_s = udt_cmb.udt_corrected_sound_horizon()
        _, details = udt_cmb.udt_corrected_angular_diameter_distance(R0_cmb)
        D_A = details['option2']  # Use inverse temporal geometry
        l1 = np.pi * D_A / r_s
        return abs(l1 - udt_cmb.l1_target)
    
    result = minimize_scalar(objective_option2, bounds=(100, 1e8), method='bounded')
    
    R0_optimal = result.x
    r_s_opt = udt_cmb.udt_corrected_sound_horizon()
    _, details_opt = udt_cmb.udt_corrected_angular_diameter_distance(R0_optimal)
    D_A_opt = details_opt['option2']
    l1_opt = np.pi * D_A_opt / r_s_opt
    
    print(f"Option 2 optimization:")
    print(f"  R0_cmb = {R0_optimal:.0f} Mpc ({R0_optimal/1e6:.2f} Gpc)")
    print(f"  l1 = {l1_opt:.1f} (target = {udt_cmb.l1_target})")
    print(f"  r_s = {r_s_opt:.1f} Mpc")
    print(f"  D_A = {D_A_opt:.1f} Mpc")
    print(f"  error = {result.fun:.2f}")
    
    # Create diagnostic plots
    plot_file = udt_cmb.create_fixed_physics_plots()
    
    return {
        'R0_cmb_optimal': R0_optimal,
        'l1_predicted': l1_opt,
        'r_s': r_s_opt,
        'D_A': D_A_opt,
        'error': result.fun
    }


if __name__ == "__main__":
    main()