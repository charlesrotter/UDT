#!/usr/bin/env python3
"""
Test Derived UDT Kernel Against SPARC Data
==========================================

Test our mathematically derived kernel:
K(r,r') = G × [1 + α × f(τ(r), τ(r'))] / |r-r'|²

where f(τ, 1) = C × (1-τ)²/[τ(3-2τ)]

This should produce v = V_scale × √(r/(r + R₀/3)) × (1 + r/R₀)

Author: Charles Rotter
Date: 2025-07-18
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from pathlib import Path

class TestDerivedKernel:
    def __init__(self):
        print("TESTING DERIVED UDT KERNEL AGAINST SPARC DATA")
        print("=" * 50)
        print("Testing mathematically derived kernel from distance equivalence")
        print("=" * 50)
        
        # Physical constants
        self.G = 4.302e-6  # km^2/s^2 per solar mass per kpc
        
        # UDT parameters to optimize
        self.R0_gal = 62.4  # kpc (from successful phenomenology)
        self.alpha = 1.0    # Coupling constant (to be optimized)
        
        print(f"Initial parameters:")
        print(f"R0_gal = {self.R0_gal} kpc")
        print(f"alpha = {self.alpha}")
        print()
    
    def tau(self, r):
        """Distance equivalence principle."""
        return self.R0_gal / (self.R0_gal + r)
    
    def f_modification(self, tau_val):
        """Our derived modification function."""
        return (1 - tau_val)**2 / (tau_val * (3 - 2*tau_val))
    
    def effective_gravity(self, r):
        """Effective gravitational constant from UDT kernel."""
        tau_val = self.tau(r)
        f_val = self.f_modification(tau_val)
        return self.G * (1 + self.alpha * f_val)
    
    def derived_velocity_profile(self, r, M_central):
        """
        Velocity profile from derived UDT kernel.
        
        For central mass M_central, the velocity is:
        v² = G_eff(r) × M_central / r
        """
        G_eff = self.effective_gravity(r)
        return np.sqrt(G_eff * M_central / r)
    
    def theoretical_velocity_profile(self, r, V_scale):
        """
        Theoretical velocity profile we derived:
        v = V_scale × √(r/(r + R₀/3)) × (1 + r/R₀)
        """
        base_profile = np.sqrt(r / (r + self.R0_gal/3))
        enhancement = (1 + r/self.R0_gal)
        return V_scale * base_profile * enhancement
    
    def load_sparc_galaxy(self, galaxy_name):
        """Load a specific SPARC galaxy."""
        sparc_dir = Path("C:/UDT/data/sparc_database")
        file_path = sparc_dir / f"{galaxy_name}_rotmod.dat"
        
        if not file_path.exists():
            return None
        
        data = []
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        rad = float(parts[0])      # kpc
                        vobs = float(parts[1])     # km/s
                        errv = float(parts[2])     # km/s
                        
                        if rad > 0 and vobs > 0 and errv > 0:
                            data.append([rad, vobs, errv])
                    except ValueError:
                        continue
        
        if len(data) < 5:
            return None
        
        data = np.array(data)
        return {
            'name': galaxy_name,
            'radius': data[:, 0],
            'velocity': data[:, 1],
            'velocity_error': data[:, 2]
        }
    
    def fit_with_derived_kernel(self, galaxy):
        """Fit galaxy using our derived kernel."""
        r = galaxy['radius']
        v_obs = galaxy['velocity']
        v_err = galaxy['velocity_error']
        
        weights = 1 / v_err**2
        
        def objective(M_central):
            if M_central <= 0:
                return 1e10
            
            try:
                v_model = self.derived_velocity_profile(r, M_central)
                chi2 = np.sum(weights * (v_obs - v_model)**2)
                return chi2
            except:
                return 1e10
        
        # Optimize central mass
        result = minimize_scalar(objective, bounds=(1e8, 1e13), method='bounded')
        
        if result.success:
            M_best = result.x
            chi2_best = result.fun
            dof = len(r) - 1
            chi2_per_dof = chi2_best / dof
            
            v_model = self.derived_velocity_profile(r, M_best)
            rms = np.sqrt(np.mean((v_obs - v_model)**2))
            
            return {
                'success': True,
                'M_central': M_best,
                'chi2_per_dof': chi2_per_dof,
                'rms': rms,
                'v_model': v_model
            }
        else:
            return {'success': False}
    
    def fit_with_theoretical_formula(self, galaxy):
        """Fit galaxy using theoretical formula."""
        r = galaxy['radius']
        v_obs = galaxy['velocity']
        v_err = galaxy['velocity_error']
        
        weights = 1 / v_err**2
        
        def objective(V_scale):
            if V_scale <= 0:
                return 1e10
            
            try:
                v_model = self.theoretical_velocity_profile(r, V_scale)
                chi2 = np.sum(weights * (v_obs - v_model)**2)
                return chi2
            except:
                return 1e10
        
        # Optimize velocity scale
        result = minimize_scalar(objective, bounds=(10, 500), method='bounded')
        
        if result.success:
            V_best = result.x
            chi2_best = result.fun
            dof = len(r) - 1
            chi2_per_dof = chi2_best / dof
            
            v_model = self.theoretical_velocity_profile(r, V_best)
            rms = np.sqrt(np.mean((v_obs - v_model)**2))
            
            return {
                'success': True,
                'V_scale': V_best,
                'chi2_per_dof': chi2_per_dof,
                'rms': rms,
                'v_model': v_model
            }
        else:
            return {'success': False}
    
    def test_single_galaxy(self, galaxy_name):
        """Test both approaches on a single galaxy."""
        print(f"TESTING GALAXY: {galaxy_name}")
        print("-" * 30)
        
        galaxy = self.load_sparc_galaxy(galaxy_name)
        if galaxy is None:
            print(f"Could not load galaxy {galaxy_name}")
            return None
        
        print(f"Loaded {len(galaxy['radius'])} data points")
        print(f"Radius range: {galaxy['radius'].min():.1f} - {galaxy['radius'].max():.1f} kpc")
        print(f"Velocity range: {galaxy['velocity'].min():.1f} - {galaxy['velocity'].max():.1f} km/s")
        print()
        
        # Test derived kernel
        print("1. DERIVED KERNEL FIT:")
        kernel_fit = self.fit_with_derived_kernel(galaxy)
        if kernel_fit['success']:
            print(f"   M_central = {kernel_fit['M_central']:.2e} M_sun")
            print(f"   chi2/DOF = {kernel_fit['chi2_per_dof']:.2f}")
            print(f"   RMS = {kernel_fit['rms']:.1f} km/s")
        else:
            print("   FAILED")
        
        # Test theoretical formula
        print("2. THEORETICAL FORMULA FIT:")
        theory_fit = self.fit_with_theoretical_formula(galaxy)
        if theory_fit['success']:
            print(f"   V_scale = {theory_fit['V_scale']:.1f} km/s")
            print(f"   chi2/DOF = {theory_fit['chi2_per_dof']:.2f}")
            print(f"   RMS = {theory_fit['rms']:.1f} km/s")
        else:
            print("   FAILED")
        
        print()
        
        # Create comparison plot
        if kernel_fit['success'] and theory_fit['success']:
            self.plot_comparison(galaxy, kernel_fit, theory_fit)
        
        return {
            'galaxy': galaxy,
            'kernel_fit': kernel_fit,
            'theory_fit': theory_fit
        }
    
    def plot_comparison(self, galaxy, kernel_fit, theory_fit):
        """Plot comparison of fits."""
        fig, axes = plt.subplots(2, 1, figsize=(10, 10))
        
        r = galaxy['radius']
        v_obs = galaxy['velocity']
        v_err = galaxy['velocity_error']
        
        # Top panel: Velocity curves
        ax1 = axes[0]
        ax1.errorbar(r, v_obs, yerr=v_err, fmt='o', color='black', 
                    alpha=0.7, label='SPARC Data')
        ax1.plot(r, kernel_fit['v_model'], 'r-', linewidth=2, 
                label=f'Derived Kernel (χ²/ν={kernel_fit["chi2_per_dof"]:.2f})')
        ax1.plot(r, theory_fit['v_model'], 'b--', linewidth=2, 
                label=f'Theoretical Formula (χ²/ν={theory_fit["chi2_per_dof"]:.2f})')
        
        ax1.set_xlabel('Radius (kpc)')
        ax1.set_ylabel('Velocity (km/s)')
        ax1.set_title(f'{galaxy["name"]} - Rotation Curve Fits')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Bottom panel: Residuals
        ax2 = axes[1]
        res_kernel = v_obs - kernel_fit['v_model']
        res_theory = v_obs - theory_fit['v_model']
        
        ax2.errorbar(r, res_kernel, yerr=v_err, fmt='ro', alpha=0.7, 
                    label='Derived Kernel Residuals')
        ax2.errorbar(r, res_theory, yerr=v_err, fmt='bo', alpha=0.7, 
                    label='Theoretical Formula Residuals')
        ax2.axhline(y=0, color='k', linestyle='-', alpha=0.3)
        
        ax2.set_xlabel('Radius (kpc)')
        ax2.set_ylabel('Residuals (km/s)')
        ax2.set_title('Fit Residuals')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f'C:/UDT/results/kernel_test_{galaxy["name"]}.png', dpi=150)
        plt.close()
        
        print(f"Saved: C:/UDT/results/kernel_test_{galaxy['name']}.png")
    
    def test_multiple_galaxies(self):
        """Test on multiple galaxies."""
        print("TESTING MULTIPLE GALAXIES")
        print("-" * 30)
        
        # Test on a variety of galaxies
        test_galaxies = [
            'NGC3198',  # Classic flat rotation curve
            'DDO154',   # Small galaxy
            'NGC2403',  # Medium galaxy
            'NGC6946',  # Large galaxy
            'UGC07089'  # Another test case
        ]
        
        results = []
        successful_kernel = 0
        successful_theory = 0
        
        for galaxy_name in test_galaxies:
            print(f"\n{'='*50}")
            result = self.test_single_galaxy(galaxy_name)
            
            if result is not None:
                results.append(result)
                if result['kernel_fit']['success']:
                    successful_kernel += 1
                if result['theory_fit']['success']:
                    successful_theory += 1
        
        print(f"\n{'='*50}")
        print("SUMMARY RESULTS")
        print("=" * 50)
        
        print(f"Galaxies tested: {len(results)}")
        print(f"Derived kernel successful: {successful_kernel}/{len(results)}")
        print(f"Theoretical formula successful: {successful_theory}/{len(results)}")
        
        if results:
            # Calculate statistics
            kernel_chi2 = [r['kernel_fit']['chi2_per_dof'] for r in results 
                          if r['kernel_fit']['success']]
            theory_chi2 = [r['theory_fit']['chi2_per_dof'] for r in results 
                          if r['theory_fit']['success']]
            
            if kernel_chi2:
                print(f"\nDERIVED KERNEL STATISTICS:")
                print(f"Mean chi2/DOF: {np.mean(kernel_chi2):.2f} ± {np.std(kernel_chi2):.2f}")
                print(f"Median chi2/DOF: {np.median(kernel_chi2):.2f}")
            
            if theory_chi2:
                print(f"\nTHEORETICAL FORMULA STATISTICS:")
                print(f"Mean chi2/DOF: {np.mean(theory_chi2):.2f} ± {np.std(theory_chi2):.2f}")
                print(f"Median chi2/DOF: {np.median(theory_chi2):.2f}")
        
        print(f"\n{'='*50}")
        print("VALIDATION CONCLUSION")
        print("=" * 50)
        
        if len(kernel_chi2) > 0 and len(theory_chi2) > 0:
            kernel_mean = np.mean(kernel_chi2)
            theory_mean = np.mean(theory_chi2)
            
            print(f"Derived kernel mean chi2/DOF: {kernel_mean:.2f}")
            print(f"Theoretical formula mean chi2/DOF: {theory_mean:.2f}")
            print(f"Difference: {abs(kernel_mean - theory_mean):.2f}")
            
            if abs(kernel_mean - theory_mean) < 0.5:
                print("\nSUCCESS: Derived kernel matches theoretical formula!")
                print("The mathematical derivation is validated by real data.")
            else:
                print("\nDISCREPANCY: Derived kernel differs from theoretical formula.")
                print("May need to adjust coupling constant alpha or kernel form.")
        
        return results

def main():
    """Run the kernel validation test."""
    tester = TestDerivedKernel()
    results = tester.test_multiple_galaxies()
    return results

if __name__ == "__main__":
    main()