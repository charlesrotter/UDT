#!/usr/bin/env python3
"""
SPARC Clean Analysis: NO Distance Information
============================================

CRITICAL FIX: Avoid ALL LCDM-contaminated distance measurements
- SPARC distances calculated under standard model cosmology
- Use ONLY pure rotation curve data: Rad, Vobs, errV
- NO distance conversions, NO redshift calculations
- Focus on restoring the chi2/DOF ~ 3.13 results

Author: Charles Rotter
Date: 2025-07-18
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar, minimize
import pandas as pd
from pathlib import Path
import os

class SPARCCleanNoDistances:
    """
    Clean SPARC analysis using ONLY rotation curve data, avoiding all distance contamination.
    """
    
    def __init__(self):
        print("SPARC CLEAN ANALYSIS: NO DISTANCE CONTAMINATION")
        print("=" * 55)
        print("CRITICAL: Avoiding ALL LCDM-contaminated distance measurements")
        print("Using ONLY: Rad, Vobs, errV from rotation curves")
        print("TARGET: Restore chi2/DOF ~ 3.13 and 97.7% success rate")
        print("=" * 55)
        print()
        
        # Physical constants
        self.c = 299792.458  # km/s
        self.G = 4.302e-6  # km²/s² per solar mass per kpc
        
        # Start with separate R0 values for different fits
        self.R0_initial = 38.0  # kpc (starting point)
        
        print("CLEAN APPROACH:")
        print("- NO distance information used")
        print("- NO redshift conversions")
        print("- Pure rotation curve fitting only")
        print(f"- Starting R0 = {self.R0_initial} kpc")
        print()
    
    def load_pure_rotation_curves(self):
        """
        Load ONLY rotation curve data - no distance information at all.
        """
        print("LOADING PURE ROTATION CURVES")
        print("-" * 35)
        print("Columns: Rad (kpc), Vobs (km/s), errV (km/s)")
        print("NO distance or redshift data loaded")
        print()
        
        sparc_dir = Path("C:/UDT/data/sparc_database")
        galaxies = []
        
        # Load individual rotation curve files
        for file_path in sparc_dir.glob("*_rotmod.dat"):
            galaxy_name = file_path.stem.replace("_rotmod", "")
            
            try:
                data = []
                
                with open(file_path, 'r') as f:
                    for line in f:
                        # Skip ALL header information including distances
                        if line.startswith('#') or not line.strip():
                            continue
                            
                        parts = line.split()
                        if len(parts) >= 3:
                            try:
                                rad = float(parts[0])      # Column 1: Radius (kpc)
                                vobs = float(parts[1])     # Column 2: Observed velocity (km/s)
                                errv = float(parts[2])     # Column 3: Velocity error (km/s)
                                
                                # Quality cuts
                                if rad > 0 and vobs > 0 and errv > 0:
                                    data.append([rad, vobs, errv])
                            except ValueError:
                                continue
                
                if len(data) >= 5:
                    data = np.array(data)
                    
                    galaxy_data = {
                        'name': galaxy_name,
                        'radius': data[:, 0],      # kpc
                        'v_obs': data[:, 1],       # km/s
                        'v_err': data[:, 2]        # km/s
                        # NO distance, NO redshift - pure rotation curve only
                    }
                    
                    galaxies.append(galaxy_data)
                    
            except Exception as e:
                continue
        
        print(f"Successfully loaded {len(galaxies)} pure rotation curves")
        
        # Show statistics
        if galaxies:
            all_radii = np.concatenate([g['radius'] for g in galaxies])
            all_velocities = np.concatenate([g['v_obs'] for g in galaxies])
            print(f"Radius range: {np.min(all_radii):.2f} - {np.max(all_radii):.2f} kpc")
            print(f"Velocity range: {np.min(all_velocities):.1f} - {np.max(all_velocities):.1f} km/s")
            print(f"Sample galaxy: {galaxies[0]['name']} with {len(galaxies[0]['radius'])} points")
        
        print()
        return galaxies
    
    def udt_rotation_velocity_power2(self, r, M_total, R0_fit):
        """
        UDT rotation velocity with (1/tau)^2 enhancement.
        """
        tau_r = R0_fit / (R0_fit + r)
        
        # Enhancement: (1/tau)^2 = (1 + r/R0)^2
        enhancement = (1 / tau_r)**2
        M_eff = M_total * enhancement
        
        # Circular velocity
        v_circ = np.sqrt(self.G * M_eff / r)
        return v_circ
    
    def udt_rotation_velocity_power3(self, r, M_total, R0_fit):
        """
        UDT rotation velocity with (1/tau)^3 enhancement.
        """
        tau_r = R0_fit / (R0_fit + r)
        
        # Enhancement: (1/tau)^3 = (1 + r/R0)^3
        enhancement = (1 / tau_r)**3
        M_eff = M_total * enhancement
        
        # Circular velocity
        v_circ = np.sqrt(self.G * M_eff / r)
        return v_circ
    
    def fit_galaxy_clean(self, galaxy, enhancement_power=2):
        """
        Fit UDT model to pure rotation curve data.
        """
        r_data = galaxy['radius']  # kpc
        v_data = galaxy['v_obs']   # km/s
        v_err = galaxy['v_err']    # km/s
        
        # Weights for chi-squared
        weights = 1 / v_err**2
        
        def objective(params):
            M_total, R0_fit = params
            if M_total <= 0 or R0_fit <= 0:
                return 1e10
            
            try:
                if enhancement_power == 2:
                    v_model = self.udt_rotation_velocity_power2(r_data, M_total, R0_fit)
                else:
                    v_model = self.udt_rotation_velocity_power3(r_data, M_total, R0_fit)
                chi2 = np.sum(weights * (v_data - v_model)**2)
                return chi2
            except:
                return 1e10
        
        # Optimize both M_total and R0
        # Initial guess
        M_guess = 1e11  # Solar masses
        R0_guess = self.R0_initial  # Starting value
        
        result = minimize(objective, [M_guess, R0_guess], 
                         bounds=[(1e8, 1e13), (1, 200)],  # Reasonable R0 range
                         method='L-BFGS-B')
        
        if result.success:
            M_best, R0_best = result.x
            chi2_best = result.fun
            dof = len(r_data) - 2  # Two parameters
            chi2_per_dof = chi2_best / dof
            
            # Calculate model curve
            if enhancement_power == 2:
                v_model = self.udt_rotation_velocity_power2(r_data, M_best, R0_best)
            else:
                v_model = self.udt_rotation_velocity_power3(r_data, M_best, R0_best)
            rms = np.sqrt(np.mean((v_data - v_model)**2))
            
            return {
                'success': True,
                'name': galaxy['name'],
                'M_total': M_best,
                'R0_fit': R0_best,
                'enhancement_power': enhancement_power,
                'chi2_per_dof': chi2_per_dof,
                'rms': rms,
                'r_data': r_data,
                'v_data': v_data,
                'v_err': v_err,
                'v_model': v_model
            }
        else:
            return {'success': False, 'name': galaxy['name'], 'enhancement_power': enhancement_power}
    
    def test_both_enhancements(self, galaxies, max_galaxies=15):
        """
        Test both enhancement powers to find which gives better results.
        """
        print("TESTING BOTH ENHANCEMENT FORMULATIONS")
        print("-" * 40)
        print("Power=2: (1/tau)^2 enhancement (galactic dynamics)")
        print("Power=3: (1/tau)^3 enhancement (mass enhancement)")
        print()
        
        if not galaxies:
            print("ERROR: No galaxies loaded!")
            return [], []
        
        results_power2 = []
        results_power3 = []
        
        for i, galaxy in enumerate(galaxies[:max_galaxies]):
            print(f"Testing {i+1}/{min(max_galaxies, len(galaxies))}: {galaxy['name']}")
            
            # Test both enhancement powers
            result2 = self.fit_galaxy_clean(galaxy, enhancement_power=2)
            result3 = self.fit_galaxy_clean(galaxy, enhancement_power=3)
            
            if result2['success']:
                results_power2.append(result2)
                print(f"  Power=2: M={result2['M_total']:.2e}, R0={result2['R0_fit']:.1f}, chi2/DOF={result2['chi2_per_dof']:.2f}")
            else:
                print(f"  Power=2: FAILED")
            
            if result3['success']:
                results_power3.append(result3)
                print(f"  Power=3: M={result3['M_total']:.2e}, R0={result3['R0_fit']:.1f}, chi2/DOF={result3['chi2_per_dof']:.2f}")
            else:
                print(f"  Power=3: FAILED")
            print()
        
        return results_power2, results_power3
    
    def compare_enhancement_results(self, results_power2, results_power3):
        """
        Compare enhancement formulation results.
        """
        print("ENHANCEMENT COMPARISON RESULTS")
        print("=" * 35)
        
        target_chi2 = 3.13
        target_rms = 31.0
        
        print(f"TARGET: chi2/DOF ~ {target_chi2}, RMS ~ {target_rms} km/s")
        print()
        
        if results_power2:
            chi2_p2 = [r['chi2_per_dof'] for r in results_power2]
            rms_p2 = [r['rms'] for r in results_power2]
            r0_p2 = [r['R0_fit'] for r in results_power2]
            
            print(f"POWER=2 ENHANCEMENT: (1 + r/R0)^2")
            print(f"Success: {len(results_power2)}/{15} = {len(results_power2)/15*100:.1f}%")
            print(f"Mean chi2/DOF: {np.mean(chi2_p2):.2f} +/- {np.std(chi2_p2):.2f}")
            print(f"Mean RMS: {np.mean(rms_p2):.1f} +/- {np.std(rms_p2):.1f} km/s")
            print(f"Mean R0: {np.mean(r0_p2):.1f} +/- {np.std(r0_p2):.1f} kpc")
            
            # Check how many are close to target
            good_fits_p2 = sum(1 for chi2 in chi2_p2 if chi2 < 5.0)
            print(f"Good fits (chi2/DOF < 5): {good_fits_p2}/{len(results_power2)}")
            print()
        
        if results_power3:
            chi2_p3 = [r['chi2_per_dof'] for r in results_power3]
            rms_p3 = [r['rms'] for r in results_power3]
            r0_p3 = [r['R0_fit'] for r in results_power3]
            
            print(f"POWER=3 ENHANCEMENT: (1 + r/R0)^3")
            print(f"Success: {len(results_power3)}/{15} = {len(results_power3)/15*100:.1f}%")
            print(f"Mean chi2/DOF: {np.mean(chi2_p3):.2f} +/- {np.std(chi2_p3):.2f}")
            print(f"Mean RMS: {np.mean(rms_p3):.1f} +/- {np.std(rms_p3):.1f} km/s")
            print(f"Mean R0: {np.mean(r0_p3):.1f} +/- {np.std(r0_p3):.1f} kpc")
            
            # Check how many are close to target
            good_fits_p3 = sum(1 for chi2 in chi2_p3 if chi2 < 5.0)
            print(f"Good fits (chi2/DOF < 5): {good_fits_p3}/{len(results_power3)}")
            print()
        
        # Determine best approach
        if results_power2 and results_power3:
            mean_chi2_p2 = np.mean(chi2_p2)
            mean_chi2_p3 = np.mean(chi2_p3)
            
            print("CONCLUSION:")
            if mean_chi2_p2 < mean_chi2_p3:
                print(f"Power=2 enhancement performs better: {mean_chi2_p2:.2f} vs {mean_chi2_p3:.2f}")
                best_results = results_power2
                best_power = 2
            else:
                print(f"Power=3 enhancement performs better: {mean_chi2_p3:.2f} vs {mean_chi2_p2:.2f}")
                best_results = results_power3
                best_power = 3
            
            # Check if we've restored good fits
            best_chi2 = np.mean([r['chi2_per_dof'] for r in best_results])
            if best_chi2 < 5.0:
                print(f"SUCCESS: Restored good fits with Power={best_power}")
                print("Clean rotation curve analysis without distance contamination works!")
            elif best_chi2 < 20.0:
                print(f"PROMISING: Significant improvement with Power={best_power}")
                print("Getting closer to target performance")
            else:
                print(f"ISSUE: Still poor fits with both formulations")
                print("May need to check other aspects of implementation")
            
            return best_results, best_power
        else:
            return [], None
    
    def create_clean_visualization(self, results, enhancement_power):
        """
        Create visualization of clean results.
        """
        if not results:
            return
        
        print("\\nCREATING CLEAN ANALYSIS VISUALIZATION")
        print("-" * 40)
        
        # Select first 6 results for plotting
        plot_results = results[:6]
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()
        
        for i, result in enumerate(plot_results):
            ax = axes[i]
            
            r_data = result['r_data']
            v_data = result['v_data']
            v_err = result['v_err']
            v_model = result['v_model']
            
            # Plot data with error bars
            ax.errorbar(r_data, v_data, yerr=v_err, fmt='o', 
                       color='blue', alpha=0.7, label='SPARC Data')
            
            # Plot UDT model
            ax.plot(r_data, v_model, 'r-', linewidth=2,
                   label=f'UDT (Power={enhancement_power}, R0={result["R0_fit"]:.1f})')
            
            ax.set_xlabel('Radius (kpc)')
            ax.set_ylabel('Velocity (km/s)')
            ax.set_title(f'{result["name"]}\\nchi2/nu={result["chi2_per_dof"]:.2f}')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        plt.suptitle(f'SPARC Clean Analysis: No Distance Contamination (Power={enhancement_power})', fontsize=14)
        plt.tight_layout()
        plt.savefig('C:/UDT/results/sparc_clean_no_distances.png', 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Saved: C:/UDT/results/sparc_clean_no_distances.png")
    
    def run_complete_clean_analysis(self):
        """
        Run complete clean analysis without distance contamination.
        """
        print("COMPLETE CLEAN ANALYSIS - NO DISTANCE CONTAMINATION")
        print("=" * 55)
        print()
        
        # Load pure rotation curves
        galaxies = self.load_pure_rotation_curves()
        
        if not galaxies:
            print("CRITICAL ERROR: No galaxies loaded!")
            return {'conclusion': 'DATA_LOADING_FAILURE'}
        
        # Test both enhancement formulations
        results_power2, results_power3 = self.test_both_enhancements(galaxies, max_galaxies=15)
        
        # Compare and select best approach
        best_results, best_power = self.compare_enhancement_results(results_power2, results_power3)
        
        # Create visualization
        if best_results:
            self.create_clean_visualization(best_results, best_power)
        
        # Final assessment
        print("\\n" + "=" * 60)
        print("FINAL ASSESSMENT: CLEAN ANALYSIS WITHOUT DISTANCE CONTAMINATION")
        print("=" * 60)
        
        if best_results:
            chi2_values = [r['chi2_per_dof'] for r in best_results]
            mean_chi2 = np.mean(chi2_values)
            
            print(f"BEST APPROACH: Power={best_power} enhancement")
            print(f"Success rate: {len(best_results)}/15")
            print(f"Mean chi2/DOF: {mean_chi2:.2f}")
            
            if mean_chi2 < 5.0:
                conclusion = "SUCCESS: Restored good fits without distance contamination"
            elif mean_chi2 < 20.0:
                conclusion = "PROMISING: Significant improvement"
            else:
                conclusion = "NEEDS_WORK: Further investigation needed"
            
            print(f"Conclusion: {conclusion}")
            
            return {
                'results': best_results,
                'enhancement_power': best_power,
                'mean_chi2_per_dof': mean_chi2,
                'conclusion': conclusion
            }
        else:
            return {'conclusion': 'FAILURE'}

def main():
    """Run clean analysis without distance contamination."""
    analyzer = SPARCCleanNoDistances()
    results = analyzer.run_complete_clean_analysis()
    return results

if __name__ == "__main__":
    main()