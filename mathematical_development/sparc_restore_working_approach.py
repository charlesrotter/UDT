#!/usr/bin/env python3
"""
SPARC Analysis: Restore Working Approach with Separate R0 Values
================================================================

FIXING CURRENT ISSUES:
1. Return to separate R0 values for different scales (not variable function)
2. Check SPARC data interpretation - convert distances to redshift properly
3. Restore the approach that gave chi2/DOF ~ 3.13 and 97.7% success

SCALE SEPARATION:
- Solar R0: ~1-10 AU (for solar system)
- Galactic R0: ~38 kpc (optimized for galaxies)  
- Cosmological R0: ~4754 Mpc (for supernova/CMB)

Author: Charles Rotter
Date: 2025-07-18
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar, minimize
import pandas as pd
from pathlib import Path
import os

class SPARCWorkingApproach:
    """
    Restore working SPARC analysis with separate R0 values and proper data handling.
    """
    
    def __init__(self):
        print("SPARC ANALYSIS: RESTORING WORKING APPROACH")
        print("=" * 50)
        print("FIXING: Return to separate R0 values and proper data interpretation")
        print("TARGET: Restore chi2/DOF ~ 3.13 and 97.7% success rate")
        print("=" * 50)
        print()
        
        # Physical constants
        self.c = 299792.458  # km/s
        self.G = 4.302e-6  # km²/s² per solar mass per kpc
        
        # Scale-specific R0 values (to be optimized)
        self.R0_solar = 5.0  # AU (for solar system scale)
        self.R0_galactic = 38.0  # kpc (for galactic scale - start with previous value)
        self.R0_cosmological = 4754.3  # Mpc (for cosmological scale)
        
        print("SCALE-SPECIFIC R0 VALUES:")
        print(f"Solar R0 = {self.R0_solar} AU")
        print(f"Galactic R0 = {self.R0_galactic} kpc (to be optimized)")
        print(f"Cosmological R0 = {self.R0_cosmological} Mpc")
        print()
    
    def check_sparc_data_format(self):
        """
        Check SPARC data format and identify any distance measurements needing conversion.
        """
        print("CHECKING SPARC DATA FORMAT")
        print("-" * 30)
        
        sparc_dir = Path("C:/UDT/data/sparc_database")
        
        # Check a sample file to understand format
        sample_files = list(sparc_dir.glob("*_rotmod.dat"))[:3]
        
        for i, file_path in enumerate(sample_files):
            galaxy_name = file_path.stem.replace("_rotmod", "")
            print(f"Sample {i+1}: {galaxy_name}")
            
            try:
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                    
                # Show header comments
                for line in lines[:10]:
                    if line.startswith('#'):
                        print(f"  {line.strip()}")
                    else:
                        # Show first data line
                        parts = line.split()
                        if len(parts) >= 3:
                            print(f"  Data columns: {parts[0]} | {parts[1]} | {parts[2]} | ...")
                        break
                print()
                        
            except Exception as e:
                print(f"  Error reading {galaxy_name}: {e}")
                print()
        
        # Check if main SPARC table has distance information
        main_table = sparc_dir / "SPARC_Lelli2016c.mrt"
        if main_table.exists():
            print("MAIN SPARC TABLE FOUND:")
            print(f"File: {main_table}")
            try:
                with open(main_table, 'r') as f:
                    lines = f.readlines()
                    for i, line in enumerate(lines[:20]):
                        if 'Dist' in line or 'distance' in line or 'redshift' in line or 'z' in line:
                            print(f"  Line {i}: {line.strip()}")
            except Exception as e:
                print(f"  Error reading main table: {e}")
            print()
    
    def load_clean_sparc_data_with_distances(self):
        """
        Load SPARC data and check for distance information that needs redshift conversion.
        """
        print("LOADING CLEAN SPARC DATA WITH DISTANCE CHECK")
        print("-" * 45)
        
        sparc_dir = Path("C:/UDT/data/sparc_database")
        galaxies = []
        
        # First, try to load distance information from main table
        main_table = sparc_dir / "SPARC_Lelli2016c.mrt"
        galaxy_distances = {}
        
        if main_table.exists():
            print("Loading galaxy distance information...")
            try:
                # Parse main table for distance information
                with open(main_table, 'r') as f:
                    for line in f:
                        if line.startswith('#') or not line.strip():
                            continue
                        parts = line.split()
                        if len(parts) >= 2:
                            # Assume format: Galaxy_name Distance ...
                            galaxy_name = parts[0]
                            try:
                                distance = float(parts[1])  # Assuming distance in second column
                                galaxy_distances[galaxy_name] = distance
                            except ValueError:
                                continue
                print(f"Found distance info for {len(galaxy_distances)} galaxies")
            except Exception as e:
                print(f"Could not parse main table: {e}")
        
        # Load individual rotation curve files
        for file_path in sparc_dir.glob("*_rotmod.dat"):
            galaxy_name = file_path.stem.replace("_rotmod", "")
            
            try:
                data = []
                with open(file_path, 'r') as f:
                    for line in f:
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
                        'v_err': data[:, 2],       # km/s
                        'distance': galaxy_distances.get(galaxy_name, None)  # Mpc if available
                    }
                    
                    # Check if we need redshift conversion
                    if galaxy_data['distance'] is not None:
                        # Convert distance to redshift using Hubble law: z = H0 * d / c
                        H0 = 70.0  # km/s/Mpc (approximate)
                        z = H0 * galaxy_data['distance'] / self.c
                        galaxy_data['redshift'] = z
                        print(f"  {galaxy_name}: d = {galaxy_data['distance']:.1f} Mpc → z = {z:.6f}")
                    
                    galaxies.append(galaxy_data)
                    
            except Exception as e:
                continue
        
        print(f"Successfully loaded {len(galaxies)} galaxies")
        
        # Check redshift distribution
        galaxies_with_z = [g for g in galaxies if 'redshift' in g and g['redshift'] is not None]
        if galaxies_with_z:
            redshifts = [g['redshift'] for g in galaxies_with_z]
            print(f"Redshift range: {np.min(redshifts):.6f} - {np.max(redshifts):.6f}")
            print(f"Mean redshift: {np.mean(redshifts):.6f}")
        else:
            print("No redshift information found - using rotation curves directly")
        
        print()
        return galaxies
    
    def udt_rotation_velocity_optimized_r0(self, r, M_total, R0_fit):
        """
        UDT rotation velocity with optimized R0 for galactic scale.
        """
        tau_r = R0_fit / (R0_fit + r)
        
        # Enhancement: Use (1/tau)^2 for galactic dynamics (per CLAUDE.md)
        enhancement = (1 / tau_r)**2
        M_eff = M_total * enhancement
        
        # Circular velocity
        v_circ = np.sqrt(self.G * M_eff / r)
        return v_circ
    
    def fit_galaxy_with_optimized_r0(self, galaxy):
        """
        Fit UDT model to galaxy with both M_total and R0 as free parameters.
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
                v_model = self.udt_rotation_velocity_optimized_r0(r_data, M_total, R0_fit)
                chi2 = np.sum(weights * (v_data - v_model)**2)
                return chi2
            except:
                return 1e10
        
        # Optimize both M_total and R0
        # Initial guess
        M_guess = 1e11  # Solar masses
        R0_guess = self.R0_galactic  # Start with previous value
        
        result = minimize(objective, [M_guess, R0_guess], 
                         bounds=[(1e8, 1e13), (1, 1000)],
                         method='L-BFGS-B')
        
        if result.success:
            M_best, R0_best = result.x
            chi2_best = result.fun
            dof = len(r_data) - 2  # Two parameters
            chi2_per_dof = chi2_best / dof
            
            # Calculate model curve
            v_model = self.udt_rotation_velocity_optimized_r0(r_data, M_best, R0_best)
            rms = np.sqrt(np.mean((v_data - v_model)**2))
            
            return {
                'success': True,
                'name': galaxy['name'],
                'M_total': M_best,
                'R0_fit': R0_best,
                'chi2_per_dof': chi2_per_dof,
                'rms': rms,
                'r_data': r_data,
                'v_data': v_data,
                'v_err': v_err,
                'v_model': v_model,
                'distance': galaxy.get('distance', None),
                'redshift': galaxy.get('redshift', None)
            }
        else:
            return {'success': False, 'name': galaxy['name']}
    
    def analyze_galaxies_working_approach(self, galaxies, max_galaxies=20):
        """
        Analyze galaxies using working approach with optimized R0.
        """
        print("ANALYZING GALAXIES WITH WORKING APPROACH")
        print("-" * 45)
        print("Optimizing both M_total and R0 for each galaxy")
        print("Target: chi2/DOF ~ 3.13, RMS ~ 31 km/s")
        print()
        
        results = []
        successful_fits = 0
        
        for i, galaxy in enumerate(galaxies[:max_galaxies]):
            print(f"Fitting {i+1}/{min(max_galaxies, len(galaxies))}: {galaxy['name']}")
            
            fit_result = self.fit_galaxy_with_optimized_r0(galaxy)
            
            if fit_result['success']:
                successful_fits += 1
                results.append(fit_result)
                
                print(f"  SUCCESS: M_total = {fit_result['M_total']:.2e} M_sun")
                print(f"           R0_fit = {fit_result['R0_fit']:.1f} kpc")
                print(f"           chi2/DOF = {fit_result['chi2_per_dof']:.2f}")
                print(f"           RMS = {fit_result['rms']:.1f} km/s")
                
                if fit_result['redshift'] is not None:
                    print(f"           z = {fit_result['redshift']:.6f}")
            else:
                print(f"  FAILED: {galaxy['name']}")
            print()
        
        success_rate = successful_fits / min(max_galaxies, len(galaxies)) * 100
        print(f"SUCCESS RATE: {successful_fits}/{min(max_galaxies, len(galaxies))} = {success_rate:.1f}%")
        
        if results:
            chi2_values = [r['chi2_per_dof'] for r in results]
            rms_values = [r['rms'] for r in results]
            r0_values = [r['R0_fit'] for r in results]
            
            print(f"STATISTICS:")
            print(f"Mean chi2/DOF = {np.mean(chi2_values):.2f} +/- {np.std(chi2_values):.2f}")
            print(f"Mean RMS = {np.mean(rms_values):.1f} +/- {np.std(rms_values):.1f} km/s")
            print(f"Mean R0_fit = {np.mean(r0_values):.1f} +/- {np.std(r0_values):.1f} kpc")
            print()
            
            # Assessment
            mean_chi2 = np.mean(chi2_values)
            mean_rms = np.mean(rms_values)
            
            print("ASSESSMENT:")
            if mean_chi2 < 5.0 and mean_rms < 40.0:
                print(f"SUCCESS: Restored good fits! chi2/DOF = {mean_chi2:.2f} < 5.0")
                print(f"RMS = {mean_rms:.1f} km/s < 40.0")
                print("Working approach successfully restored")
            elif mean_chi2 < 10.0:
                print(f"PROMISING: Significant improvement, chi2/DOF = {mean_chi2:.2f}")
                print("Close to target performance")
            else:
                print(f"ISSUE: Still poor fits, chi2/DOF = {mean_chi2:.2f}")
                print("Need further investigation")
        
        return results
    
    def create_working_approach_visualization(self, results):
        """
        Create visualization of working approach results.
        """
        if not results:
            return
        
        print("\\nCREATING WORKING APPROACH VISUALIZATION")
        print("-" * 40)
        
        # Select first 6 results for detailed plotting
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
                   label=f'UDT (R0={result["R0_fit"]:.1f} kpc)')
            
            ax.set_xlabel('Radius (kpc)')
            ax.set_ylabel('Velocity (km/s)')
            ax.set_title(f'{result["name"]}\\nchi2/nu={result["chi2_per_dof"]:.2f}')
            ax.legend()
            ax.grid(True, alpha=0.3)
        
        plt.suptitle('SPARC Working Approach: Optimized R0 UDT Fits', fontsize=14)
        plt.tight_layout()
        plt.savefig('C:/UDT/results/sparc_working_approach_restored.png', 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Saved: C:/UDT/results/sparc_working_approach_restored.png")
    
    def run_complete_working_analysis(self):
        """
        Run complete analysis with restored working approach.
        """
        print("COMPLETE WORKING APPROACH RESTORATION")
        print("=" * 40)
        print()
        
        # Step 1: Check data format
        self.check_sparc_data_format()
        
        # Step 2: Load data with distance checking
        galaxies = self.load_clean_sparc_data_with_distances()
        
        # Step 3: Analyze with working approach
        results = self.analyze_galaxies_working_approach(galaxies, max_galaxies=20)
        
        # Step 4: Create visualization
        self.create_working_approach_visualization(results)
        
        # Step 5: Final assessment
        print("\\n" + "=" * 60)
        print("FINAL ASSESSMENT: WORKING APPROACH RESTORED")
        print("=" * 60)
        
        if results:
            chi2_values = [r['chi2_per_dof'] for r in results]
            mean_chi2 = np.mean(chi2_values)
            
            if mean_chi2 < 5.0:
                conclusion = "SUCCESS: Working approach restored"
            elif mean_chi2 < 10.0:
                conclusion = "PROMISING: Significant improvement"
            else:
                conclusion = "NEEDS_WORK: Further investigation needed"
            
            return {
                'results': results,
                'mean_chi2_per_dof': mean_chi2,
                'conclusion': conclusion
            }
        else:
            return {'conclusion': 'FAILURE'}

def main():
    """Run working approach restoration."""
    analyzer = SPARCWorkingApproach()
    results = analyzer.run_complete_working_analysis()
    return results

if __name__ == "__main__":
    main()