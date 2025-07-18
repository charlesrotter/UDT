#!/usr/bin/env python3
"""
SPARC Analysis with Correct Enhancement Formula
==============================================

DEBUGGING: Test different enhancement factor formulations to match previous success
- Enhancement Option 1: (1 + r/R0)^2 (from CLAUDE.md galactic dynamics)
- Enhancement Option 2: (1 + r/R0)^3 (from mass enhancement formula)
- Find which formulation gave chi2/DOF ~ 3.13 results

Author: Charles Rotter
Date: 2025-07-18
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
import pandas as pd
from pathlib import Path
import os

class SPARCEnhancementFormulaTesting:
    """
    Test different enhancement formulations to match previous successful results.
    """
    
    def __init__(self):
        print("SPARC ENHANCEMENT FORMULA TESTING")
        print("=" * 40)
        print("DEBUGGING: Find correct enhancement factor that gave good results")
        print("Testing: (1+r/R0)^2 vs (1+r/R0)^3 formulations")
        print("=" * 40)
        print()
        
        # Physical constants
        self.c = 299792.458  # km/s
        self.G = 4.302e-6  # km²/s² per solar mass per kpc
        
        # Use constant R0 for galactic scale (per corrected approach)
        self.R0_galactic = 38.0  # kpc
        
        print("TESTING PARAMETERS:")
        print(f"R0_galactic = {self.R0_galactic} kpc (constant)")
        print(f"G = {self.G} km²/s² per M_sun per kpc")
        print()
    
    def load_clean_sparc_data(self):
        """Load ONLY clean observational SPARC data."""
        print("LOADING CLEAN SPARC DATA")
        print("-" * 25)
        
        sparc_dir = Path("C:/UDT/data/sparc_database")
        galaxies = []
        
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
                    galaxies.append({
                        'name': galaxy_name,
                        'radius': data[:, 0],
                        'v_obs': data[:, 1],
                        'v_err': data[:, 2]
                    })
                    
            except Exception as e:
                continue
        
        print(f"Successfully loaded {len(galaxies)} galaxies")
        return galaxies
    
    def udt_rotation_velocity_enhancement2(self, r, M_total):
        """
        UDT rotation velocity with (1 + r/R0)^2 enhancement.
        """
        R0 = self.R0_galactic
        tau_r = R0 / (R0 + r)
        
        # Enhancement: (1/tau)^2 = (1 + r/R0)^2
        enhancement = (1 / tau_r)**2
        M_eff = M_total * enhancement
        
        # Circular velocity
        v_circ = np.sqrt(self.G * M_eff / r)
        return v_circ
    
    def udt_rotation_velocity_enhancement3(self, r, M_total):
        """
        UDT rotation velocity with (1 + r/R0)^3 enhancement.
        """
        R0 = self.R0_galactic
        tau_r = R0 / (R0 + r)
        
        # Enhancement: (1/tau)^3 = (1 + r/R0)^3
        enhancement = (1 / tau_r)**3
        M_eff = M_total * enhancement
        
        # Circular velocity
        v_circ = np.sqrt(self.G * M_eff / r)
        return v_circ
    
    def fit_galaxy_with_enhancement(self, galaxy, enhancement_power=2):
        """
        Fit UDT model to galaxy with specified enhancement power.
        """
        r_data = galaxy['radius']  # kpc
        v_data = galaxy['v_obs']   # km/s
        v_err = galaxy['v_err']    # km/s
        
        # Weights for chi-squared
        weights = 1 / v_err**2
        
        def objective(M_total):
            if M_total <= 0:
                return 1e10
            
            try:
                if enhancement_power == 2:
                    v_model = self.udt_rotation_velocity_enhancement2(r_data, M_total)
                else:
                    v_model = self.udt_rotation_velocity_enhancement3(r_data, M_total)
                chi2 = np.sum(weights * (v_data - v_model)**2)
                return chi2
            except:
                return 1e10
        
        # Optimize total mass only
        result = minimize_scalar(objective, bounds=(1e8, 1e13), method='bounded')
        
        if result.success:
            M_best = result.x
            chi2_best = result.fun
            dof = len(r_data) - 1  # One parameter (M_total)
            chi2_per_dof = chi2_best / dof
            
            # Calculate model curve
            if enhancement_power == 2:
                v_model = self.udt_rotation_velocity_enhancement2(r_data, M_best)
            else:
                v_model = self.udt_rotation_velocity_enhancement3(r_data, M_best)
            rms = np.sqrt(np.mean((v_data - v_model)**2))
            
            return {
                'success': True,
                'name': galaxy['name'],
                'M_total': M_best,
                'enhancement_power': enhancement_power,
                'chi2_per_dof': chi2_per_dof,
                'rms': rms,
                'r_data': r_data,
                'v_data': v_data,
                'v_err': v_err,
                'v_model': v_model
            }
        else:
            return {'success': False, 'name': galaxy['name']}
    
    def test_enhancement_formulations(self, galaxies, max_galaxies=10):
        """
        Test both enhancement formulations on sample galaxies.
        """
        print("TESTING ENHANCEMENT FORMULATIONS")
        print("-" * 35)
        print("Comparing (1+r/R0)^2 vs (1+r/R0)^3 enhancements")
        print()
        
        results_power2 = []
        results_power3 = []
        
        for i, galaxy in enumerate(galaxies[:max_galaxies]):
            print(f"Testing galaxy {i+1}/{min(max_galaxies, len(galaxies))}: {galaxy['name']}")
            
            # Test power=2 enhancement
            result2 = self.fit_galaxy_with_enhancement(galaxy, enhancement_power=2)
            if result2['success']:
                results_power2.append(result2)
                print(f"  Power=2: chi2/DOF = {result2['chi2_per_dof']:.2f}, RMS = {result2['rms']:.1f} km/s")
            else:
                print(f"  Power=2: FAILED")
            
            # Test power=3 enhancement
            result3 = self.fit_galaxy_with_enhancement(galaxy, enhancement_power=3)
            if result3['success']:
                results_power3.append(result3)
                print(f"  Power=3: chi2/DOF = {result3['chi2_per_dof']:.2f}, RMS = {result3['rms']:.1f} km/s")
            else:
                print(f"  Power=3: FAILED")
            print()
        
        return results_power2, results_power3
    
    def compare_enhancement_results(self, results_power2, results_power3):
        """
        Compare results from different enhancement formulations.
        """
        print("ENHANCEMENT FORMULATION COMPARISON")
        print("=" * 40)
        
        if results_power2:
            chi2_power2 = [r['chi2_per_dof'] for r in results_power2]
            rms_power2 = [r['rms'] for r in results_power2]
            
            print(f"POWER=2 ENHANCEMENT: (1 + r/R0)^2")
            print(f"Success rate: {len(results_power2)}/10")
            print(f"Mean chi2/DOF: {np.mean(chi2_power2):.2f} +/- {np.std(chi2_power2):.2f}")
            print(f"Mean RMS: {np.mean(rms_power2):.1f} +/- {np.std(rms_power2):.1f} km/s")
            print()
        
        if results_power3:
            chi2_power3 = [r['chi2_per_dof'] for r in results_power3]
            rms_power3 = [r['rms'] for r in results_power3]
            
            print(f"POWER=3 ENHANCEMENT: (1 + r/R0)^3")
            print(f"Success rate: {len(results_power3)}/10")
            print(f"Mean chi2/DOF: {np.mean(chi2_power3):.2f} +/- {np.std(chi2_power3):.2f}")
            print(f"Mean RMS: {np.mean(rms_power3):.1f} +/- {np.std(rms_power3):.1f} km/s")
            print()
        
        # Determine which is closer to target
        target_chi2 = 3.13
        target_rms = 31.0
        
        print("COMPARISON WITH TARGET:")
        print(f"Target: chi2/DOF ~ {target_chi2}, RMS ~ {target_rms} km/s")
        
        if results_power2:
            diff2_chi2 = abs(np.mean(chi2_power2) - target_chi2)
            diff2_rms = abs(np.mean(rms_power2) - target_rms)
            print(f"Power=2: |Delta chi2| = {diff2_chi2:.2f}, |Delta RMS| = {diff2_rms:.1f}")
        
        if results_power3:
            diff3_chi2 = abs(np.mean(chi2_power3) - target_chi2)
            diff3_rms = abs(np.mean(rms_power3) - target_rms)
            print(f"Power=3: |Delta chi2| = {diff3_chi2:.2f}, |Delta RMS| = {diff3_rms:.1f}")
        
        print()
        
        # Conclusion
        if results_power2 and results_power3:
            mean_chi2_2 = np.mean(chi2_power2)
            mean_chi2_3 = np.mean(chi2_power3)
            
            if mean_chi2_2 < mean_chi2_3 and mean_chi2_2 < 10:
                print("CONCLUSION: Power=2 enhancement performs better")
                print("This matches CLAUDE.md: 'Enhancement: 1/tau^2 = (1 + r/R0)^2 for galactic dynamics'")
                recommended = 2
            elif mean_chi2_3 < mean_chi2_2 and mean_chi2_3 < 10:
                print("CONCLUSION: Power=3 enhancement performs better")
                print("This matches mass enhancement: rho_eff = rho_matter * (1 + r/R0)^3")
                recommended = 3
            else:
                print("CONCLUSION: Both formulations show poor fits")
                print("Issue may be elsewhere in the implementation")
                recommended = None
            
            return recommended
        else:
            print("FAILURE: Could not test enhancement formulations")
            return None
    
    def run_enhancement_testing(self):
        """
        Run complete enhancement formulation testing.
        """
        print("COMPLETE ENHANCEMENT FORMULATION TESTING")
        print("=" * 45)
        print()
        
        # Load data
        galaxies = self.load_clean_sparc_data()
        
        if not galaxies:
            print("ERROR: No clean SPARC data loaded")
            return None
        
        # Test both formulations
        results_power2, results_power3 = self.test_enhancement_formulations(galaxies, max_galaxies=10)
        
        # Compare results
        recommended = self.compare_enhancement_results(results_power2, results_power3)
        
        return {
            'results_power2': results_power2,
            'results_power3': results_power3,
            'recommended_power': recommended
        }

def main():
    """Run enhancement formulation testing."""
    tester = SPARCEnhancementFormulaTesting()
    results = tester.run_enhancement_testing()
    return results

if __name__ == "__main__":
    main()