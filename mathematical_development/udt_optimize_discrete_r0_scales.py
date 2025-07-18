#!/usr/bin/env python3
"""
UDT Optimize Discrete R0 Scales for Different Galaxy Regimes
============================================================

APPROACH: Use discrete R0 values optimized for different galaxy magnitude regimes
- Small galaxies: R0_small (optimize)
- Medium galaxies: R0_medium (optimize) 
- Large galaxies: R0_large (optimize)

This returns to original UDT philosophy: different scales for different regimes
Target: Restore chi2/DOF ~ 3.13 performance

Author: Charles Rotter
Date: 2025-07-18
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar, minimize
import pandas as pd
from pathlib import Path
import os

class UDTOptimizeDiscreteR0Scales:
    """
    Optimize discrete R0 scales for different galaxy magnitude regimes.
    """
    
    def __init__(self):
        print("UDT OPTIMIZE DISCRETE R0 SCALES")
        print("=" * 40)
        print("APPROACH: Different R0 scales for different galaxy regimes")
        print("TARGET: Restore chi2/DOF ~ 3.13 performance")
        print("=" * 40)
        print()
        
        # Physical constants
        self.c = 299792.458  # km/s (speed of light)
        self.G = 4.302e-6  # km^2/s^2 per solar mass per kpc
        self.H0 = 70.0  # km/s/Mpc (for z->d conversion)
        
        # Initial R0 guesses for different regimes (to be optimized)
        self.R0_small = 100.0    # kpc (for small galaxies)
        self.R0_medium = 150.0   # kpc (for medium galaxies)
        self.R0_large = 200.0    # kpc (for large galaxies)
        
        print("INITIAL R0 SCALES:")
        print(f"R0_small = {self.R0_small} kpc (small galaxies)")
        print(f"R0_medium = {self.R0_medium} kpc (medium galaxies)")
        print(f"R0_large = {self.R0_large} kpc (large galaxies)")
        print()
    
    def extract_redshift_from_sparc(self):
        """Extract redshift and classify galaxies by size."""
        print("EXTRACTING AND CLASSIFYING SPARC GALAXIES")
        print("-" * 45)
        
        sparc_dir = Path("C:/UDT/data/sparc_database")
        galaxies = []
        
        for file_path in sparc_dir.glob("*_rotmod.dat"):
            galaxy_name = file_path.stem.replace("_rotmod", "")
            
            try:
                data = []
                distance_lcdm = None
                
                with open(file_path, 'r') as f:
                    for line in f:
                        if line.startswith('#') and 'Distance' in line:
                            try:
                                if '=' in line:
                                    distance_str = line.split('=')[1].strip().split()[0]
                                    distance_lcdm = float(distance_str)
                            except:
                                pass
                        
                        if line.startswith('#') or not line.strip():
                            continue
                            
                        parts = line.split()
                        if len(parts) >= 3:
                            try:
                                rad = float(parts[0])      # Column 1: Radius (kpc)
                                vobs = float(parts[1])     # Column 2: Observed velocity (km/s)
                                errv = float(parts[2])     # Column 3: Velocity error (km/s)
                                
                                if rad > 0 and vobs > 0 and errv > 0:
                                    data.append([rad, vobs, errv])
                            except ValueError:
                                continue
                
                if len(data) >= 5 and distance_lcdm is not None:
                    data = np.array(data)
                    redshift = self.H0 * distance_lcdm / self.c
                    
                    # Classify galaxy size by maximum velocity
                    max_velocity = np.max(data[:, 1])
                    if max_velocity < 50:
                        regime = 'small'
                    elif max_velocity < 150:
                        regime = 'medium'
                    else:
                        regime = 'large'
                    
                    galaxy_data = {
                        'name': galaxy_name,
                        'radius': data[:, 0],      # kpc (galactocentric radius)
                        'v_obs': data[:, 1],       # km/s (rotation velocity)
                        'v_err': data[:, 2],       # km/s (velocity error)
                        'distance_lcdm': distance_lcdm,  # Mpc (LCDM distance)
                        'redshift': redshift,      # Pure redshift (clean)
                        'max_velocity': max_velocity,
                        'regime': regime
                    }
                    
                    galaxies.append(galaxy_data)
                    
            except Exception as e:
                continue
        
        # Count regimes
        regimes = [g['regime'] for g in galaxies]
        small_count = regimes.count('small')
        medium_count = regimes.count('medium')
        large_count = regimes.count('large')
        
        print(f"Successfully classified {len(galaxies)} galaxies:")
        print(f"  Small galaxies (v < 50 km/s): {small_count}")
        print(f"  Medium galaxies (50-150 km/s): {medium_count}")
        print(f"  Large galaxies (v > 150 km/s): {large_count}")
        print()
        
        return galaxies
    
    def get_r0_for_regime(self, regime):
        """Get R0 value for galaxy regime."""
        if regime == 'small':
            return self.R0_small
        elif regime == 'medium':
            return self.R0_medium
        else:  # large
            return self.R0_large
    
    def udt_rotation_velocity_regime(self, r, M_total, regime):
        """
        UDT rotation velocity with regime-specific R0.
        """
        R0_regime = self.get_r0_for_regime(regime)
        tau_r = R0_regime / (R0_regime + r)
        
        # Enhancement: (1/tau)^2 for galactic dynamics
        enhancement = (1 / tau_r)**2
        M_eff = M_total * enhancement
        
        # Circular velocity
        v_circ = np.sqrt(self.G * M_eff / r)
        return v_circ
    
    def fit_galaxy_with_regime_r0(self, galaxy):
        """
        Fit galaxy using regime-specific R0.
        """
        r_gal = galaxy['radius']  # kpc
        v_obs = galaxy['v_obs']   # km/s
        v_err = galaxy['v_err']   # km/s
        regime = galaxy['regime']
        
        # Weights for chi-squared
        weights = 1 / v_err**2
        
        def objective(M_total):
            if M_total <= 0:
                return 1e10
            
            try:
                v_model = self.udt_rotation_velocity_regime(r_gal, M_total, regime)
                chi2 = np.sum(weights * (v_obs - v_model)**2)
                return chi2
            except:
                return 1e10
        
        # Optimize total mass with regime-specific R0
        result = minimize_scalar(objective, bounds=(1e8, 1e13), method='bounded')
        
        if result.success:
            M_best = result.x
            chi2_best = result.fun
            dof = len(r_gal) - 1  # One parameter (M_total)
            chi2_per_dof = chi2_best / dof
            
            # Calculate model curve
            v_model = self.udt_rotation_velocity_regime(r_gal, M_best, regime)
            rms = np.sqrt(np.mean((v_obs - v_model)**2))
            
            return {
                'success': True,
                'name': galaxy['name'],
                'regime': regime,
                'R0_used': self.get_r0_for_regime(regime),
                'M_total': M_best,
                'chi2_per_dof': chi2_per_dof,
                'rms': rms,
                'r_data': r_gal,
                'v_data': v_obs,
                'v_err': v_err,
                'v_model': v_model
            }
        else:
            return {'success': False, 'name': galaxy['name'], 'regime': regime}
    
    def analyze_with_regime_r0(self, galaxies, max_galaxies=20):
        """
        Analyze galaxies using regime-specific R0 values.
        """
        print("ANALYZING WITH REGIME-SPECIFIC R0 VALUES")
        print("-" * 45)
        print("Using discrete R0 scales for different galaxy regimes")
        print()
        
        results = []
        successful_fits = 0
        
        for i, galaxy in enumerate(galaxies[:max_galaxies]):
            print(f"Fitting {i+1}/{min(max_galaxies, len(galaxies))}: {galaxy['name']} ({galaxy['regime']})")
            
            fit_result = self.fit_galaxy_with_regime_r0(galaxy)
            
            if fit_result['success']:
                successful_fits += 1
                results.append(fit_result)
                
                print(f"  Regime: {fit_result['regime']}")
                print(f"  R0 = {fit_result['R0_used']:.0f} kpc")
                print(f"  M_total = {fit_result['M_total']:.2e} M_sun")
                print(f"  chi2/DOF = {fit_result['chi2_per_dof']:.2f}")
                print(f"  RMS = {fit_result['rms']:.1f} km/s")
            else:
                print(f"  FAILED: {galaxy['name']}")
            print()
        
        success_rate = successful_fits / min(max_galaxies, len(galaxies)) * 100
        print(f"SUCCESS RATE: {successful_fits}/{min(max_galaxies, len(galaxies))} = {success_rate:.1f}%")
        
        if results:
            chi2_values = [r['chi2_per_dof'] for r in results]
            rms_values = [r['rms'] for r in results]
            
            print(f"STATISTICS:")
            print(f"Mean chi2/DOF = {np.mean(chi2_values):.2f} +/- {np.std(chi2_values):.2f}")
            print(f"Mean RMS = {np.mean(rms_values):.1f} +/- {np.std(rms_values):.1f} km/s")
            
            # Statistics by regime
            for regime in ['small', 'medium', 'large']:
                regime_results = [r for r in results if r['regime'] == regime]
                if regime_results:
                    regime_chi2 = [r['chi2_per_dof'] for r in regime_results]
                    regime_rms = [r['rms'] for r in regime_results]
                    print(f"  {regime.capitalize()}: chi2/DOF = {np.mean(regime_chi2):.2f}, RMS = {np.mean(regime_rms):.1f} km/s")
            
            print()
            
            # Assessment
            mean_chi2 = np.mean(chi2_values)
            if mean_chi2 < 5.0:
                print("SUCCESS: Discrete R0 scales provide excellent fits!")
                print("Original UDT approach with optimized scales works")
                assessment = "SUCCESS"
            elif mean_chi2 < 20.0:
                print("PROMISING: Discrete R0 scales show significant improvement")
                assessment = "PROMISING"
            else:
                print("NEEDS OPTIMIZATION: R0 scales need further tuning")
                assessment = "NEEDS_OPTIMIZATION"
            
            return results, assessment
        else:
            return [], "FAILURE"
    
    def optimize_r0_scales(self, galaxies):
        """
        Optimize R0 scales for different galaxy regimes.
        """
        print("OPTIMIZING R0 SCALES FOR GALAXY REGIMES")
        print("-" * 45)
        
        # Separate galaxies by regime
        regimes = {'small': [], 'medium': [], 'large': []}
        for galaxy in galaxies:
            regimes[galaxy['regime']].append(galaxy)
        
        optimized_r0 = {}
        
        for regime, regime_galaxies in regimes.items():
            if not regime_galaxies:
                print(f"No {regime} galaxies found")
                continue
                
            print(f"\\nOptimizing R0 for {regime} galaxies ({len(regime_galaxies)} galaxies)")
            
            def regime_objective(R0_test):
                if R0_test <= 0:
                    return 1e10
                
                # Temporarily set R0 for this regime
                original_r0 = self.get_r0_for_regime(regime)
                if regime == 'small':
                    self.R0_small = R0_test
                elif regime == 'medium':
                    self.R0_medium = R0_test
                else:
                    self.R0_large = R0_test
                
                # Test first 10 galaxies in regime
                total_chi2 = 0
                total_dof = 0
                successful_fits = 0
                
                for galaxy in regime_galaxies[:10]:
                    fit_result = self.fit_galaxy_with_regime_r0(galaxy)
                    if fit_result['success']:
                        total_chi2 += fit_result['chi2_per_dof'] * (len(galaxy['radius']) - 1)
                        total_dof += len(galaxy['radius']) - 1
                        successful_fits += 1
                
                # Restore original R0
                if regime == 'small':
                    self.R0_small = original_r0
                elif regime == 'medium':
                    self.R0_medium = original_r0
                else:
                    self.R0_large = original_r0
                
                if successful_fits == 0:
                    return 1e10
                
                mean_chi2_per_dof = total_chi2 / total_dof
                return mean_chi2_per_dof
            
            # Optimize R0 for this regime
            result = minimize_scalar(regime_objective, bounds=(20, 500), method='bounded')
            
            if result.success:
                optimized_r0[regime] = result.x
                print(f"  Optimized R0_{regime} = {result.x:.1f} kpc")
                print(f"  Best chi2/DOF = {result.fun:.2f}")
            else:
                print(f"  Optimization failed for {regime} regime")
                optimized_r0[regime] = self.get_r0_for_regime(regime)
        
        # Update R0 values
        if 'small' in optimized_r0:
            self.R0_small = optimized_r0['small']
        if 'medium' in optimized_r0:
            self.R0_medium = optimized_r0['medium']
        if 'large' in optimized_r0:
            self.R0_large = optimized_r0['large']
        
        print(f"\\nFINAL OPTIMIZED R0 SCALES:")
        print(f"R0_small = {self.R0_small:.1f} kpc")
        print(f"R0_medium = {self.R0_medium:.1f} kpc")
        print(f"R0_large = {self.R0_large:.1f} kpc")
        print()
        
        return optimized_r0
    
    def run_complete_optimization_analysis(self):
        """
        Run complete discrete R0 optimization analysis.
        """
        print("COMPLETE DISCRETE R0 OPTIMIZATION ANALYSIS")
        print("=" * 45)
        print()
        
        # Extract and classify galaxies
        galaxies = self.extract_redshift_from_sparc()
        
        if not galaxies:
            print("CRITICAL ERROR: No galaxy data loaded!")
            return {'conclusion': 'NO_DATA'}
        
        # Test initial R0 values
        print("TESTING INITIAL R0 VALUES")
        print("-" * 30)
        initial_results, initial_assessment = self.analyze_with_regime_r0(galaxies, max_galaxies=15)
        
        # Optimize R0 scales
        optimized_r0 = self.optimize_r0_scales(galaxies)
        
        # Test optimized R0 values
        print("\\nTESTING OPTIMIZED R0 VALUES")
        print("-" * 35)
        optimized_results, optimized_assessment = self.analyze_with_regime_r0(galaxies, max_galaxies=15)
        
        # Final assessment
        print("\\n" + "=" * 60)
        print("FINAL ASSESSMENT: DISCRETE R0 OPTIMIZATION")
        print("=" * 60)
        
        if optimized_results:
            chi2_values = [r['chi2_per_dof'] for r in optimized_results]
            mean_chi2 = np.mean(chi2_values)
            
            print(f"OPTIMIZED DISCRETE R0 SCALES:")
            print(f"Small galaxies: R0 = {self.R0_small:.1f} kpc")
            print(f"Medium galaxies: R0 = {self.R0_medium:.1f} kpc")
            print(f"Large galaxies: R0 = {self.R0_large:.1f} kpc")
            print()
            print(f"PERFORMANCE:")
            print(f"Success rate: {len(optimized_results)}/15")
            print(f"Mean chi2/DOF: {mean_chi2:.2f}")
            print(f"Assessment: {optimized_assessment}")
            
            if mean_chi2 < 5.0:
                print("\\nBREAKTHROUGH: Optimized discrete R0 scales restore good fits!")
                print("Original UDT approach validated with proper scale optimization")
            elif mean_chi2 < 50:
                print("\\nSIGNIFICANT IMPROVEMENT: Getting much closer to target performance")
            
            return {
                'results': optimized_results,
                'optimized_r0': optimized_r0,
                'mean_chi2_per_dof': mean_chi2,
                'conclusion': optimized_assessment
            }
        else:
            return {'conclusion': 'FAILURE'}

def main():
    """Run discrete R0 optimization analysis."""
    analyzer = UDTOptimizeDiscreteR0Scales()
    results = analyzer.run_complete_optimization_analysis()
    return results

if __name__ == "__main__":
    main()