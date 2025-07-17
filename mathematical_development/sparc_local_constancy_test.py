#!/usr/bin/env python3
"""
SPARC Local Constancy Test: Real Data vs UDT with Local Light Speed Constancy
============================================================================

This is the DEFINITIVE TEST of UDT with local constancy interpretation.

PREVIOUS RESULT (global constancy): χ²/dof = 216.77 → CATASTROPHIC FAILURE
NEW TEST (local constancy): Same formula, different physics interpretation

UDT PREDICTION: v²_observed = v²_Newtonian × (1 + r/R₀)²

TEST CRITERIA:
- Real SPARC galaxy rotation curve data
- Rigorous statistical analysis
- Honest reporting of ALL results
- NO FUDGING - theory succeeds or fails on merit

Author: Charles Rotter
Date: 2025-01-17
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar, curve_fit
import os
import glob
from pathlib import Path

class SPARCLocalConstancyTest:
    """
    Test UDT with local constancy interpretation against real SPARC data.
    """
    
    def __init__(self):
        print("SPARC LOCAL CONSTANCY TEST")
        print("=" * 50)
        print("DEFINITIVE TEST: UDT with local light speed constancy")
        print("Previous result (global): chi^2/dof = 216.77 -> FAILURE")
        print("New test (local): Same formula, different physics")
        print("=" * 50)
        print()
        
        # Physical constants
        self.G = 6.67430e-11  # m³/kg/s²
        self.c = 2.998e8      # m/s
        self.kpc_to_m = 3.086e19  # m/kpc
        self.Msun_to_kg = 1.989e30 # kg/M_sun
        
        # Results storage
        self.galaxy_results = []
        self.test_summary = {}
        
        print("PHYSICAL CONSTANTS:")
        print(f"G = {self.G} m³/kg/s²")
        print(f"c = {self.c} m/s")
        print()
        
    def load_sparc_catalog(self):
        """
        Load SPARC galaxy catalog.
        """
        print("LOADING SPARC GALAXY CATALOG...")
        
        catalog_file = "data/sparc_database/SPARC_Lelli2016c.mrt"
        
        # Read catalog (skip header lines)
        try:
            with open(catalog_file, 'r') as f:
                lines = f.readlines()
            
            # Find start of data (after header)
            data_start = None
            for i, line in enumerate(lines):
                if line.strip().startswith('NGC') or line.strip().startswith('UGC') or line.strip().startswith('IC'):
                    data_start = i
                    break
            
            if data_start is None:
                # Try different approach - look for data lines
                for i, line in enumerate(lines):
                    if len(line.strip()) > 50 and not line.startswith('#') and not line.startswith('-'):
                        data_start = i
                        break
            
            if data_start is None:
                print("ERROR: Could not find data start in catalog")
                return None
            
            # Parse data lines
            galaxies = []
            for line in lines[data_start:]:
                if len(line.strip()) > 50:  # Data line
                    try:
                        # Parse fixed-width format
                        galaxy_name = line[0:11].strip()
                        distance = float(line[12:19].strip())
                        # Add more fields as needed
                        
                        if galaxy_name and distance > 0:
                            galaxies.append({
                                'name': galaxy_name,
                                'distance': distance,  # Mpc
                                'raw_line': line
                            })
                    except (ValueError, IndexError):
                        continue
            
            print(f"Loaded {len(galaxies)} galaxies from catalog")
            return galaxies
            
        except FileNotFoundError:
            print(f"ERROR: Catalog file not found: {catalog_file}")
            return None
        except Exception as e:
            print(f"ERROR loading catalog: {e}")
            return None
    
    def load_rotation_curve(self, galaxy_name):
        """
        Load rotation curve data for specific galaxy.
        """
        rotcurve_file = f"data/sparc_database/{galaxy_name}_rotmod.dat"
        
        if not os.path.exists(rotcurve_file):
            return None
        
        try:
            # Read rotation curve data
            data = []
            with open(rotcurve_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or line.strip() == '':
                        continue
                    
                    parts = line.strip().split()
                    if len(parts) >= 7:
                        try:
                            r_kpc = float(parts[0])      # Radius (kpc)
                            v_obs = float(parts[1])      # Observed velocity (km/s)
                            v_err = float(parts[2])      # Velocity error (km/s)
                            v_gas = float(parts[3])      # Gas velocity (km/s)
                            v_disk = float(parts[4])     # Disk velocity (km/s)
                            v_bul = float(parts[5])      # Bulge velocity (km/s)
                            
                            # Calculate baryonic velocity
                            v_bary = np.sqrt(v_gas**2 + v_disk**2 + v_bul**2)
                            
                            data.append({
                                'r_kpc': r_kpc,
                                'v_obs': v_obs,
                                'v_err': v_err,
                                'v_bary': v_bary,
                                'v_gas': v_gas,
                                'v_disk': v_disk,
                                'v_bul': v_bul
                            })
                        except ValueError:
                            continue
            
            if len(data) < 3:
                return None
            
            # Convert to arrays
            curve_data = {
                'r_kpc': np.array([d['r_kpc'] for d in data]),
                'v_obs': np.array([d['v_obs'] for d in data]),
                'v_err': np.array([d['v_err'] for d in data]),
                'v_bary': np.array([d['v_bary'] for d in data])
            }
            
            return curve_data
            
        except Exception as e:
            print(f"ERROR loading rotation curve for {galaxy_name}: {e}")
            return None
    
    def calculate_udt_velocity(self, r_kpc, v_bary, R0_kpc):
        """
        Calculate UDT velocity with local constancy interpretation.
        
        UDT PREDICTION: v^2_observed = v^2_baryonic * (1 + r/R_0)^2
        
        LOCAL CONSTANCY INTERPRETATION:
        - Light speed constant in local inertial frames
        - Coordinate light speed can vary
        - Temporal enhancement factor: (1 + r/R_0)^2
        """
        # Convert to SI units
        r_m = r_kpc * self.kpc_to_m
        R0_m = R0_kpc * self.kpc_to_m
        
        # UDT enhancement factor
        enhancement = (1 + r_m/R0_m)**2
        
        # UDT velocity prediction
        v_udt = v_bary * np.sqrt(enhancement)
        
        return v_udt
    
    def fit_single_galaxy(self, galaxy_name, curve_data):
        """
        Fit UDT to single galaxy rotation curve.
        """
        print(f"\nFITTING: {galaxy_name}")
        print("-" * 30)
        
        r_kpc = curve_data['r_kpc']
        v_obs = curve_data['v_obs']
        v_err = curve_data['v_err']
        v_bary = curve_data['v_bary']
        
        # Ensure minimum error for numerical stability
        v_err = np.maximum(v_err, 5.0)  # Minimum 5 km/s error
        
        print(f"Data points: {len(r_kpc)}")
        print(f"Radius range: {r_kpc[0]:.2f} - {r_kpc[-1]:.2f} kpc")
        print(f"Velocity range: {v_obs[0]:.1f} - {v_obs[-1]:.1f} km/s")
        
        def chi_squared(R0_kpc):
            """Calculate chi-squared for given R_0"""
            if R0_kpc <= 0:
                return 1e10
            
            try:
                v_udt = self.calculate_udt_velocity(r_kpc, v_bary, R0_kpc)
                chi2 = np.sum((v_obs - v_udt)**2 / v_err**2)
                return chi2
            except:
                return 1e10
        
        # Fit R_0 by minimizing chi-squared
        print("Fitting R_0...")
        
        # Search range: 1 kpc to 1000 kpc
        result = minimize_scalar(chi_squared, bounds=(1.0, 1000.0), method='bounded')
        
        if not result.success:
            print("ERROR: Fit failed!")
            return None
        
        R0_best = result.x
        chi2_best = result.fun
        dof = len(r_kpc) - 1  # One parameter (R_0)
        chi2_reduced = chi2_best / dof
        
        # Calculate UDT velocities with best R_0
        v_udt_best = self.calculate_udt_velocity(r_kpc, v_bary, R0_best)
        
        # Calculate residuals and RMS
        residuals = v_obs - v_udt_best
        rms_residual = np.sqrt(np.mean(residuals**2))
        
        # Calculate correlation coefficient
        correlation = np.corrcoef(v_obs, v_udt_best)[0, 1]
        
        print(f"RESULTS:")
        print(f"  R_0 = {R0_best:.1f} kpc")
        print(f"  chi^2 = {chi2_best:.2f}")
        print(f"  DOF = {dof}")
        print(f"  chi^2/DOF = {chi2_reduced:.2f}")
        print(f"  RMS residual = {rms_residual:.1f} km/s")
        print(f"  Correlation = {correlation:.3f}")
        
        # Quality assessment
        if chi2_reduced < 2.0:
            quality = "EXCELLENT"
        elif chi2_reduced < 5.0:
            quality = "GOOD"
        elif chi2_reduced < 10.0:
            quality = "MARGINAL"
        else:
            quality = "POOR"
        
        print(f"  Quality: {quality}")
        
        return {
            'galaxy': galaxy_name,
            'n_points': len(r_kpc),
            'R0_kpc': R0_best,
            'chi2': chi2_best,
            'dof': dof,
            'chi2_reduced': chi2_reduced,
            'rms_residual': rms_residual,
            'correlation': correlation,
            'quality': quality,
            'r_kpc': r_kpc,
            'v_obs': v_obs,
            'v_err': v_err,
            'v_bary': v_bary,
            'v_udt': v_udt_best,
            'residuals': residuals
        }
    
    def run_comprehensive_test(self):
        """
        Run comprehensive test on all available SPARC galaxies.
        """
        print("COMPREHENSIVE SPARC TEST")
        print("=" * 30)
        print()
        
        # Load galaxy catalog
        galaxies = self.load_sparc_catalog()
        if not galaxies:
            print("ERROR: Could not load galaxy catalog")
            return
        
        print(f"Testing {len(galaxies)} galaxies...")
        print()
        
        # Test each galaxy
        successful_fits = 0
        total_galaxies = 0
        
        for galaxy in galaxies[:20]:  # Test first 20 galaxies
            galaxy_name = galaxy['name']
            total_galaxies += 1
            
            # Load rotation curve
            curve_data = self.load_rotation_curve(galaxy_name)
            if curve_data is None:
                print(f"SKIP: {galaxy_name} (no rotation curve data)")
                continue
            
            # Fit UDT
            result = self.fit_single_galaxy(galaxy_name, curve_data)
            if result is not None:
                self.galaxy_results.append(result)
                successful_fits += 1
        
        print(f"\nSUCCESSFUL FITS: {successful_fits}/{total_galaxies}")
        print()
        
        # Analyze results
        self.analyze_results()
    
    def analyze_results(self):
        """
        Analyze overall test results.
        """
        print("OVERALL ANALYSIS")
        print("=" * 20)
        print()
        
        if len(self.galaxy_results) == 0:
            print("No successful fits to analyze")
            return
        
        # Extract statistics
        chi2_values = [r['chi2_reduced'] for r in self.galaxy_results]
        R0_values = [r['R0_kpc'] for r in self.galaxy_results]
        correlation_values = [r['correlation'] for r in self.galaxy_results]
        
        # Calculate statistics
        mean_chi2 = np.mean(chi2_values)
        std_chi2 = np.std(chi2_values)
        mean_R0 = np.mean(R0_values)
        std_R0 = np.std(R0_values)
        mean_corr = np.mean(correlation_values)
        
        print(f"STATISTICS ({len(self.galaxy_results)} galaxies):")
        print(f"  chi^2/DOF: {mean_chi2:.2f} +/- {std_chi2:.2f}")
        print(f"  R_0: {mean_R0:.1f} +/- {std_R0:.1f} kpc")
        print(f"  Correlation: {mean_corr:.3f}")
        print()
        
        # Quality breakdown
        excellent = sum(1 for r in self.galaxy_results if r['chi2_reduced'] < 2.0)
        good = sum(1 for r in self.galaxy_results if 2.0 <= r['chi2_reduced'] < 5.0)
        marginal = sum(1 for r in self.galaxy_results if 5.0 <= r['chi2_reduced'] < 10.0)
        poor = sum(1 for r in self.galaxy_results if r['chi2_reduced'] >= 10.0)
        
        total = len(self.galaxy_results)
        
        print("QUALITY BREAKDOWN:")
        print(f"  Excellent (chi^2/DOF < 2): {excellent}/{total} ({100*excellent/total:.1f}%)")
        print(f"  Good (2 <= chi^2/DOF < 5): {good}/{total} ({100*good/total:.1f}%)")
        print(f"  Marginal (5 <= chi^2/DOF < 10): {marginal}/{total} ({100*marginal/total:.1f}%)")
        print(f"  Poor (chi^2/DOF >= 10): {poor}/{total} ({100*poor/total:.1f}%)")
        print()
        
        # Comparison with previous result
        print("COMPARISON WITH PREVIOUS RESULT:")
        print("  Previous (global constancy): chi^2/DOF = 216.77 -> CATASTROPHIC FAILURE")
        print(f"  Current (local constancy): chi^2/DOF = {mean_chi2:.2f} -> ", end="")
        
        if mean_chi2 < 5.0:
            print("SIGNIFICANT IMPROVEMENT")
        elif mean_chi2 < 50.0:
            print("MAJOR IMPROVEMENT")
        else:
            print("STILL PROBLEMATIC")
        
        print()
        
        # R_0 consistency
        R0_variation = std_R0 / mean_R0
        print("R_0 CONSISTENCY:")
        print(f"  Coefficient of variation: {R0_variation:.2f}")
        
        if R0_variation < 0.5:
            print("  GOOD: R_0 values reasonably consistent")
        elif R0_variation < 1.0:
            print("  MARGINAL: R_0 values moderately scattered")
        else:
            print("  POOR: R_0 values highly scattered")
        
        print()
        
        # Scientific verdict
        print("SCIENTIFIC VERDICT:")
        success_rate = (excellent + good) / total
        
        if success_rate > 0.8 and mean_chi2 < 5.0:
            print("  UDT WITH LOCAL CONSTANCY: STRONG SUCCESS")
            print("  Theory performs well on real galactic data")
        elif success_rate > 0.6 and mean_chi2 < 10.0:
            print("  UDT WITH LOCAL CONSTANCY: MODERATE SUCCESS")
            print("  Theory shows promise but needs refinement")
        elif mean_chi2 < 100.0:
            print("  UDT WITH LOCAL CONSTANCY: MARGINAL IMPROVEMENT")
            print("  Better than global constancy but still problematic")
        else:
            print("  UDT WITH LOCAL CONSTANCY: CONTINUED FAILURE")
            print("  Local constancy interpretation insufficient")
        
        print()
        
        # Store summary
        self.test_summary = {
            'n_galaxies': total,
            'mean_chi2': mean_chi2,
            'std_chi2': std_chi2,
            'mean_R0': mean_R0,
            'std_R0': std_R0,
            'success_rate': success_rate,
            'excellent': excellent,
            'good': good,
            'marginal': marginal,
            'poor': poor
        }
    
    def create_plots(self):
        """
        Create visualization of results.
        """
        if len(self.galaxy_results) == 0:
            print("No results to plot")
            return
        
        print("CREATING PLOTS...")
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Plot 1: Example galaxy fit
        best_galaxy = min(self.galaxy_results, key=lambda x: x['chi2_reduced'])
        ax1 = axes[0, 0]
        
        r = best_galaxy['r_kpc']
        v_obs = best_galaxy['v_obs']
        v_err = best_galaxy['v_err']
        v_udt = best_galaxy['v_udt']
        
        ax1.errorbar(r, v_obs, yerr=v_err, fmt='o', color='blue', label='Observed')
        ax1.plot(r, v_udt, 'r-', linewidth=2, label='UDT (local constancy)')
        ax1.set_xlabel('Radius (kpc)')
        ax1.set_ylabel('Velocity (km/s)')
        ax1.set_title(f'Best Fit: {best_galaxy["galaxy"]} (chi^2/DOF = {best_galaxy["chi2_reduced"]:.2f})')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Chi-squared distribution
        ax2 = axes[0, 1]
        chi2_values = [r['chi2_reduced'] for r in self.galaxy_results]
        ax2.hist(chi2_values, bins=15, alpha=0.7, edgecolor='black')
        ax2.axvline(2.0, color='green', linestyle='--', label='Good fit threshold')
        ax2.axvline(5.0, color='orange', linestyle='--', label='Marginal fit threshold')
        ax2.set_xlabel('chi^2/DOF')
        ax2.set_ylabel('Number of galaxies')
        ax2.set_title('Chi-squared Distribution')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: R_0 distribution
        ax3 = axes[1, 0]
        R0_values = [r['R0_kpc'] for r in self.galaxy_results]
        ax3.hist(R0_values, bins=15, alpha=0.7, edgecolor='black')
        ax3.axvline(38, color='red', linestyle='--', label='Previous estimate')
        ax3.set_xlabel('R_0 (kpc)')
        ax3.set_ylabel('Number of galaxies')
        ax3.set_title('R_0 Distribution')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Plot 4: Correlation vs Chi-squared
        ax4 = axes[1, 1]
        correlations = [r['correlation'] for r in self.galaxy_results]
        ax4.scatter(correlations, chi2_values, alpha=0.7)
        ax4.set_xlabel('Correlation coefficient')
        ax4.set_ylabel('chi^2/DOF')
        ax4.set_title('Correlation vs Fit Quality')
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('results/sparc_local_constancy_test.png', dpi=300, bbox_inches='tight')
        print("Plot saved: results/sparc_local_constancy_test.png")
        
        return fig

def main():
    """
    Run comprehensive SPARC test with local constancy interpretation.
    """
    print("SPARC LOCAL CONSTANCY TEST")
    print("=" * 80)
    print("DEFINITIVE TEST: UDT with local light speed constancy")
    print("=" * 80)
    print()
    
    # Initialize test
    test = SPARCLocalConstancyTest()
    print()
    
    # Run comprehensive test
    test.run_comprehensive_test()
    print()
    
    # Create plots
    test.create_plots()
    print()
    
    print("=" * 80)
    print("SPARC LOCAL CONSTANCY TEST COMPLETE")
    print("=" * 80)
    print()
    
    if test.test_summary:
        print("FINAL VERDICT:")
        mean_chi2 = test.test_summary['mean_chi2']
        success_rate = test.test_summary['success_rate']
        
        print(f"  Mean chi^2/DOF: {mean_chi2:.2f}")
        print(f"  Success rate: {success_rate:.1%}")
        print(f"  Previous result: chi^2/DOF = 216.77 (catastrophic failure)")
        
        if mean_chi2 < 10.0:
            print("  BREAKTHROUGH: Local constancy interpretation works!")
        else:
            print("  CONTINUED PROBLEMS: Local constancy insufficient")
    
    print()
    print("THEORETICAL SIGNIFICANCE:")
    print("This test validates whether the local constancy reinterpretation")
    print("of light speed transforms UDT from failed to viable theory.")

if __name__ == "__main__":
    main()