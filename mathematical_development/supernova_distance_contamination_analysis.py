#!/usr/bin/env python3
"""
Supernova Distance Contamination Analysis
==========================================

CRITICAL INSIGHT: Supernova data may be contaminated by LCDM distance assumptions!

The key issue: Most supernova analyses convert redshift to distance using LCDM cosmology.
If UDT predicts a fundamentally different universe size/structure, then:
1. The "distances" in supernova catalogs are wrong for UDT
2. The magnitudes may be correct, but distances are LCDM-contaminated
3. We need to work backwards from redshift + magnitude to test UDT

This analysis investigates whether the poor supernova fits are due to:
- Fundamental problems with UDT cosmology
- OR contamination from LCDM distance assumptions

Author: Charles Rotter
Date: 2025-07-18
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
import sys
import os

# Add project root to path
sys.path.append(os.path.abspath('.'))

class SupernovaDistanceContamination:
    def __init__(self):
        print("SUPERNOVA DISTANCE CONTAMINATION ANALYSIS")
        print("=" * 45)
        
        # Physical constants
        self.c = 299792.458  # km/s
        
        # Load supernova data
        self.load_supernova_data()
        
    def load_supernova_data(self):
        """Load supernova data from previous analysis."""
        print("\nLoading supernova data...")
        
        # Load Pantheon+ data (more reliable than CSP)
        pantheon_file = "data/Pantheon_SH0ES.dat"
        try:
            # Read the data with proper column names
            self.df = pd.read_csv(pantheon_file, delim_whitespace=True, comment='#')
            print(f"Loaded {len(self.df)} supernovae from Pantheon+")
            
            # Filter to low redshift for local universe analysis
            self.df = self.df[self.df['zCMB'] <= 0.1]
            print(f"After z <= 0.1 cut: {len(self.df)} supernovae")
            print(f"Redshift range: {self.df['zCMB'].min():.4f} - {self.df['zCMB'].max():.4f}")
            
        except FileNotFoundError:
            print(f"Could not find {pantheon_file}")
            print("Using synthetic data for demonstration")
            self.create_synthetic_data()
    
    def create_synthetic_data(self):
        """Create synthetic supernova data for testing."""
        print("Creating synthetic supernova data...")
        
        # Generate synthetic redshifts
        z = np.random.uniform(0.01, 0.1, 100)
        
        # Generate magnitudes using UDT distance relation
        R0_true = 4000  # Mpc
        M_B_true = -18.6  # Absolute magnitude
        
        # UDT distance: d_L = z * R0
        d_L_udt = z * R0_true
        
        # Distance modulus: m - M = 5*log10(d_L/10 pc)
        # d_L in Mpc, so d_L/10 pc = d_L * 1e6 / 10 = d_L * 1e5
        distance_modulus = 5 * np.log10(d_L_udt * 1e5)
        mB = M_B_true + distance_modulus
        
        # Add noise
        mB += np.random.normal(0, 0.1, len(z))
        
        # Create DataFrame
        self.df = pd.DataFrame({
            'zCMB': z,
            'mB': mB,
            'mBERR': np.full(len(z), 0.1)
        })
        
        print(f"Created {len(self.df)} synthetic supernovae")
    
    def lcdm_distance(self, z, H0=70, Omega_m=0.3, Omega_Lambda=0.7):
        """Calculate LCDM luminosity distance."""
        # Simplified LCDM distance for low redshift
        # d_L ≈ (c*z/H0) * [1 + z*(1 - q0)/2]
        # where q0 = Omega_m/2 - Omega_Lambda ≈ -0.55
        
        q0 = Omega_m/2 - Omega_Lambda
        d_H = self.c / H0  # Hubble distance in Mpc
        
        # Low redshift approximation
        d_L = d_H * z * (1 + z * (1 - q0) / 2)
        
        return d_L
    
    def udt_distance(self, z, R0):
        """Calculate UDT luminosity distance."""
        # UDT: d_L = z * R0 (pure temporal geometry)
        return z * R0
    
    def analyze_distance_contamination(self):
        """Analyze how LCDM distance assumptions affect supernova analysis."""
        print("\nDISTANCE CONTAMINATION ANALYSIS")
        print("-" * 33)
        
        # Get redshifts and magnitudes
        z = self.df['zCMB'].values
        mB = self.df['mB'].values
        
        # Calculate distances using both models
        d_L_lcdm = self.lcdm_distance(z)
        
        # Test different UDT R0 values
        R0_values = [2000, 3000, 4000, 5000, 6000]  # Mpc
        
        print("Distance comparison at different redshifts:")
        print("z\t\tLCDM d_L\tUDT d_L (R0=4000)")
        print("-" * 50)
        
        test_redshifts = [0.01, 0.03, 0.05, 0.07, 0.1]
        for z_test in test_redshifts:
            d_lcdm = self.lcdm_distance(z_test)
            d_udt = self.udt_distance(z_test, 4000)
            ratio = d_udt / d_lcdm
            print(f"{z_test:.3f}\t\t{d_lcdm:.1f} Mpc\t\t{d_udt:.1f} Mpc\t({ratio:.2f}x)")
        
        print()
        
        # Calculate what happens if we use wrong distances
        print("CONTAMINATION EFFECT ANALYSIS:")
        print("If supernova catalogs use LCDM distances but UDT is correct...")
        print()
        
        # Scenario: True UDT universe, but we measure using LCDM
        R0_true = 4000  # Mpc (true UDT parameter)
        M_B_true = -18.6  # Absolute magnitude
        
        # True UDT distances
        d_L_true_udt = self.udt_distance(z, R0_true)
        
        # But observers think they're at LCDM distances
        d_L_assumed_lcdm = self.lcdm_distance(z)
        
        # The true distance modulus (what we should measure)
        mu_true = 5 * np.log10(d_L_true_udt * 1e5)
        
        # But we calculate distance modulus assuming LCDM
        mu_assumed = 5 * np.log10(d_L_assumed_lcdm * 1e5)
        
        # The magnitude we actually observe
        mB_observed = M_B_true + mu_true
        
        # But we interpret it using LCDM distances
        M_B_inferred = mB_observed - mu_assumed
        
        print(f"True UDT distances: {np.mean(d_L_true_udt):.1f} ± {np.std(d_L_true_udt):.1f} Mpc")
        print(f"Assumed LCDM distances: {np.mean(d_L_assumed_lcdm):.1f} ± {np.std(d_L_assumed_lcdm):.1f} Mpc")
        print(f"Distance ratio (UDT/LCDM): {np.mean(d_L_true_udt/d_L_assumed_lcdm):.2f}")
        print()
        
        print(f"True absolute magnitude: {M_B_true:.2f}")
        print(f"Inferred absolute magnitude (using LCDM): {np.mean(M_B_inferred):.2f} ± {np.std(M_B_inferred):.2f}")
        print(f"Systematic error: {np.mean(M_B_inferred) - M_B_true:.2f} mag")
        
        return d_L_true_udt, d_L_assumed_lcdm, M_B_inferred
    
    def test_pure_redshift_analysis(self):
        """Test UDT using pure redshift-magnitude relation."""
        print("\nPURE REDSHIFT-MAGNITUDE ANALYSIS")
        print("-" * 34)
        
        print("Testing UDT directly from redshift and magnitude:")
        print("If UDT is correct: m = M + 5*log10(z*R0*1e5)")
        print("Rearranging: R0 = 10^((m-M)/5) / (z*1e5)")
        print()
        
        # Get data
        z = self.df['zCMB'].values
        mB = self.df['mB'].values
        
        # Assume standard absolute magnitude
        M_B_standard = -18.6
        
        # Calculate R0 for each supernova
        R0_individual = 10**((mB - M_B_standard)/5) / (z * 1e5)
        
        print(f"Individual R0 values:")
        print(f"  Mean: {np.mean(R0_individual):.1f} Mpc")
        print(f"  Median: {np.median(R0_individual):.1f} Mpc")
        print(f"  Std: {np.std(R0_individual):.1f} Mpc")
        print(f"  Range: {np.min(R0_individual):.1f} - {np.max(R0_individual):.1f} Mpc")
        print()
        
        # Calculate scatter
        R0_mean = np.mean(R0_individual)
        residuals = R0_individual - R0_mean
        rms_scatter = np.sqrt(np.mean(residuals**2))
        
        print(f"R0 scatter analysis:")
        print(f"  RMS scatter: {rms_scatter:.1f} Mpc")
        print(f"  Fractional scatter: {rms_scatter/R0_mean:.3f}")
        print()
        
        # Test if this is consistent with UDT
        if rms_scatter/R0_mean < 0.1:
            print("RESULT: LOW SCATTER - Consistent with UDT")
        elif rms_scatter/R0_mean < 0.3:
            print("RESULT: MODERATE SCATTER - Possibly consistent with UDT")
        else:
            print("RESULT: HIGH SCATTER - Not consistent with UDT")
        
        return R0_individual
    
    def plot_contamination_effects(self):
        """Plot the effects of distance contamination."""
        print("\nCreating contamination analysis plots...")
        
        # Test redshift range
        z_test = np.linspace(0.01, 0.1, 100)
        
        # Calculate distances
        d_L_lcdm = self.lcdm_distance(z_test)
        d_L_udt_2000 = self.udt_distance(z_test, 2000)
        d_L_udt_4000 = self.udt_distance(z_test, 4000)
        d_L_udt_6000 = self.udt_distance(z_test, 6000)
        
        # Create plots
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Panel 1: Distance vs redshift
        ax1 = axes[0, 0]
        ax1.plot(z_test, d_L_lcdm, 'k-', linewidth=2, label='LCDM')
        ax1.plot(z_test, d_L_udt_2000, 'r--', linewidth=2, label='UDT (R0=2000)')
        ax1.plot(z_test, d_L_udt_4000, 'g--', linewidth=2, label='UDT (R0=4000)')
        ax1.plot(z_test, d_L_udt_6000, 'b--', linewidth=2, label='UDT (R0=6000)')
        ax1.set_xlabel('Redshift z')
        ax1.set_ylabel('Luminosity Distance (Mpc)')
        ax1.set_title('Distance-Redshift Relations')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Panel 2: Distance ratio
        ax2 = axes[0, 1]
        ax2.plot(z_test, d_L_udt_2000/d_L_lcdm, 'r-', linewidth=2, label='UDT/LCDM (R0=2000)')
        ax2.plot(z_test, d_L_udt_4000/d_L_lcdm, 'g-', linewidth=2, label='UDT/LCDM (R0=4000)')
        ax2.plot(z_test, d_L_udt_6000/d_L_lcdm, 'b-', linewidth=2, label='UDT/LCDM (R0=6000)')
        ax2.axhline(y=1, color='k', linestyle=':', alpha=0.5)
        ax2.set_xlabel('Redshift z')
        ax2.set_ylabel('Distance Ratio')
        ax2.set_title('UDT/LCDM Distance Ratios')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Panel 3: Magnitude contamination
        ax3 = axes[1, 0]
        M_B_true = -18.6
        
        # True magnitudes (UDT universe)
        mB_true_udt = M_B_true + 5 * np.log10(d_L_udt_4000 * 1e5)
        
        # Magnitudes interpreted with LCDM
        M_B_inferred = mB_true_udt - 5 * np.log10(d_L_lcdm * 1e5)
        
        ax3.plot(z_test, M_B_inferred, 'r-', linewidth=2, label='Inferred M_B (using LCDM)')
        ax3.axhline(y=M_B_true, color='k', linestyle='--', label='True M_B')
        ax3.set_xlabel('Redshift z')
        ax3.set_ylabel('Absolute Magnitude')
        ax3.set_title('Magnitude Contamination Effect')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Panel 4: R0 scatter from data
        ax4 = axes[1, 1]
        if hasattr(self, 'df') and len(self.df) > 0:
            z_data = self.df['zCMB'].values
            mB_data = self.df['mB'].values
            
            # Calculate R0 for each supernova
            M_B_standard = -18.6
            R0_individual = 10**((mB_data - M_B_standard)/5) / (z_data * 1e5)
            
            ax4.scatter(z_data, R0_individual, alpha=0.6, s=20)
            ax4.axhline(y=np.mean(R0_individual), color='r', linestyle='--', 
                       label=f'Mean R0 = {np.mean(R0_individual):.0f} Mpc')
            ax4.set_xlabel('Redshift z')
            ax4.set_ylabel('R0 (Mpc)')
            ax4.set_title('R0 Scatter from Data')
            ax4.legend()
            ax4.grid(True, alpha=0.3)
        else:
            ax4.text(0.5, 0.5, 'No data available', transform=ax4.transAxes, 
                    ha='center', va='center')
        
        plt.tight_layout()
        plt.savefig('C:/UDT/results/supernova_contamination_analysis.png', dpi=150)
        plt.close()
        
        print("Contamination analysis plot saved to: C:/UDT/results/supernova_contamination_analysis.png")
    
    def run_complete_analysis(self):
        """Run complete contamination analysis."""
        print("COMPLETE CONTAMINATION ANALYSIS")
        print("=" * 35)
        
        # Analyze distance contamination
        d_L_true_udt, d_L_assumed_lcdm, M_B_inferred = self.analyze_distance_contamination()
        
        # Test pure redshift analysis
        R0_individual = self.test_pure_redshift_analysis()
        
        # Create plots
        self.plot_contamination_effects()
        
        print("\n" + "=" * 50)
        print("CONTAMINATION ANALYSIS CONCLUSIONS")
        print("=" * 50)
        
        print("\n1. DISTANCE CONTAMINATION EFFECTS:")
        if hasattr(self, 'df') and len(self.df) > 0:
            z = self.df['zCMB'].values
            d_L_lcdm = self.lcdm_distance(z)
            d_L_udt = self.udt_distance(z, 4000)
            
            mean_ratio = np.mean(d_L_udt / d_L_lcdm)
            print(f"   UDT/LCDM distance ratio: {mean_ratio:.2f}")
            
            if mean_ratio > 1.5:
                print("   SIGNIFICANT: UDT predicts much larger distances")
            elif mean_ratio < 0.7:
                print("   SIGNIFICANT: UDT predicts much smaller distances")
            else:
                print("   MODERATE: UDT and LCDM distances are comparable")
        
        print("\n2. MAGNITUDE CONTAMINATION:")
        if 'M_B_inferred' in locals():
            systematic_error = np.mean(M_B_inferred) - (-18.6)
            print(f"   Systematic error in absolute magnitude: {systematic_error:.2f} mag")
            
            if abs(systematic_error) > 0.1:
                print("   SIGNIFICANT: Large systematic error from distance assumptions")
            else:
                print("   MINOR: Small systematic error from distance assumptions")
        
        print("\n3. R0 CONSISTENCY:")
        if 'R0_individual' in locals():
            fractional_scatter = np.std(R0_individual) / np.mean(R0_individual)
            print(f"   R0 fractional scatter: {fractional_scatter:.3f}")
            
            if fractional_scatter < 0.1:
                print("   EXCELLENT: Low scatter supports UDT")
            elif fractional_scatter < 0.3:
                print("   GOOD: Moderate scatter, possibly consistent with UDT")
            else:
                print("   POOR: High scatter, not consistent with UDT")
        
        print("\n4. OVERALL ASSESSMENT:")
        print("   The poor supernova fits may be due to:")
        print("   a) LCDM distance contamination in the data")
        print("   b) Fundamental problems with UDT cosmology")
        print("   c) Missing evolutionary effects in UDT")
        print("   d) Different physics at cosmological scales")
        print()
        print("   RECOMMENDATION: Analyze supernova data using pure")
        print("   redshift-magnitude relations without distance assumptions.")

def main():
    """Run contamination analysis."""
    analysis = SupernovaDistanceContamination()
    analysis.run_complete_analysis()

if __name__ == "__main__":
    main()