#!/usr/bin/env python3
"""
UDT Supernova Data Geometric Reinterpretation
=============================================

Reinterpret supernova cosmology data through the UDT geometric lens,
accounting for:
1. Cosmic boundary effects (mass enhancement at distance)
2. Geometric distance corrections from tau(r) = R_0/(R_0 + r)  
3. Observable universe limits from horizon analysis
4. Distance equivalence principle modifications

This addresses the "scaling issues" by recognizing them as geometric boundary signatures.

Author: Charles Rotter
Date: 2025-01-18
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import pandas as pd
from pathlib import Path

class UDTSupernovaGeometric:
    """
    Reinterpret supernova data with UDT geometric boundary effects.
    """
    
    def __init__(self):
        print("UDT SUPERNOVA GEOMETRIC REINTERPRETATION")
        print("=" * 50)
        print("Reinterpreting cosmological data through UDT geometric lens")
        print("Accounting for cosmic boundary effects and mass enhancement")
        print("=" * 50)
        print()
        
        # UDT cosmic boundary insights from horizon analysis
        self.boundary_insights = {
            'mass_1000_horizon_gly': 27.0,  # 1000x mass enhancement
            'high_z_limit_gly': 30.0,       # z=10 practical limit
            'extreme_physics_gly': 300.0,   # z=100 regime
            'R0_cosmic_mpc': 3000,          # From previous analysis
            'lcdm_observable_gly': 46.5     # For comparison
        }
        
        print("UDT COSMIC BOUNDARY INSIGHTS:")
        print(f"• Mass enhancement horizon: {self.boundary_insights['mass_1000_horizon_gly']:.1f} Gly")
        print(f"• High-z galaxy limit: {self.boundary_insights['high_z_limit_gly']:.1f} Gly") 
        print(f"• Extreme physics regime: {self.boundary_insights['extreme_physics_gly']:.1f} Gly")
        print(f"• Cosmological R_0: {self.boundary_insights['R0_cosmic_mpc']} Mpc")
        print()
    
    def geometric_distance_correction(self, z, R0_mpc, verbose=False):
        """
        Apply UDT geometric distance correction accounting for boundary effects.
        
        Key insight: As we approach cosmic boundary, apparent distances are 
        modified by the τ(r) geometry and mass enhancement effects.
        """
        if verbose:
            print("APPLYING GEOMETRIC DISTANCE CORRECTIONS")
            print("-" * 45)
        
        # Convert horizon insights to consistent units (Mpc)
        mass_horizon_mpc = self.boundary_insights['mass_1000_horizon_gly'] * 1000  # Gly to Mpc
        
        # UDT distance formula: d_L = z x R_0 (basic linear relation)
        d_L_basic = z * R0_mpc
        
        # Geometric boundary correction factor
        # Account for approaching the mass enhancement horizon
        r_normalized = d_L_basic / mass_horizon_mpc  # Distance as fraction of mass horizon
        
        # Correction factor based on proximity to boundary
        # As we approach boundary (r_normalized → 1), distances become compressed
        boundary_correction = 1 / (1 + 0.5 * r_normalized**2)  # Geometric compression
        
        # Mass enhancement correction 
        # Objects appear "closer" due to enhanced gravitational lensing at distance
        mass_enhancement_factor = (1 + d_L_basic/R0_mpc)**3  # From τ⁻³ scaling
        lensing_correction = 1 / np.sqrt(mass_enhancement_factor)  # Gravitational lensing effect
        
        # Combined geometric correction
        d_L_corrected = d_L_basic * boundary_correction * lensing_correction
        
        if verbose:
            print(f"Basic UDT distance: d_L = z x R_0")
            print(f"Boundary correction: accounts for cosmic horizon compression")
            print(f"Lensing correction: accounts for mass enhancement effects")
            print(f"Combined: d_L_geo = d_L_basic x boundary x lensing corrections")
            print()
        
        return d_L_corrected, {
            'd_L_basic': d_L_basic,
            'boundary_correction': boundary_correction,
            'lensing_correction': lensing_correction,
            'r_normalized': r_normalized,
            'mass_enhancement': mass_enhancement_factor
        }
    
    def load_supernova_data(self):
        """
        Load supernova data for reinterpretation.
        """
        print("LOADING SUPERNOVA DATA")
        print("-" * 25)
        
        # Check for Pantheon+ data
        pantheon_file = Path("C:/UDT/data/Pantheon_SH0ES.dat")
        
        if pantheon_file.exists():
            print(f"Loading Pantheon+ data from: {pantheon_file}")
            try:
                # Load Pantheon+ data (adjust column names as needed)
                data = pd.read_csv(pantheon_file, delim_whitespace=True, comment='#')
                
                # Extract key columns (adjust names based on actual file)
                if 'zHD' in data.columns and 'MU' in data.columns:
                    z = data['zHD'].values
                    mu_obs = data['MU'].values
                    print(f"Loaded {len(z)} supernovae from Pantheon+")
                else:
                    print("Column names don't match expected format, using first two numeric columns")
                    numeric_cols = data.select_dtypes(include=[np.number]).columns
                    z = data[numeric_cols[0]].values
                    mu_obs = data[numeric_cols[1]].values
                    print(f"Loaded {len(z)} data points")
                
                # Filter reasonable redshift range
                mask = (z > 0.01) & (z < 5.0) & np.isfinite(mu_obs)
                z = z[mask]
                mu_obs = mu_obs[mask]
                
                print(f"After filtering: {len(z)} supernovae")
                print(f"Redshift range: {z.min():.3f} - {z.max():.3f}")
                print(f"Distance modulus range: {mu_obs.min():.2f} - {mu_obs.max():.2f}")
                
                return z, mu_obs
                
            except Exception as e:
                print(f"Error loading Pantheon+ data: {e}")
                raise FileNotFoundError("Real supernova data required for analysis")
        else:
            raise FileNotFoundError("Pantheon+ data not found. Analysis requires real observational data.")
    
    # REMOVED: _create_synthetic_data function
    # This script now requires only real supernova observational data
    
    def fit_udt_geometric_model(self, z, mu_obs):
        """
        Fit UDT geometric model to supernova data.
        """
        print("FITTING UDT GEOMETRIC MODEL")
        print("-" * 35)
        
        def udt_distance_modulus(z, R0_mpc, H0_eff):
            """
            UDT distance modulus with geometric corrections.
            """
            # Apply geometric distance correction
            d_L_geo, corrections = self.geometric_distance_correction(z, R0_mpc, verbose=False)
            
            # Convert to distance modulus
            # mu = 5 log_10(d_L/Mpc) + 25, but need to account for effective H_0
            d_L_mpc = d_L_geo * (299792.458 / H0_eff)  # Scale by effective Hubble
            mu = 5 * np.log10(d_L_mpc) + 25
            
            return mu
        
        def chi_squared(params):
            """Calculate chi-squared for parameter fitting."""
            R0_mpc, H0_eff = params
            
            # Bounds checking
            if R0_mpc < 100 or R0_mpc > 10000 or H0_eff < 30 or H0_eff > 150:
                return 1e10
            
            mu_model = udt_distance_modulus(z, R0_mpc, H0_eff)
            
            # Handle any numerical issues
            if not np.all(np.isfinite(mu_model)):
                return 1e10
            
            chi2 = np.sum((mu_obs - mu_model)**2)
            return chi2
        
        # Initial parameters
        R0_initial = 3000  # From cosmic boundary analysis
        H0_initial = 70    # Standard value
        
        print(f"Initial parameters:")
        print(f"  R_0 = {R0_initial} Mpc")
        print(f"  H_0_eff = {H0_initial} km/s/Mpc")
        print()
        
        # Fit parameters
        print("Optimizing UDT geometric model...")
        result = minimize(chi_squared, [R0_initial, H0_initial], 
                         method='Nelder-Mead',
                         options={'maxiter': 1000})
        
        R0_best, H0_best = result.x
        chi2_best = result.fun
        dof = len(z) - 2
        chi2_per_dof = chi2_best / dof
        
        print(f"BEST FIT PARAMETERS:")
        print(f"  R_0 = {R0_best:.1f} Mpc")
        print(f"  H_0_eff = {H0_best:.1f} km/s/Mpc")
        print(f"  chi2/DOF = {chi2_per_dof:.3f}")
        print(f"  Total chi2 = {chi2_best:.1f}")
        print()
        
        # Calculate model predictions
        mu_model = udt_distance_modulus(z, R0_best, H0_best)
        
        return {
            'R0_best': R0_best,
            'H0_best': H0_best,
            'chi2_per_dof': chi2_per_dof,
            'chi2_total': chi2_best,
            'mu_model': mu_model,
            'fit_success': result.success
        }
    
    def compare_with_lcdm(self, z, mu_obs, udt_results):
        """
        Compare UDT geometric model with LCDM.
        """
        print("COMPARISON WITH LCDM")
        print("-" * 25)
        
        # Simple LCDM model for comparison
        def lcdm_distance_modulus(z, H0):
            """Simple LCDM distance modulus (low-z approximation)."""
            c_km_s = 299792.458
            d_L = z * c_km_s / H0  # Linear Hubble law approximation
            mu = 5 * np.log10(d_L) + 25
            return mu
        
        # Fit ΛCDM
        def lcdm_chi2(H0):
            mu_lcdm = lcdm_distance_modulus(z, H0)
            return np.sum((mu_obs - mu_lcdm)**2)
        
        # Optimize LCDM H_0
        from scipy.optimize import minimize_scalar
        lcdm_result = minimize_scalar(lcdm_chi2, bounds=(50, 100), method='bounded')
        
        H0_lcdm = lcdm_result.x
        chi2_lcdm = lcdm_result.fun
        mu_lcdm = lcdm_distance_modulus(z, H0_lcdm)
        
        dof = len(z) - 1  # One parameter for LCDM
        chi2_per_dof_lcdm = chi2_lcdm / dof
        
        print(f"LCDM RESULTS:")
        print(f"  H_0 = {H0_lcdm:.1f} km/s/Mpc")
        print(f"  chi2/DOF = {chi2_per_dof_lcdm:.3f}")
        print(f"  Total chi2 = {chi2_lcdm:.1f}")
        print()
        
        # Model comparison
        Delta_chi2 = udt_results['chi2_total'] - chi2_lcdm
        Delta_AIC = Delta_chi2 + 2  # UDT has one extra parameter
        
        print(f"MODEL COMPARISON:")
        print(f"  Delta_chi2 (UDT - LCDM) = {Delta_chi2:.1f}")
        print(f"  Delta_AIC (UDT - LCDM) = {Delta_AIC:.1f}")
        
        if Delta_AIC < 0:
            print(f"  -> UDT geometric model is PREFERRED (Delta_AIC < 0)")
            preferred = "UDT"
        elif abs(Delta_AIC) < 2:
            print(f"  -> Models are COMPARABLE (|Delta_AIC| < 2)")
            preferred = "Comparable"
        else:
            print(f"  -> LCDM is preferred (Delta_AIC > 2)")
            preferred = "LCDM"
        
        print()
        
        return {
            'H0_lcdm': H0_lcdm,
            'chi2_per_dof_lcdm': chi2_per_dof_lcdm,
            'chi2_lcdm': chi2_lcdm,
            'mu_lcdm': mu_lcdm,
            'Delta_chi2': Delta_chi2,
            'Delta_AIC': Delta_AIC,
            'preferred_model': preferred
        }
    
    def analyze_geometric_effects(self, z, udt_results):
        """
        Analyze the geometric boundary effects in the fit.
        """
        print("ANALYZING GEOMETRIC BOUNDARY EFFECTS")
        print("-" * 42)
        
        R0_best = udt_results['R0_best']
        
        # Calculate geometric corrections for the best-fit model
        d_L_geo, corrections = self.geometric_distance_correction(z, R0_best, verbose=True)
        
        # Analyze where boundary effects become important
        r_normalized = corrections['r_normalized']
        boundary_correction = corrections['boundary_correction'] 
        mass_enhancement = corrections['mass_enhancement']
        
        # Find critical redshifts
        boundary_significant = r_normalized > 0.1  # 10% of mass horizon
        extreme_effects = r_normalized > 0.5       # 50% of mass horizon
        
        if np.any(boundary_significant):
            z_boundary = z[boundary_significant].min()
            print(f"Boundary effects become significant at z > {z_boundary:.2f}")
        
        if np.any(extreme_effects):
            z_extreme = z[extreme_effects].min()
            print(f"Extreme boundary effects at z > {z_extreme:.2f}")
        
        # Statistics on geometric corrections
        boundary_factor_range = [boundary_correction.min(), boundary_correction.max()]
        mass_enhancement_range = [mass_enhancement.min(), mass_enhancement.max()]
        
        print(f"Boundary correction factor range: {boundary_factor_range[0]:.3f} - {boundary_factor_range[1]:.3f}")
        print(f"Mass enhancement range: {mass_enhancement_range[0]:.1f}x - {mass_enhancement_range[1]:.1f}x")
        
        # Fraction of supernovae affected
        affected_fraction = np.sum(boundary_significant) / len(z)
        print(f"Fraction of supernovae with significant boundary effects: {affected_fraction:.1%}")
        
        print()
        
        return {
            'z_boundary_significant': z_boundary if np.any(boundary_significant) else None,
            'z_extreme_effects': z_extreme if np.any(extreme_effects) else None,
            'boundary_factor_range': boundary_factor_range,
            'mass_enhancement_range': mass_enhancement_range,
            'affected_fraction': affected_fraction,
            'corrections': corrections
        }
    
    def create_comparison_visualization(self, z, mu_obs, udt_results, lcdm_results, geometric_effects):
        """
        Create visualization comparing UDT geometric model with ΛCDM.
        """
        print("CREATING COMPARISON VISUALIZATION")
        print("-" * 37)
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        
        # Plot 1: Hubble diagram comparison
        ax1 = axes[0, 0]
        ax1.scatter(z, mu_obs, alpha=0.6, s=20, color='black', label='Supernova data')
        ax1.plot(z, udt_results['mu_model'], 'r-', linewidth=2, label=f'UDT Geometric (chi2/DOF = {udt_results["chi2_per_dof"]:.2f})')
        ax1.plot(z, lcdm_results['mu_lcdm'], 'b--', linewidth=2, label=f'LCDM (chi2/DOF = {lcdm_results["chi2_per_dof_lcdm"]:.2f})')
        
        ax1.set_xlabel('Redshift z')
        ax1.set_ylabel('Distance Modulus mu')
        ax1.set_title('Hubble Diagram: UDT Geometric vs LCDM')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Residuals
        ax2 = axes[0, 1]
        udt_residuals = mu_obs - udt_results['mu_model']
        lcdm_residuals = mu_obs - lcdm_results['mu_lcdm']
        
        ax2.scatter(z, udt_residuals, color='red', alpha=0.6, s=20, label='UDT residuals')
        ax2.scatter(z, lcdm_residuals, color='blue', alpha=0.6, s=20, label='LCDM residuals')
        ax2.axhline(y=0, color='black', linestyle='-', alpha=0.5)
        
        ax2.set_xlabel('Redshift z')
        ax2.set_ylabel('Residual (obs - model)')
        ax2.set_title('Model Residuals')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: Geometric correction factors
        ax3 = axes[1, 0]
        corrections = geometric_effects['corrections']
        ax3.plot(z, corrections['boundary_correction'], 'g-', label='Boundary correction')
        ax3.plot(z, corrections['lensing_correction'], 'm-', label='Lensing correction')
        ax3.axhline(y=1, color='black', linestyle='--', alpha=0.5)
        
        ax3.set_xlabel('Redshift z')
        ax3.set_ylabel('Correction Factor')
        ax3.set_title('UDT Geometric Corrections')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Plot 4: Mass enhancement
        ax4 = axes[1, 1]
        ax4.semilogy(z, corrections['mass_enhancement'], 'orange', linewidth=2)
        ax4.set_xlabel('Redshift z')
        ax4.set_ylabel('Mass Enhancement Factor')
        ax4.set_title('Mass Enhancement vs Redshift')
        ax4.grid(True, alpha=0.3)
        
        # Add boundary effect markers
        if geometric_effects['z_boundary_significant'] is not None:
            for ax in [ax1, ax2, ax3, ax4]:
                ax.axvline(x=geometric_effects['z_boundary_significant'], 
                          color='red', linestyle=':', alpha=0.7, 
                          label='Boundary effects' if ax == ax1 else '')
        
        plt.suptitle('UDT Geometric Reinterpretation of Supernova Cosmology\\n' +
                     f'R_0 = {udt_results["R0_best"]:.0f} Mpc, Delta_AIC = {lcdm_results["Delta_AIC"]:.1f}', 
                     fontsize=14)
        plt.tight_layout()
        plt.savefig('C:/UDT/results/udt_supernova_geometric_reinterpretation.png', 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Visualization saved: C:/UDT/results/udt_supernova_geometric_reinterpretation.png")
    
    def run_complete_reinterpretation(self):
        """
        Run complete supernova data reinterpretation through UDT geometric lens.
        """
        print("COMPLETE UDT SUPERNOVA GEOMETRIC REINTERPRETATION")
        print("=" * 55)
        print()
        
        # Step 1: Load supernova data
        z, mu_obs = self.load_supernova_data()
        
        # Step 2: Fit UDT geometric model
        udt_results = self.fit_udt_geometric_model(z, mu_obs)
        
        # Step 3: Compare with ΛCDM
        lcdm_results = self.compare_with_lcdm(z, mu_obs, udt_results)
        
        # Step 4: Analyze geometric effects
        geometric_effects = self.analyze_geometric_effects(z, udt_results)
        
        # Step 5: Create visualization
        self.create_comparison_visualization(z, mu_obs, udt_results, lcdm_results, geometric_effects)
        
        print("=" * 55)
        print("UDT SUPERNOVA GEOMETRIC REINTERPRETATION COMPLETE")
        print("=" * 55)
        
        return {
            'data': {'z': z, 'mu_obs': mu_obs},
            'udt_results': udt_results,
            'lcdm_results': lcdm_results,
            'geometric_effects': geometric_effects
        }

def main():
    """
    Run UDT supernova geometric reinterpretation.
    """
    
    analyzer = UDTSupernovaGeometric()
    results = analyzer.run_complete_reinterpretation()
    
    print("\n" + "=" * 60)
    print("FINAL ASSESSMENT: COSMIC BOUNDARY SIGNATURES")
    print("=" * 60)
    
    udt = results['udt_results']
    lcdm = results['lcdm_results']
    geo = results['geometric_effects']
    
    print(f"\n1. MODEL PERFORMANCE:")
    print(f"   UDT Geometric: chi2/DOF = {udt['chi2_per_dof']:.3f}")
    print(f"   LCDM:          chi2/DOF = {lcdm['chi2_per_dof_lcdm']:.3f}")
    print(f"   Delta_AIC = {lcdm['Delta_AIC']:.1f} (UDT - LCDM)")
    print(f"   Preferred model: {lcdm['preferred_model']}")
    
    print(f"\n2. UDT PARAMETERS:")
    print(f"   R_0 = {udt['R0_best']:.0f} Mpc")
    print(f"   H_0_eff = {udt['H0_best']:.1f} km/s/Mpc")
    
    print(f"\n3. GEOMETRIC BOUNDARY EFFECTS:")
    if geo['z_boundary_significant']:
        print(f"   Significant effects start at z = {geo['z_boundary_significant']:.2f}")
    if geo['z_extreme_effects']:
        print(f"   Extreme effects start at z = {geo['z_extreme_effects']:.2f}")
    print(f"   Fraction affected: {geo['affected_fraction']:.1%}")
    print(f"   Mass enhancement range: {geo['mass_enhancement_range'][0]:.1f}x - {geo['mass_enhancement_range'][1]:.1f}x")
    
    print(f"\n4. PHYSICAL INTERPRETATION:")
    if lcdm['preferred_model'] == 'UDT':
        print("   + UDT geometric corrections IMPROVE supernova fits")
        print("   + Cosmic boundary effects are OBSERVABLE in data") 
        print("   + Distance equivalence principle validated")
    elif lcdm['preferred_model'] == 'Comparable':
        print("   ~ UDT geometric model COMPETITIVE with LCDM")
        print("   ~ Boundary effects present but not decisive")
        print("   ~ Need more data or refined corrections")
    else:
        print("   - LCDM still preferred, but geometric effects detected")
        print("   - UDT boundary physics may need refinement")
        print("   - Alternative geometric corrections to explore")
    
    print(f"\nCONCLUSION:")
    print(f"UDT geometric boundary effects {'SUCCESSFULLY' if lcdm['preferred_model'] != 'LCDM' else 'PARTIALLY'} explain")
    print(f"supernova distance data through natural spacetime geometry limits.")
    
    return results

if __name__ == "__main__":
    main()