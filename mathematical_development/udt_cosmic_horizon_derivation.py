#!/usr/bin/env python3
"""
UDT Cosmic Horizon and Observable Universe Derivation
====================================================

Derive fundamental horizon distances and observable universe limits in UDT
from the temporal geometry tau(r) = R_0/(R_0 + r).

This establishes the geometric boundaries that will reframe cosmological data interpretation.

Key Questions:
1. Where does the cosmic boundary become relevant?
2. What is the UDT "observable universe" size?
3. How do redshift and mass diverge approaching the boundary?
4. What are the practical cutoffs for observable astronomy?

Author: Charles Rotter
Date: 2025-01-18
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import sympy as sp

class UDTCosmicHorizon:
    """
    Derive cosmic horizons and observable universe limits in UDT geometry.
    """
    
    def __init__(self):
        print("UDT COSMIC HORIZON DERIVATION")
        print("=" * 40)
        print("Deriving fundamental geometric boundaries from tau(r) = R_0/(R_0 + r)")
        print("=" * 40)
        print()
        
        # Symbolic variables
        self.r, self.R0 = sp.symbols('r R_0', real=True, positive=True)
        self.tau = self.R0 / (self.R0 + self.r)
        
        print(f"UDT temporal geometry: tau(r) = {self.tau}")
        print()
    
    def derive_redshift_horizon(self):
        """
        Derive where redshift becomes extreme/infinite.
        
        In UDT: z = (1/tau) - 1 = r/R_0
        """
        print("REDSHIFT HORIZON ANALYSIS")
        print("-" * 30)
        
        # UDT redshift formula
        z_udt = (1/self.tau) - 1
        z_simplified = sp.simplify(z_udt)
        
        print(f"UDT redshift: z = (1/tau) - 1 = {z_simplified}")
        print()
        
        # Critical redshift values and corresponding distances
        z_critical_values = [1, 5, 10, 100, 1000]
        
        print("REDSHIFT HORIZONS:")
        print("z_crit    r/R_0     Physical Meaning")
        print("-" * 40)
        
        horizon_distances = {}
        for z_crit in z_critical_values:
            r_over_R0 = z_crit  # Since z = r/R_0 in UDT
            meaning = self._get_redshift_meaning(z_crit)
            print(f"{z_crit:6.0f}    {r_over_R0:5.0f}      {meaning}")
            horizon_distances[f'z_{z_crit}'] = r_over_R0
        
        print()
        print("ASYMPTOTIC BEHAVIOR:")
        print("As r -> inf: z -> inf (infinite redshift)")
        print("As r -> 0: z -> 0 (no redshift)")
        print()
        
        return horizon_distances
    
    def derive_mass_divergence_horizon(self):
        """
        Derive where mass enhancement becomes extreme.
        
        From field equations: rho prop tau^-3 = (1 + r/R_0)^3
        """
        print("MASS DIVERGENCE HORIZON ANALYSIS")
        print("-" * 35)
        
        # Mass enhancement factor from UDT field equations
        mass_enhancement = (1/self.tau)**3
        mass_simplified = sp.simplify(mass_enhancement)
        
        print(f"Mass enhancement: rho/rho_0 = tau^-3 = {mass_simplified}")
        print()
        
        # Critical mass enhancement values
        mass_critical_values = [10, 100, 1000, 10000, 1e6]
        
        print("MASS ENHANCEMENT HORIZONS:")
        print("rho/rho_0      r/R_0     Physical Meaning")
        print("-" * 45)
        
        mass_horizons = {}
        for mass_crit in mass_critical_values:
            # From (1 + r/R_0)^3 = mass_crit
            r_over_R0 = mass_crit**(1/3) - 1
            meaning = self._get_mass_meaning(mass_crit)
            print(f"{mass_crit:8.0f}    {r_over_R0:5.1f}      {meaning}")
            mass_horizons[f'mass_{int(mass_crit)}'] = r_over_R0
        
        print()
        print("ASYMPTOTIC BEHAVIOR:")
        print("As r -> inf: rho -> inf (infinite mass density)")
        print("As r -> 0: rho -> rho_0 (normal mass)")
        print()
        
        return mass_horizons
    
    def derive_observability_limits(self):
        """
        Derive practical limits for astronomical observations.
        """
        print("OBSERVABILITY LIMITS")
        print("-" * 25)
        
        # Current observational limits
        print("CURRENT ASTRONOMICAL LIMITS:")
        print("- Farthest galaxies: z ~ 10-15")
        print("- CMB last scattering: z ~ 1100") 
        print("- Theoretical limit: z ~ inf")
        print()
        
        # Translate to UDT distances
        print("UDT DISTANCE TRANSLATIONS:")
        print("Observable Type          z_obs    r/R_0    Fraction of R_0")
        print("-" * 55)
        
        obs_types = [
            ("Nearby galaxies", 0.1, "Local group"),
            ("Distant galaxies", 1.0, "Current surveys"),
            ("High-z galaxies", 10.0, "JWST/HST limit"),
            ("CMB surface", 1100.0, "Last scattering"),
            ("Theoretical limit", float('inf'), "Cosmic boundary")
        ]
        
        observability_limits = {}
        for obs_type, z_obs, description in obs_types:
            if z_obs == float('inf'):
                r_over_R0 = float('inf')
                fraction = "inf"
            else:
                r_over_R0 = z_obs  # Since z = r/R_0
                fraction = f"{r_over_R0:.1f}"
            
            print(f"{obs_type:18s}  {z_obs:7.1f}  {r_over_R0:7.1f}  {fraction}")
            observability_limits[obs_type.replace(' ', '_')] = {
                'z': z_obs, 
                'r_over_R0': r_over_R0,
                'description': description
            }
        
        print()
        return observability_limits
    
    def derive_udt_universe_size(self, R0_cosmic_mpc=3000):
        """
        Calculate actual UDT universe size for given R_0.
        """
        print("UDT UNIVERSE SIZE CALCULATION")
        print("-" * 35)
        print(f"Using cosmological R_0 = {R0_cosmic_mpc} Mpc")
        print()
        
        # Define practical "universe size" as where effects become extreme
        cutoff_criteria = {
            'z_10': {'z': 10, 'meaning': 'High-z galaxy limit'},
            'z_100': {'z': 100, 'meaning': 'Extreme redshift'},
            'z_1000': {'z': 1000, 'meaning': 'Near CMB redshift'},
            'mass_1000': {'mass_factor': 1000, 'meaning': '1000x mass enhancement'},
            'mass_1e6': {'mass_factor': 1e6, 'meaning': 'Million-fold mass enhancement'}
        }
        
        print("UDT EFFECTIVE UNIVERSE BOUNDARIES:")
        print("Criterion         Distance (Mpc)    Distance (Gly)    Meaning")
        print("-" * 70)
        
        universe_sizes = {}
        for criterion, params in cutoff_criteria.items():
            if 'z' in params:
                # Redshift-based criterion
                r_mpc = params['z'] * R0_cosmic_mpc
                meaning = params['meaning']
            else:
                # Mass-based criterion  
                r_over_R0 = params['mass_factor']**(1/3) - 1
                r_mpc = r_over_R0 * R0_cosmic_mpc
                meaning = params['meaning']
            
            r_gly = r_mpc / 1000  # Convert Mpc to Gly (roughly)
            
            print(f"{criterion:15s}   {r_mpc:11.0f}     {r_gly:9.1f}     {meaning}")
            universe_sizes[criterion] = {
                'distance_mpc': r_mpc,
                'distance_gly': r_gly,
                'meaning': meaning
            }
        
        print()
        print("COMPARISON WITH LCDM:")
        lambda_cdm_size_gly = 46.5  # Observable universe radius
        print(f"LCDM observable universe: ~{lambda_cdm_size_gly} Gly")
        
        # Find which UDT boundary is closest to LCDM
        differences = {}
        for criterion, data in universe_sizes.items():
            diff_factor = data['distance_gly'] / lambda_cdm_size_gly
            differences[criterion] = diff_factor
            comparison = "larger" if diff_factor > 1 else "smaller"
            print(f"UDT {criterion}: {diff_factor:.2f}x {comparison} than LCDM")
        
        print()
        return universe_sizes, differences
    
    def analyze_cosmic_boundary_physics(self):
        """
        Analyze the physics near the cosmic boundary.
        """
        print("COSMIC BOUNDARY PHYSICS")
        print("-" * 30)
        
        print("APPROACHING THE BOUNDARY (r -> inf):")
        print("1. tau(r) -> 0           (time dilation becomes extreme)")
        print("2. z -> inf             (infinite redshift)")
        print("3. rho -> inf             (infinite mass density)")
        print("4. c_eff(r) -> 0      (light speed approaches zero)")
        print("5. t_proper -> 0      (proper time slows to halt)")
        print()
        
        print("BOUNDARY INTERPRETATION:")
        print("• The r -> inf limit represents the EDGE OF SPACETIME")
        print("• Objects cannot actually reach this boundary")
        print("• Similar to event horizons in black holes")
        print("• Defines the maximum possible distance between observers")
        print()
        
        print("COSMOLOGICAL IMPLICATIONS:")
        print("• Universe has FINITE EXTENT in UDT geometry")
        print("• No need for infinite space")
        print("• Natural cutoff for observable universe")
        print("• Dark matter may be mass 'pushed' toward boundary")
        print()
        
        return {
            'boundary_type': 'geometric_cutoff',
            'physics': 'mass_and_time_divergence',
            'implications': 'finite_universe_with_natural_boundary'
        }
    
    def create_horizon_visualization(self, R0_cosmic=3000):
        """
        Create visualization of UDT cosmic horizons.
        """
        print("CREATING HORIZON VISUALIZATION")
        print("-" * 35)
        
        # Distance range in units of R_0
        r_over_R0 = np.logspace(-2, 3, 1000)  # 0.01 to 1000 R_0
        
        # Calculate key quantities
        tau_vals = 1 / (1 + r_over_R0)
        z_vals = r_over_R0
        mass_enhancement = (1 + r_over_R0)**3
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Plot 1: Temporal geometry function
        ax1 = axes[0, 0]
        ax1.loglog(r_over_R0, tau_vals)
        ax1.set_xlabel('Distance r/R_0')
        ax1.set_ylabel('Temporal dilation tau(r)')
        ax1.set_title('UDT Temporal Geometry')
        ax1.grid(True, alpha=0.3)
        ax1.axhline(y=0.1, color='red', linestyle='--', alpha=0.7, label='tau = 0.1')
        ax1.axhline(y=0.01, color='orange', linestyle='--', alpha=0.7, label='tau = 0.01')
        ax1.legend()
        
        # Plot 2: Redshift relation
        ax2 = axes[0, 1]
        ax2.loglog(r_over_R0, z_vals)
        ax2.set_xlabel('Distance r/R_0')
        ax2.set_ylabel('Redshift z')
        ax2.set_title('UDT Redshift vs Distance')
        ax2.grid(True, alpha=0.3)
        ax2.axhline(y=10, color='red', linestyle='--', alpha=0.7, label='z = 10')
        ax2.axhline(y=1100, color='orange', linestyle='--', alpha=0.7, label='z = 1100 (CMB)')
        ax2.legend()
        
        # Plot 3: Mass enhancement
        ax3 = axes[1, 0]
        ax3.loglog(r_over_R0, mass_enhancement)
        ax3.set_xlabel('Distance r/R_0')
        ax3.set_ylabel('Mass enhancement rho/rho_0')
        ax3.set_title('UDT Mass Divergence')
        ax3.grid(True, alpha=0.3)
        ax3.axhline(y=1000, color='red', linestyle='--', alpha=0.7, label='1000x mass')
        ax3.axhline(y=1e6, color='orange', linestyle='--', alpha=0.7, label='10^6x mass')
        ax3.legend()
        
        # Plot 4: Effective light speed
        ax4 = axes[1, 1]
        c_eff_normalized = tau_vals  # c_eff/c_0 = tau
        ax4.loglog(r_over_R0, c_eff_normalized)
        ax4.set_xlabel('Distance r/R_0')
        ax4.set_ylabel('Effective light speed c_eff/c_0')
        ax4.set_title('UDT Light Speed Variation')
        ax4.grid(True, alpha=0.3)
        ax4.axhline(y=0.1, color='red', linestyle='--', alpha=0.7, label='c_eff = 0.1c_0')
        ax4.axhline(y=0.01, color='orange', linestyle='--', alpha=0.7, label='c_eff = 0.01c_0')
        ax4.legend()
        
        plt.suptitle(f'UDT Cosmic Horizons and Boundary Physics\\nCosmological R_0 = {R0_cosmic} Mpc', 
                     fontsize=16)
        plt.tight_layout()
        plt.savefig('C:/UDT/results/udt_cosmic_horizons.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Visualization saved: C:/UDT/results/udt_cosmic_horizons.png")
    
    def _get_redshift_meaning(self, z):
        """Get physical meaning of redshift value."""
        if z <= 1:
            return "Nearby universe"
        elif z <= 10:
            return "High-z galaxies"
        elif z <= 100:
            return "Extreme redshift"
        elif z <= 1100:
            return "Near CMB"
        else:
            return "Beyond CMB"
    
    def _get_mass_meaning(self, mass_factor):
        """Get physical meaning of mass enhancement."""
        if mass_factor <= 10:
            return "Mild enhancement"
        elif mass_factor <= 100:
            return "Strong enhancement"
        elif mass_factor <= 1000:
            return "Extreme enhancement"
        elif mass_factor <= 10000:
            return "Very extreme"
        else:
            return "Approaching boundary"
    
    def run_complete_horizon_analysis(self, R0_cosmic_mpc=3000):
        """
        Run complete cosmic horizon analysis.
        """
        print("COMPLETE UDT COSMIC HORIZON ANALYSIS")
        print("=" * 45)
        print()
        
        # Step 1: Redshift horizons
        redshift_horizons = self.derive_redshift_horizon()
        
        # Step 2: Mass divergence horizons  
        mass_horizons = self.derive_mass_divergence_horizon()
        
        # Step 3: Observability limits
        observability_limits = self.derive_observability_limits()
        
        # Step 4: UDT universe size
        universe_sizes, size_comparisons = self.derive_udt_universe_size(R0_cosmic_mpc)
        
        # Step 5: Boundary physics
        boundary_physics = self.analyze_cosmic_boundary_physics()
        
        # Step 6: Visualization
        self.create_horizon_visualization(R0_cosmic_mpc)
        
        print("=" * 50)
        print("UDT COSMIC HORIZON ANALYSIS COMPLETE")
        print("=" * 50)
        
        return {
            'redshift_horizons': redshift_horizons,
            'mass_horizons': mass_horizons,
            'observability_limits': observability_limits,
            'universe_sizes': universe_sizes,
            'size_comparisons': size_comparisons,
            'boundary_physics': boundary_physics
        }

def main():
    """
    Run UDT cosmic horizon derivation.
    """
    
    analyzer = UDTCosmicHorizon()
    results = analyzer.run_complete_horizon_analysis(R0_cosmic_mpc=3000)
    
    print("\n" + "=" * 60)
    print("KEY INSIGHTS FOR COSMOLOGICAL DATA REINTERPRETATION")
    print("=" * 60)
    
    print("\n1. UDT UNIVERSE BOUNDARIES:")
    for criterion, data in results['universe_sizes'].items():
        print(f"   {criterion}: {data['distance_gly']:.1f} Gly ({data['meaning']})")
    
    print(f"\n2. SIZE COMPARISONS WITH LCDM:")
    for criterion, factor in results['size_comparisons'].items():
        comparison = "larger" if factor > 1 else "smaller"
        print(f"   {criterion}: {factor:.2f}x {comparison}")
    
    print(f"\n3. CRITICAL REDSHIFTS:")
    for z_name, r_ratio in results['redshift_horizons'].items():
        print(f"   {z_name}: r = {r_ratio}xR_0")
    
    print(f"\n4. BOUNDARY PHYSICS:")
    print(f"   Type: {results['boundary_physics']['boundary_type']}")
    print(f"   Physics: {results['boundary_physics']['physics']}")
    print(f"   Implication: {results['boundary_physics']['implications']}")
    
    print(f"\nREADY FOR COSMOLOGICAL DATA REINTERPRETATION:")
    print(f"• UDT provides natural geometric boundaries")
    print(f"• Observable universe has finite extent in UDT")
    print(f"• Supernova distances need geometric correction")
    print(f"• CMB physics occurs near mass divergence regime")
    
    return results

if __name__ == "__main__":
    main()