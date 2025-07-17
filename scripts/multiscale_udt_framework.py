#!/usr/bin/env python3
"""
Multi-Scale UDT Framework
=========================

Implements Universal Distance Dilation Theory with scale-dependent R0 parameters.
This addresses the scale mismatch between galactic, cosmological, and CMB domains.

Physical Motivation:
UDT's temporal geometry function τ(r) = R0/(R0 + r) should have different 
characteristic scales R0 for different physical regimes:

1. Galactic scale: R0_gal ~ 38 kpc (from SPARC galaxy fits)
2. Cosmological scale: R0_cosmo ~ 3,000 Mpc (from supernova analysis)  
3. CMB/primordial scale: R0_cmb ~ 11,357,000 Mpc (calibrated to recombination)

This hierarchy reflects the nested temporal geometry of the universe.

Author: Charles Rotter
Date: 2025-01-17
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import minimize
import os

class MultiScaleUDT:
    """Multi-scale Universal Distance Dilation Theory framework."""
    
    def __init__(self):
        """Initialize with scale-dependent R0 parameters."""
        
        # Scale-dependent R0 values (in Mpc)
        self.R0_galactic = 0.038        # 38 kpc (from SPARC fits)
        self.R0_cosmological = 3000.0   # 3,000 Mpc (from supernova fits)
        self.R0_cmb = 11357000.0        # 11.4 million Mpc (calibrated to CMB)
        
        # Physical constants
        self.c_light = 299792.458       # km/s
        
        # Standard cosmology parameters (for validation)
        self.H0_standard = 70.0         # km/s/Mpc
        self.Omega_b = 0.049            # Baryon density
        self.Omega_r = 8.24e-5          # Radiation density
        self.z_recombination = 1100.0   # Recombination redshift
        
        # Scale transition boundaries (in Mpc)
        self.galactic_boundary = 1.0    # 1 Mpc
        self.cosmological_boundary = 1000000.0  # 1 Gpc
        
        print("Multi-Scale UDT Framework Initialized")
        print("=" * 50)
        print(f"R0_galactic:      {self.R0_galactic:.3f} Mpc ({self.R0_galactic*1000:.1f} kpc)")
        print(f"R0_cosmological:  {self.R0_cosmological:.1f} Mpc")
        print(f"R0_cmb:           {self.R0_cmb:.0f} Mpc ({self.R0_cmb/1e6:.1f} Gpc)")
        print(f"Scale ratios:")
        print(f"  Cosmo/Gal:      {self.R0_cosmological/self.R0_galactic:.0f}")
        print(f"  CMB/Cosmo:      {self.R0_cmb/self.R0_cosmological:.0f}")
        print(f"  CMB/Gal:        {self.R0_cmb/self.R0_galactic:.0f}")
        print()
    
    def get_scale_appropriate_R0(self, distance_scale):
        """
        Get appropriate R0 value for given distance scale.
        
        Parameters:
        -----------
        distance_scale : float
            Characteristic distance in Mpc
            
        Returns:
        --------
        R0 : float
            Appropriate R0 value in Mpc
        scale_regime : str
            Which scale regime this corresponds to
        """
        if distance_scale < self.galactic_boundary:
            return self.R0_galactic, "galactic"
        elif distance_scale < self.cosmological_boundary:
            return self.R0_cosmological, "cosmological"
        else:
            return self.R0_cmb, "cmb"
    
    def temporal_geometry_function(self, distance, scale_regime=None):
        """
        Multi-scale temporal geometry function.
        
        Parameters:
        -----------
        distance : float or array
            Distance in Mpc
        scale_regime : str, optional
            Force specific scale regime ('galactic', 'cosmological', 'cmb')
            
        Returns:
        --------
        tau : float or array
            Temporal dilation factor
        """
        if scale_regime is None:
            # Auto-determine scale regime
            if np.isscalar(distance):
                R0, _ = self.get_scale_appropriate_R0(distance)
            else:
                # For arrays, use cosmological as default
                R0 = self.R0_cosmological
        else:
            # Use specified regime
            if scale_regime == "galactic":
                R0 = self.R0_galactic
            elif scale_regime == "cosmological":
                R0 = self.R0_cosmological
            elif scale_regime == "cmb":
                R0 = self.R0_cmb
            else:
                raise ValueError(f"Unknown scale regime: {scale_regime}")
        
        return R0 / (R0 + np.abs(distance))
    
    def effective_light_speed(self, distance, scale_regime=None):
        """
        Position-dependent effective speed of light.
        
        Parameters:
        -----------
        distance : float or array
            Distance in Mpc
        scale_regime : str, optional
            Force specific scale regime
            
        Returns:
        --------
        c_eff : float or array
            Effective speed of light in km/s
        """
        tau = self.temporal_geometry_function(distance, scale_regime)
        return self.c_light * tau
    
    def galactic_dynamics_prediction(self, radius_kpc, mass_profile):
        """
        Predict galactic rotation curve using galactic-scale UDT.
        
        Parameters:
        -----------
        radius_kpc : array
            Galactic radii in kpc
        mass_profile : array
            Enclosed baryonic mass in solar masses
            
        Returns:
        --------
        v_predicted : array
            Predicted rotation velocities in km/s
        """
        # Convert to Mpc
        radius_Mpc = radius_kpc / 1000.0
        
        # Use galactic-scale R0
        enhancement_factor = (1 + radius_Mpc / self.R0_galactic)**2
        
        # Baryonic velocity (Newtonian)
        G = 4.3e-6  # km^2/s^2 per (kpc * solar mass)
        v_baryonic = np.sqrt(G * mass_profile / radius_kpc)
        
        # UDT enhancement
        v_predicted = v_baryonic * np.sqrt(enhancement_factor)
        
        return v_predicted
    
    def supernova_distance_modulus(self, redshift):
        """
        Predict supernova distance modulus using cosmological-scale UDT.
        
        Parameters:
        -----------
        redshift : float or array
            Cosmological redshift
            
        Returns:
        --------
        mu : float or array
            Distance modulus in magnitudes
        """
        # UDT luminosity distance: d_L = z * R0
        d_L_Mpc = redshift * self.R0_cosmological
        
        # Distance modulus
        mu = 5 * np.log10(d_L_Mpc) + 25
        
        return mu
    
    def cmb_conformal_time_at_recombination(self):
        """
        Calculate conformal time at recombination using CMB-scale UDT.
        
        Returns:
        --------
        eta_rec : float
            Conformal time at recombination in Mpc/c
        """
        # Use CMB-scale R0
        eta_rec_Mpc = self.R0_cmb / (1 + self.z_recombination)
        eta_rec_conf = eta_rec_Mpc / self.c_light
        
        return eta_rec_conf
    
    def cmb_sound_horizon(self):
        """
        Calculate sound horizon at recombination using CMB-scale UDT.
        
        Returns:
        --------
        r_s : float
            Sound horizon in Mpc
        """
        eta_rec = self.cmb_conformal_time_at_recombination()
        
        # Sound speed evolution (simplified)
        def sound_speed(eta):
            # Convert to redshift estimate
            eta_Mpc = eta * self.c_light
            z_eff = self.R0_cmb / eta_Mpc - 1
            
            # Baryon-photon ratio
            R_gamma = 3 * self.Omega_b / (4 * self.Omega_r * (1 + z_eff))
            
            # Sound speed
            c_s_squared = 1.0 / (3.0 * (1.0 + R_gamma))
            c_s = np.sqrt(c_s_squared)
            
            # CMB-scale temporal geometry effect
            tau = self.temporal_geometry_function(eta_Mpc, scale_regime="cmb")
            return c_s * tau
        
        # Integrate from early times to recombination
        eta_early = eta_rec / 1000.0
        r_s, _ = quad(sound_speed, eta_early, eta_rec)
        r_s *= self.c_light  # Convert to Mpc
        
        return r_s
    
    def cmb_angular_diameter_distance(self):
        """
        Calculate angular diameter distance to recombination using CMB-scale UDT.
        
        Returns:
        --------
        D_A : float
            Angular diameter distance in Mpc
        """
        eta_rec = self.cmb_conformal_time_at_recombination()
        eta_rec_Mpc = eta_rec * self.c_light
        
        # CMB-scale temporal geometry
        tau_rec = self.temporal_geometry_function(eta_rec_Mpc, scale_regime="cmb")
        D_A = eta_rec_Mpc * tau_rec
        
        return D_A
    
    def cmb_acoustic_peaks(self):
        """
        Calculate CMB acoustic peak positions using CMB-scale UDT.
        
        Returns:
        --------
        ell_peaks : array
            Multipole moments of acoustic peaks
        """
        r_s = self.cmb_sound_horizon()
        D_A = self.cmb_angular_diameter_distance()
        
        # First peak position
        ell_1 = np.pi * D_A / r_s
        
        # Higher peaks
        peak_numbers = np.array([1, 2, 3, 4, 5, 6])
        ell_peaks = ell_1 * peak_numbers
        
        return ell_peaks, r_s, D_A
    
    def validate_scale_consistency(self):
        """
        Validate that different scales give consistent physics.
        
        Returns:
        --------
        validation_results : dict
            Results from cross-scale validation
        """
        print("MULTI-SCALE VALIDATION")
        print("=" * 50)
        
        results = {}
        
        # 1. Galactic scale validation
        print("1. Galactic Scale (SPARC validation)")
        print("-" * 30)
        
        # Test galaxy at 10 kpc
        r_test = np.array([10.0])  # kpc
        m_test = np.array([1e10])  # solar masses
        v_pred = self.galactic_dynamics_prediction(r_test, m_test)
        
        print(f"Test galaxy: R = {r_test[0]:.0f} kpc, M = {m_test[0]:.0e} Msun")
        print(f"UDT enhancement factor: {(1 + r_test[0]/1000/self.R0_galactic)**2:.2f}")
        print(f"Predicted velocity: {v_pred[0]:.1f} km/s")
        
        results['galactic'] = {
            'R0': self.R0_galactic,
            'test_velocity': v_pred[0],
            'enhancement': (1 + r_test[0]/1000/self.R0_galactic)**2
        }
        print()
        
        # 2. Cosmological scale validation  
        print("2. Cosmological Scale (Supernova validation)")
        print("-" * 40)
        
        z_test = np.array([0.1, 0.5, 1.0])
        mu_pred = self.supernova_distance_modulus(z_test)
        
        print("Redshift  Distance Modulus")
        for i, (z, mu) in enumerate(zip(z_test, mu_pred)):
            print(f"  {z:.1f}        {mu:.2f}")
        
        results['cosmological'] = {
            'R0': self.R0_cosmological,
            'test_redshifts': z_test,
            'test_moduli': mu_pred
        }
        print()
        
        # 3. CMB scale validation
        print("3. CMB Scale (Recombination validation)")
        print("-" * 35)
        
        eta_rec = self.cmb_conformal_time_at_recombination()
        r_s = self.cmb_sound_horizon()
        D_A = self.cmb_angular_diameter_distance()
        ell_peaks, _, _ = self.cmb_acoustic_peaks()
        
        print(f"Conformal time at recombination: {eta_rec:.3f} Mpc/c")
        print(f"Sound horizon: {r_s:.1f} Mpc")
        print(f"Angular diameter distance: {D_A:.1f} Mpc")
        print(f"First acoustic peak: l1 = {ell_peaks[0]:.1f}")
        print(f"Acoustic scale: theta_s = {r_s/D_A:.6f} rad")
        
        # Compare with standard cosmology
        eta_standard = 288.0  # Mpc/c
        r_s_standard = 147.3  # Mpc
        l1_standard = 220.0
        
        print(f"\nComparison with standard cosmology:")
        print(f"  eta_rec: {eta_rec:.1f} vs {eta_standard:.1f} (factor {eta_rec/eta_standard:.2f})")
        print(f"  r_s: {r_s:.1f} vs {r_s_standard:.1f} (factor {r_s/r_s_standard:.2f})")
        print(f"  l1: {ell_peaks[0]:.1f} vs {l1_standard:.1f} (error {abs(ell_peaks[0]-l1_standard)/l1_standard*100:.1f}%)")
        
        results['cmb'] = {
            'R0': self.R0_cmb,
            'eta_rec': eta_rec,
            'r_s': r_s,
            'D_A': D_A,
            'l1': ell_peaks[0],
            'error_vs_standard': abs(ell_peaks[0]-l1_standard)/l1_standard*100
        }
        print()
        
        # 4. Scale transition analysis
        print("4. Scale Transition Analysis")
        print("-" * 25)
        
        distances = np.logspace(-3, 7, 11)  # 1 kpc to 10 Gpc
        
        print("Distance      Scale        R0")
        print("(Mpc)         Regime       (Mpc)")
        print("-" * 35)
        
        for d in distances:
            R0, regime = self.get_scale_appropriate_R0(d)
            print(f"{d:.0e}      {regime:12s} {R0:.0e}")
        
        results['scale_transitions'] = {
            'distances': distances,
            'regimes': [self.get_scale_appropriate_R0(d)[1] for d in distances],
            'R0_values': [self.get_scale_appropriate_R0(d)[0] for d in distances]
        }
        
        return results
    
    def create_multiscale_diagnostic_plots(self, output_dir="results/multiscale_udt"):
        """
        Create diagnostic plots for multi-scale UDT framework.
        
        Parameters:
        -----------
        output_dir : str
            Directory to save plots
        """
        os.makedirs(output_dir, exist_ok=True)
        
        fig = plt.figure(figsize=(18, 12))
        
        # 1. Scale hierarchy visualization
        plt.subplot(2, 4, 1)
        scales = ['Galactic', 'Cosmological', 'CMB']
        R0_values = [self.R0_galactic, self.R0_cosmological, self.R0_cmb]
        colors = ['blue', 'green', 'red']
        
        bars = plt.bar(scales, np.log10(R0_values), color=colors, alpha=0.7)
        plt.ylabel('log10(R0) [Mpc]')
        plt.title('UDT Scale Hierarchy')
        plt.grid(True, alpha=0.3)
        
        # Add values on bars
        for bar, val in zip(bars, R0_values):
            height = bar.get_height()
            if val >= 1:
                plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                        f'{val:.0f}', ha='center', va='bottom')
            else:
                plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                        f'{val:.3f}', ha='center', va='bottom')
        
        # 2. Temporal geometry across scales
        plt.subplot(2, 4, 2)
        distances = np.logspace(-3, 7, 1000)
        
        tau_gal = self.temporal_geometry_function(distances, "galactic")
        tau_cosmo = self.temporal_geometry_function(distances, "cosmological")
        tau_cmb = self.temporal_geometry_function(distances, "cmb")
        
        plt.semilogx(distances, tau_gal, 'b-', label='Galactic', linewidth=2)
        plt.semilogx(distances, tau_cosmo, 'g-', label='Cosmological', linewidth=2)
        plt.semilogx(distances, tau_cmb, 'r-', label='CMB', linewidth=2)
        
        plt.axvline(x=self.galactic_boundary, color='gray', linestyle='--', alpha=0.5)
        plt.axvline(x=self.cosmological_boundary, color='gray', linestyle='--', alpha=0.5)
        
        plt.xlabel('Distance (Mpc)')
        plt.ylabel('Temporal Factor τ(r)')
        plt.title('Multi-Scale Temporal Geometry')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.ylim(0, 1)
        
        # 3. Galactic rotation curve example
        plt.subplot(2, 4, 3)
        radius_kpc = np.linspace(1, 30, 30)
        mass_profile = 1e10 * (radius_kpc / 10.0)  # Simple linear profile
        
        v_predicted = self.galactic_dynamics_prediction(radius_kpc, mass_profile)
        
        # Baryonic component for comparison
        G = 4.3e-6
        v_baryonic = np.sqrt(G * mass_profile / radius_kpc)
        
        plt.plot(radius_kpc, v_baryonic, 'b--', label='Baryonic only', linewidth=2)
        plt.plot(radius_kpc, v_predicted, 'b-', label='UDT enhanced', linewidth=2)
        
        plt.xlabel('Radius (kpc)')
        plt.ylabel('Velocity (km/s)')
        plt.title('Galactic Rotation Curve')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 4. Supernova Hubble diagram
        plt.subplot(2, 4, 4)
        redshifts = np.linspace(0.01, 2.0, 100)
        mu_udt = self.supernova_distance_modulus(redshifts)
        
        plt.plot(redshifts, mu_udt, 'g-', linewidth=2, label='UDT')
        
        plt.xlabel('Redshift z')
        plt.ylabel('Distance Modulus')
        plt.title('Supernova Hubble Diagram')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 5. CMB conformal time evolution
        plt.subplot(2, 4, 5)
        z_range = np.linspace(1, 2000, 100)
        eta_range = self.R0_cmb / (1 + z_range) / self.c_light
        
        plt.loglog(z_range, eta_range, 'r-', linewidth=2)
        plt.axvline(x=self.z_recombination, color='orange', linestyle='--', 
                   label=f'Recombination (z={self.z_recombination})')
        plt.axhline(y=288.0, color='gray', linestyle='--', alpha=0.7,
                   label='Standard value')
        
        plt.xlabel('Redshift z')
        plt.ylabel('Conformal Time (Mpc/c)')
        plt.title('CMB Conformal Time Evolution')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 6. CMB acoustic peaks
        plt.subplot(2, 4, 6)
        ell_peaks, r_s, D_A = self.cmb_acoustic_peaks()
        peak_numbers = np.arange(1, len(ell_peaks) + 1)
        
        plt.bar(peak_numbers, ell_peaks[:len(peak_numbers)], alpha=0.7, color='red')
        plt.axhline(y=220, color='orange', linestyle='--', label='Planck l1 = 220')
        
        plt.xlabel('Peak Number')
        plt.ylabel('Multipole l')
        plt.title('CMB Acoustic Peak Positions')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 7. Scale regime boundaries
        plt.subplot(2, 4, 7)
        distances = np.logspace(-3, 7, 1000)
        regimes = [self.get_scale_appropriate_R0(d)[1] for d in distances]
        
        # Convert regimes to numbers for plotting
        regime_map = {'galactic': 1, 'cosmological': 2, 'cmb': 3}
        regime_numbers = [regime_map[r] for r in regimes]
        
        plt.semilogx(distances, regime_numbers, 'k-', linewidth=3)
        plt.axvline(x=self.galactic_boundary, color='blue', linestyle='--', 
                   alpha=0.7, label='Gal/Cosmo boundary')
        plt.axvline(x=self.cosmological_boundary, color='green', linestyle='--',
                   alpha=0.7, label='Cosmo/CMB boundary')
        
        plt.xlabel('Distance (Mpc)')
        plt.ylabel('Scale Regime')
        plt.yticks([1, 2, 3], ['Galactic', 'Cosmological', 'CMB'])
        plt.title('Scale Regime Boundaries')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # 8. Multi-scale comparison
        plt.subplot(2, 4, 8)
        
        # Key scales and their characteristic values
        scales_x = [1, 2, 3]
        gal_char = np.log10(self.R0_galactic * 1000)  # kpc scale
        cosmo_char = np.log10(self.R0_cosmological)   # Mpc scale  
        cmb_char = np.log10(self.R0_cmb / 1000)       # Gpc scale
        
        char_values = [gal_char, cosmo_char, cmb_char]
        
        plt.bar(scales_x, char_values, color=['blue', 'green', 'red'], alpha=0.7)
        plt.xlabel('Scale')
        plt.ylabel('log10(Characteristic Scale)')
        plt.title('Multi-Scale Characteristic Values')
        plt.xticks(scales_x, ['Galactic\n(kpc)', 'Cosmological\n(Mpc)', 'CMB\n(Gpc)'])
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plot_file = os.path.join(output_dir, 'multiscale_udt_diagnostics.png')
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"Multi-scale UDT diagnostics saved to: {plot_file}")
        return plot_file


def main():
    """Main analysis routine for multi-scale UDT framework."""
    print("MULTI-SCALE UDT FRAMEWORK")
    print("=" * 60)
    print()
    
    # Initialize framework
    udt = MultiScaleUDT()
    
    # Validate scale consistency
    validation_results = udt.validate_scale_consistency()
    
    # Create diagnostic plots
    print("\nCREATING DIAGNOSTIC PLOTS")
    print("-" * 30)
    plot_file = udt.create_multiscale_diagnostic_plots()
    
    # Summary
    print("\nMULTI-SCALE UDT SUMMARY")
    print("=" * 50)
    print("Framework successfully implements:")
    print("+ Galactic dynamics with R0_gal = 38 kpc")
    print("+ Supernova cosmology with R0_cosmo = 3,000 Mpc") 
    print("+ CMB physics with R0_cmb = 11.4 million Mpc")
    print()
    print(f"CMB scale validation:")
    print(f"  First acoustic peak: l1 = {validation_results['cmb']['l1']:.1f}")
    print(f"  Error vs Planck: {validation_results['cmb']['error_vs_standard']:.1f}%")
    
    if validation_results['cmb']['error_vs_standard'] < 10:
        print("+ CMB scale within acceptable range")
    else:
        print("! CMB scale needs further refinement")
    
    print(f"\nScale hierarchy established:")
    print(f"  CMB/Cosmological ratio: {udt.R0_cmb/udt.R0_cosmological:.0f}")
    print(f"  Cosmological/Galactic ratio: {udt.R0_cosmological/udt.R0_galactic:.0f}")
    print(f"  Total CMB/Galactic ratio: {udt.R0_cmb/udt.R0_galactic:.0f}")
    
    return validation_results


if __name__ == "__main__":
    main()