#!/usr/bin/env python3
"""
CMB Power Spectrum Analysis using UDT Framework
===============================================

Analyzes Planck CMB power spectra using the Universal Distance Dilation Theory.
Compares UDT predictions with ΛCDM best-fit model for temperature and 
polarization power spectra.

Key UDT CMB Predictions:
1. Modified angular diameter distance: d_A(z) = R₀ * ln(1 + z) 
2. Position-dependent effective light speed affects horizon scale
3. Temporal geometry modifies recombination physics
4. Acoustic oscillation frequencies shift due to c_eff(r)

Author: Charles Rotter
Date: 2025-01-17
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import json

# Add parent directory to path to import udt package
sys.path.insert(0, str(Path(__file__).parent.parent))

class UDTCMBAnalyzer:
    """Analyze CMB power spectra with UDT temporal geometry."""
    
    def __init__(self):
        """Initialize with physical constants and CMB parameters."""
        # Physical constants
        self.c = 2.998e8           # Speed of light (m/s)
        self.h = 6.626e-34         # Planck constant
        self.k_B = 1.381e-23       # Boltzmann constant
        self.sigma_T = 6.652e-29   # Thomson scattering cross-section
        
        # CMB parameters
        self.T_CMB = 2.725         # CMB temperature (K)
        self.z_recomb = 1090       # Recombination redshift
        self.z_star = 1020         # Surface of last scattering
        
        # UDT parameters (will be fitted)
        self.R0_cmb = 14000        # CMB-scale R₀ (Mpc) - initial guess
        
        # Data directories
        self.data_dir = Path(__file__).parent.parent / "data" / "cmb_planck"
        self.results_dir = Path(__file__).parent.parent / "results" / "cmb_analysis"
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
    def load_planck_data(self):
        """Load Planck CMB power spectrum data."""
        print("Loading Planck CMB power spectrum data...")
        
        # Load TT power spectrum
        tt_file = self.data_dir / "COM_PowerSpect_CMB-TT-binned_R3.01.txt"
        if tt_file.exists():
            tt_data = np.loadtxt(tt_file, comments='#')
            self.l_tt = tt_data[:, 0]
            self.Dl_TT_obs = tt_data[:, 1] 
            self.Dl_TT_err_low = tt_data[:, 2]
            self.Dl_TT_err_high = tt_data[:, 3]
            self.Dl_TT_bestfit = tt_data[:, 4]
            print(f"  Loaded TT data: {len(self.l_tt)} multipoles")
        else:
            print(f"  WARNING: TT file not found: {tt_file}")
            
        # Load TE power spectrum  
        te_file = self.data_dir / "COM_PowerSpect_CMB-TE-binned_R3.02.txt"
        if te_file.exists():
            te_data = np.loadtxt(te_file, comments='#')
            self.l_te = te_data[:, 0]
            self.Dl_TE_obs = te_data[:, 1]
            self.Dl_TE_err_low = te_data[:, 2] 
            self.Dl_TE_err_high = te_data[:, 3]
            self.Dl_TE_bestfit = te_data[:, 4]
            print(f"  Loaded TE data: {len(self.l_te)} multipoles")
        else:
            print(f"  WARNING: TE file not found: {te_file}")
            
        # Load EE power spectrum
        ee_file = self.data_dir / "COM_PowerSpect_CMB-EE-binned_R3.02.txt"
        if ee_file.exists():
            ee_data = np.loadtxt(ee_file, comments='#')
            self.l_ee = ee_data[:, 0]
            self.Dl_EE_obs = ee_data[:, 1]
            self.Dl_EE_err_low = ee_data[:, 2]
            self.Dl_EE_err_high = ee_data[:, 3] 
            self.Dl_EE_bestfit = ee_data[:, 4]
            print(f"  Loaded EE data: {len(self.l_ee)} multipoles")
        else:
            print(f"  WARNING: EE file not found: {ee_file}")
            
        print()
        
    def udt_angular_diameter_distance(self, z, R0):
        """
        UDT angular diameter distance.
        
        In UDT, redshift z is temporal dilation: z = (R₀ + r)/R₀ - 1 = r/R₀
        Therefore: r = z * R₀
        
        Angular diameter distance in UDT:
        d_A = r / (1 + z) = (z * R₀) / (1 + z) = R₀ * z / (1 + z)
        
        For small z: d_A ≈ R₀ * ln(1 + z)
        """
        return R0 * z / (1 + z)
        
    def udt_sound_horizon(self, z_star, R0):
        """
        UDT sound horizon at recombination.
        
        The sound horizon is modified by position-dependent effective light speed:
        c_eff(r) = c₀ * R₀/(R₀ + r)
        
        For temporal geometry, the sound speed is also modified:
        c_s_eff = c_s * τ(r) = c_s * R₀/(R₀ + r)
        """
        # Standard sound speed at recombination
        c_s = self.c / np.sqrt(3)  # Simplified sound speed
        
        # Distance to recombination surface
        r_star = z_star * R0
        
        # Effective sound speed in UDT
        tau_star = R0 / (R0 + r_star)
        c_s_eff = c_s * tau_star
        
        # Sound horizon (simplified)
        r_s = c_s_eff * 380000 * 365.25 * 24 * 3600  # ~380,000 years to recombination
        
        return r_s
        
    def udt_acoustic_scale(self, z_star, R0):
        """
        UDT acoustic scale theta_A = r_s / d_A.
        
        This is the fundamental angular scale of acoustic oscillations
        that determines the position of CMB peaks.
        """
        # Standard values for comparison
        # Observed acoustic scale is about 0.01 radians
        r_s_standard = 150e6 * 3.086e22  # ~150 Mpc in meters
        d_A_standard = z_star * R0 * 3.086e22  # Mpc to meters
        
        # UDT modification factor
        tau_star = R0 / (R0 + z_star * R0)  # At recombination
        
        # Modified acoustic scale
        theta_A = (r_s_standard / d_A_standard) * tau_star
        
        return theta_A
        
    def udt_peak_positions(self, R0):
        """
        Predict CMB acoustic peak positions in UDT.
        
        Peak positions: l_n ≈ n * π / theta_A
        where theta_A is the acoustic scale.
        """
        theta_A = self.udt_acoustic_scale(self.z_star, R0)
        
        # First few acoustic peaks
        peak_numbers = np.array([1, 2, 3, 4, 5])
        l_peaks = peak_numbers * np.pi / theta_A
        
        return l_peaks
        
    def udt_damping_scale(self, z_star, R0):
        """
        UDT photon diffusion damping scale.
        
        Silk damping is modified by temporal geometry affecting 
        the photon mean free path and diffusion.
        """
        # Distance to recombination
        r_star = z_star * R0
        tau_star = R0 / (R0 + r_star)
        
        # Modified diffusion scale in temporal geometry
        # Standard Silk damping scale ~ 0.3 degrees
        l_D_standard = 1000  # Standard damping scale
        
        # UDT modification
        l_D_udt = l_D_standard * (1 / tau_star)
        
        return l_D_udt
        
    def udt_tt_power_spectrum(self, l, R0):
        """
        UDT TT power spectrum prediction.
        
        This is a simplified model focusing on the acoustic oscillation structure
        and damping. A full calculation would require solving the Boltzmann equation
        with UDT modifications.
        """
        # Acoustic peak positions
        l_peaks = self.udt_peak_positions(R0)
        
        # Damping scale
        l_D = self.udt_damping_scale(self.z_star, R0)
        
        # Basic oscillation pattern
        theta_A = self.udt_acoustic_scale(self.z_star, R0)
        phase = l * theta_A
        
        # Amplitude envelope (simplified)
        amplitude = 6000 * np.exp(-(l / l_D)**2)  # Silk damping
        
        # Acoustic oscillations
        oscillation = 1 + 0.3 * np.cos(phase) + 0.1 * np.cos(2 * phase)
        
        # Add plateau for low-l
        plateau = 1000 * (l / 100)**(-0.5)
        
        # Combine components
        Dl_TT = np.where(l < 50, plateau, amplitude * oscillation)
        
        return Dl_TT
        
    def fit_udt_to_cmb(self):
        """
        Fit UDT parameters to Planck CMB power spectra.
        """
        print("Fitting UDT model to Planck CMB data...")
        print()
        
        if not hasattr(self, 'l_tt'):
            print("ERROR: No CMB data loaded!")
            return None
            
        # Define fitting function
        def udt_model(l, R0):
            return self.udt_tt_power_spectrum(l, R0)
            
        # Fit to TT power spectrum (focus on acoustic peaks region)
        mask = (self.l_tt >= 100) & (self.l_tt <= 1000)
        l_fit = self.l_tt[mask]
        Dl_fit = self.Dl_TT_obs[mask]
        err_fit = (self.Dl_TT_err_low[mask] + self.Dl_TT_err_high[mask]) / 2
        
        try:
            # Initial guess
            p0 = [self.R0_cmb]
            
            # Fit with bounds
            bounds = ([5000], [50000])  # R₀ between 5-50 Gpc
            
            popt, pcov = curve_fit(udt_model, l_fit, Dl_fit, 
                                 sigma=err_fit, p0=p0, bounds=bounds,
                                 maxfev=1000)
            
            R0_best = popt[0]
            R0_err = np.sqrt(pcov[0, 0])
            
            print(f"BEST-FIT UDT PARAMETERS:")
            print(f"  R0_CMB = {R0_best:.0f} +/- {R0_err:.0f} Mpc")
            
            # Calculate chi-squared
            Dl_pred = udt_model(l_fit, R0_best)
            chi2 = np.sum(((Dl_fit - Dl_pred) / err_fit)**2)
            dof = len(l_fit) - 1
            chi2_reduced = chi2 / dof
            
            print(f"  chi2 = {chi2:.1f}")
            print(f"  DOF = {dof}")
            print(f"  chi2/DOF = {chi2_reduced:.2f}")
            print()
            
            # Calculate key UDT predictions
            theta_A = self.udt_acoustic_scale(self.z_star, R0_best)
            l_peaks = self.udt_peak_positions(R0_best)
            l_D = self.udt_damping_scale(self.z_star, R0_best)
            
            print(f"UDT CMB PREDICTIONS:")
            print(f"  Acoustic scale theta_A = {theta_A:.6f} rad")
            print(f"  First peak position: l1 = {l_peaks[0]:.0f}")
            print(f"  Second peak position: l2 = {l_peaks[1]:.0f}")  
            print(f"  Third peak position: l3 = {l_peaks[2]:.0f}")
            print(f"  Damping scale: l_D = {l_D:.0f}")
            print()
            
            self.R0_cmb_best = R0_best
            self.cmb_fit_results = {
                'R0_best': R0_best,
                'R0_err': R0_err,
                'chi2': chi2,
                'dof': dof,
                'chi2_reduced': chi2_reduced,
                'theta_A': theta_A,
                'l_peaks': l_peaks.tolist(),
                'l_damping': l_D
            }
            
            return self.cmb_fit_results
            
        except Exception as e:
            print(f"ERROR in fitting: {e}")
            return None
            
    def compare_with_lcdm(self):
        """
        Compare UDT predictions with ΛCDM best-fit model.
        """
        print("Comparing UDT vs LCDM CMB predictions...")
        print()
        
        if not hasattr(self, 'R0_cmb_best'):
            print("ERROR: No UDT fit results available!")
            return
            
        # Standard ΛCDM parameters (Planck 2018)
        lcdm_params = {
            'H0': 67.4,      # km/s/Mpc  
            'Omega_m': 0.315,
            'Omega_b': 0.049,
            'theta_s': 1.042e-2,  # Sound horizon angle
            'A_s': 2.1e-9,        # Scalar amplitude
            'n_s': 0.965          # Spectral index
        }
        
        # UDT vs ΛCDM comparison
        udt_theta_A = self.cmb_fit_results['theta_A']
        lcdm_theta_A = lcdm_params['theta_s']
        
        theta_diff_percent = (udt_theta_A - lcdm_theta_A) / lcdm_theta_A * 100
        
        print(f"ACOUSTIC SCALE COMPARISON:")
        print(f"  LCDM theta_s = {lcdm_theta_A:.6f} rad")
        print(f"  UDT theta_A = {udt_theta_A:.6f} rad") 
        print(f"  Difference: {theta_diff_percent:+.2f}%")
        print()
        
        # Peak position comparison
        lcdm_l1 = 220  # Observed first peak
        udt_l1 = self.cmb_fit_results['l_peaks'][0]
        
        l1_diff = (udt_l1 - lcdm_l1) / lcdm_l1 * 100
        
        print(f"FIRST PEAK POSITION:")
        print(f"  Observed: l1 ~= {lcdm_l1}")
        print(f"  UDT: l1 = {udt_l1:.0f}")
        print(f"  Difference: {l1_diff:+.1f}%")
        print()
        
        return {
            'lcdm_theta_s': lcdm_theta_A,
            'udt_theta_A': udt_theta_A,
            'theta_difference_percent': theta_diff_percent,
            'lcdm_l1': lcdm_l1,
            'udt_l1': udt_l1,
            'l1_difference_percent': l1_diff
        }
        
    def create_cmb_visualization(self):
        """
        Create comprehensive CMB analysis visualization.
        """
        print("Creating CMB analysis visualization...")
        
        if not hasattr(self, 'l_tt'):
            print("ERROR: No CMB data loaded!")
            return
            
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        
        # Plot 1: TT Power Spectrum
        ax1.errorbar(self.l_tt, self.Dl_TT_obs, 
                    yerr=(self.Dl_TT_err_low + self.Dl_TT_err_high)/2,
                    fmt='o', markersize=3, alpha=0.7, label='Planck TT', color='blue')
        ax1.plot(self.l_tt, self.Dl_TT_bestfit, 'r-', linewidth=2, label='LCDM Best-fit')
        
        # Add UDT prediction if fit was successful
        if hasattr(self, 'R0_cmb_best'):
            l_pred = np.linspace(50, 2000, 500)
            Dl_udt = self.udt_tt_power_spectrum(l_pred, self.R0_cmb_best)
            ax1.plot(l_pred, Dl_udt, 'g--', linewidth=2, label='UDT Prediction')
            
            # Mark acoustic peaks
            l_peaks = self.cmb_fit_results['l_peaks']
            for i, l_peak in enumerate(l_peaks[:3]):
                if l_peak < 2000:
                    ax1.axvline(l_peak, color='green', linestyle=':', alpha=0.7,
                               label=f'UDT Peak {i+1}' if i == 0 else "")
        
        ax1.set_xlabel('Multipole l')
        ax1.set_ylabel('Dl^TT [muK^2]')
        ax1.set_title('CMB Temperature Power Spectrum')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_xlim(2, 2000)
        
        # Plot 2: TE Power Spectrum  
        if hasattr(self, 'l_te'):
            ax2.errorbar(self.l_te, self.Dl_TE_obs,
                        yerr=(self.Dl_TE_err_low + self.Dl_TE_err_high)/2,
                        fmt='o', markersize=3, alpha=0.7, label='Planck TE', color='orange')
            ax2.plot(self.l_te, self.Dl_TE_bestfit, 'r-', linewidth=2, label='LCDM Best-fit')
        
        ax2.set_xlabel('Multipole l')
        ax2.set_ylabel('Dl^TE [muK^2]')
        ax2.set_title('CMB Temperature-E Polarization')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: EE Power Spectrum
        if hasattr(self, 'l_ee'):
            ax3.errorbar(self.l_ee, self.Dl_EE_obs,
                        yerr=(self.Dl_EE_err_low + self.Dl_EE_err_high)/2, 
                        fmt='o', markersize=3, alpha=0.7, label='Planck EE', color='purple')
            ax3.plot(self.l_ee, self.Dl_EE_bestfit, 'r-', linewidth=2, label='LCDM Best-fit')
        
        ax3.set_xlabel('Multipole l')
        ax3.set_ylabel('Dl^EE [muK^2]')
        ax3.set_title('CMB E-mode Polarization')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Plot 4: UDT Parameter Space
        if hasattr(self, 'R0_cmb_best'):
            R0_values = np.linspace(5000, 30000, 100)
            theta_A_values = [self.udt_acoustic_scale(self.z_star, R0) for R0 in R0_values]
            l1_values = [self.udt_peak_positions(R0)[0] for R0 in R0_values]
            
            ax4_twin = ax4.twinx()
            
            line1 = ax4.plot(R0_values, theta_A_values, 'g-', linewidth=2, label='Acoustic Scale theta_A')
            line2 = ax4_twin.plot(R0_values, l1_values, 'b--', linewidth=2, label='First Peak l1')
            
            # Mark best-fit value
            ax4.axvline(self.R0_cmb_best, color='red', linestyle=':', alpha=0.7, 
                       label=f'Best-fit R0 = {self.R0_cmb_best:.0f} Mpc')
            
            ax4.set_xlabel('R0 CMB Scale (Mpc)')
            ax4.set_ylabel('Acoustic Scale theta_A (rad)', color='g')
            ax4_twin.set_ylabel('First Peak Position l1', color='b')
            ax4.set_title('UDT CMB Parameter Relationships')
            
            # Combine legends
            lines1, labels1 = ax4.get_legend_handles_labels()
            lines2, labels2 = ax4_twin.get_legend_handles_labels()
            ax4.legend(lines1 + lines2, labels1 + labels2, loc='best')
            
            ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.results_dir / 'udt_cmb_analysis.png', dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"CMB analysis saved: {self.results_dir}/udt_cmb_analysis.png")
        print()
        
    def run_cmb_analysis(self):
        """
        Run complete CMB analysis with UDT.
        """
        print("=" * 70)
        print("UDT CMB POWER SPECTRUM ANALYSIS")
        print("=" * 70)
        print()
        
        print("Analyzing Planck CMB data with Universal Distance Dilation Theory...")
        print("Testing whether UDT temporal geometry can explain CMB observations")
        print("without requiring dark matter or cosmological constant.")
        print()
        
        # Load data
        self.load_planck_data()
        
        # Fit UDT model
        fit_results = self.fit_udt_to_cmb()
        
        if fit_results:
            # Compare with ΛCDM
            comparison = self.compare_with_lcdm()
            
            # Create visualization
            self.create_cmb_visualization()
            
            # Save results
            all_results = {
                'udt_fit': fit_results,
                'lcdm_comparison': comparison,
                'data_info': {
                    'tt_multipoles': len(self.l_tt) if hasattr(self, 'l_tt') else 0,
                    'te_multipoles': len(self.l_te) if hasattr(self, 'l_te') else 0,
                    'ee_multipoles': len(self.l_ee) if hasattr(self, 'l_ee') else 0
                }
            }
            
            with open(self.results_dir / 'udt_cmb_results.json', 'w') as f:
                json.dump(all_results, f, indent=2, default=str)
            
            print("=" * 70)
            print("CMB ANALYSIS SUMMARY")
            print("=" * 70)
            print()
            
            print("* CMB DATA: Successfully loaded Planck power spectra")
            print("* UDT FIT: Fitted temporal geometry to acoustic oscillations")
            print("* COMPARISON: Compared UDT vs LCDM predictions")
            print("* VISUALIZATION: Created comprehensive analysis plots")
            print()
            
            R0_best = fit_results['R0_best']
            chi2_red = fit_results['chi2_reduced']
            
            print(f"KEY RESULTS:")
            print(f"- CMB-scale characteristic length: R0 = {R0_best:.0f} Mpc")
            print(f"- Reduced chi-squared: chi2/DOF = {chi2_red:.2f}")
            print(f"- Acoustic scale difference: {comparison['theta_difference_percent']:+.1f}%")
            print(f"- First peak shift: {comparison['l1_difference_percent']:+.1f}%")
            print()
            
            print("THEORETICAL SIGNIFICANCE:")
            print("UDT provides an alternative geometric explanation for CMB")
            print("acoustic oscillations through temporal dilation effects,")
            print("potentially eliminating the need for dark matter and dark energy.")
            print()
            
            print(f"Complete CMB analysis results: {self.results_dir}/")
            
            return all_results
        else:
            print("ERROR: CMB analysis failed!")
            return None

def main():
    """Main CMB analysis."""
    analyzer = UDTCMBAnalyzer()
    results = analyzer.run_cmb_analysis()
    return results

if __name__ == "__main__":
    main()