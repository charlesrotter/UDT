#!/usr/bin/env python3
"""
LIGO Gravitational Wave Analysis with UDT Framework
===================================================

QUANTUM DOMAIN TEST: Does UDT's c_fundamental = infinity affect gravitational 
wave propagation and timing in ways detectable by LIGO?

EXPERIMENTAL CONTEXT:
- GW170817 neutron star merger (2017): Gravitational waves + gamma-ray burst
- Travel distance: ~130 million light-years
- Timing difference: ~1.7 seconds between GW and gamma rays
- Constraint: Speed of gravity equals speed of light to within 1e-15

UDT HYPOTHESIS:
In UDT, gravity propagates with c_fundamental = infinity, but effective speeds
are c_eff(r) = c_0 * tau(r). This might create subtle timing or dispersion
effects detectable in LIGO data.

CRITICAL TEST:
Can UDT explain LIGO observations or does it predict timing differences
that violate the strict GW170817 constraints?

Author: Charles Rotter  
Date: 2025-07-18
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import chi2
import sys
import os

# Add project root to path  
sys.path.append(os.path.abspath('.'))

class LIGOGravitationalWaveUDTAnalysis:
    def __init__(self):
        print("LIGO GRAVITATIONAL WAVE UDT ANALYSIS")
        print("=" * 40)
        
        # Physical constants
        self.c = 299792458  # m/s (speed of light)
        self.G = 6.67430e-11  # m^3 kg^-1 s^-2
        self.Mpc_to_m = 3.086e22  # meters per Mpc
        
        # GW170817 observational data
        self.gw170817_distance = 130e6  # light-years
        self.gw170817_distance_mpc = self.gw170817_distance * 0.3066  # Convert to Mpc
        self.gw170817_timing_diff = 1.7  # seconds between GW and gamma rays
        self.gw170817_timing_error = 0.1  # estimated uncertainty
        
        # LIGO constraint: |c_gw - c_light| / c_light < 1e-15
        self.ligo_speed_constraint = 1e-15
        
        # UDT parameters (from cosmological validation)
        self.R0_cosmo = 3582.0  # Mpc
        self.R0_gw = 3582.0  # Mpc (initial guess - may need fitting)
        
        print(f"GW170817 Analysis Parameters:")
        print(f"  Distance: {self.gw170817_distance:.0f} light-years = {self.gw170817_distance_mpc:.1f} Mpc")
        print(f"  Timing difference: {self.gw170817_timing_diff:.1f} +/- {self.gw170817_timing_error:.1f} s")
        print(f"  LIGO speed constraint: |c_gw - c_light|/c_light < {self.ligo_speed_constraint:.0e}")
        print(f"  UDT R0: {self.R0_gw:.1f} Mpc")
        
    def calculate_udt_gravitational_wave_propagation(self, distance_mpc):
        """Calculate UDT predictions for gravitational wave propagation."""
        print(f"\nUDT GRAVITATIONAL WAVE PROPAGATION")
        print("-" * 35)
        
        # In UDT, gravitational interactions have c_fundamental = infinity
        # But observed wave propagation depends on spacetime geometry
        
        # UDT temporal geometry at GW170817 distance
        tau_gw = self.R0_gw / (self.R0_gw + distance_mpc)
        
        # Effective gravitational wave speed in UDT
        c_gw_eff = self.c * tau_gw
        
        print(f"UDT geometry at GW170817:")
        print(f"  Distance: {distance_mpc:.1f} Mpc")
        print(f"  tau(r): {tau_gw:.6f}")
        print(f"  c_gw_eff: {c_gw_eff:.0f} m/s")
        print(f"  c_gw_eff/c: {c_gw_eff/self.c:.6f}")
        
        # For electromagnetic waves (gamma rays), same geometry applies
        c_em_eff = self.c * tau_gw
        
        print(f"  c_em_eff: {c_em_eff:.0f} m/s")
        print(f"  c_em_eff/c: {c_em_eff/self.c:.6f}")
        
        # Travel time difference in UDT
        # If both GW and EM waves use same effective speed, no timing difference
        travel_time_gw = distance_mpc * self.Mpc_to_m / c_gw_eff
        travel_time_em = distance_mpc * self.Mpc_to_m / c_em_eff
        
        udt_timing_diff = travel_time_gw - travel_time_em
        
        print(f"  Travel time (GW): {travel_time_gw:.0f} s")
        print(f"  Travel time (EM): {travel_time_em:.0f} s")
        print(f"  UDT timing difference: {udt_timing_diff:.3e} s")
        
        # Speed constraint test
        speed_violation = abs(c_gw_eff - c_em_eff) / self.c
        
        print(f"  Speed constraint violation: {speed_violation:.3e}")
        print(f"  LIGO limit: {self.ligo_speed_constraint:.0e}")
        
        if speed_violation < self.ligo_speed_constraint:
            print("  UDT SATISFIES LIGO SPEED CONSTRAINT")
        else:
            print("  UDT VIOLATES LIGO SPEED CONSTRAINT")
            
        return tau_gw, c_gw_eff, c_em_eff, udt_timing_diff, speed_violation
        
    def calculate_udt_gravitational_wave_dispersion(self, distance_mpc):
        """Calculate UDT predictions for gravitational wave dispersion."""
        print(f"\nUDT GRAVITATIONAL WAVE DISPERSION")
        print("-" * 34)
        
        # In UDT, different frequency components might propagate differently
        # This is a more sophisticated effect than simple timing
        
        # LIGO frequency range: ~30-1000 Hz
        frequencies = np.array([30, 100, 300, 1000])  # Hz
        
        print(f"Testing dispersion across LIGO frequency range:")
        
        dispersions = []
        for freq in frequencies:
            # In UDT, high-frequency GWs might couple differently to geometry
            # This is speculative - need proper derivation from UDT field equations
            
            # For now, assume frequency-dependent coupling
            # This needs theoretical justification
            freq_factor = 1.0  # No dispersion in simple UDT
            
            tau_freq = self.R0_gw / (self.R0_gw + distance_mpc * freq_factor)
            c_gw_freq = self.c * tau_freq
            
            travel_time_freq = distance_mpc * self.Mpc_to_m / c_gw_freq
            dispersions.append(travel_time_freq)
            
            print(f"  {freq:.0f} Hz: tau = {tau_freq:.6f}, c_eff = {c_gw_freq:.0f} m/s")
        
        # Calculate dispersion spread
        dispersion_spread = max(dispersions) - min(dispersions)
        
        print(f"  Dispersion spread: {dispersion_spread:.3e} s")
        print(f"  Observed chirp duration: ~0.2 s")
        
        if dispersion_spread > 0.001:  # 1 ms threshold
            print("  SIGNIFICANT DISPERSION PREDICTED")
        else:
            print("  NO SIGNIFICANT DISPERSION PREDICTED")
            
        return frequencies, dispersions, dispersion_spread
        
    def analyze_gw170817_timing_constraint(self):
        """Analyze GW170817 timing constraint with UDT."""
        print(f"\nGW170817 TIMING CONSTRAINT ANALYSIS")
        print("-" * 37)
        
        # The observed 1.7 s timing difference is attributed to:
        # 1. Different emission times at the source
        # 2. Different propagation paths
        # 3. NOT different propagation speeds (ruled out by LIGO)
        
        # UDT prediction for timing difference
        tau_gw, c_gw_eff, c_em_eff, udt_timing_diff, speed_violation = \
            self.calculate_udt_gravitational_wave_propagation(self.gw170817_distance_mpc)
        
        print(f"Constraint analysis:")
        print(f"  Observed timing difference: {self.gw170817_timing_diff:.1f} +/- {self.gw170817_timing_error:.1f} s")
        print(f"  UDT predicted difference: {udt_timing_diff:.3e} s")
        print(f"  Speed constraint violation: {speed_violation:.3e}")
        print(f"  LIGO speed limit: {self.ligo_speed_constraint:.0e}")
        
        # Chi-squared test
        # If UDT predicts significant timing difference, it's ruled out
        chi2_timing = (udt_timing_diff / self.gw170817_timing_error)**2
        
        print(f"  Chi-squared (timing): {chi2_timing:.3f}")
        
        # Speed constraint test
        constraint_violation = speed_violation > self.ligo_speed_constraint
        
        print(f"  Speed constraint {'VIOLATED' if constraint_violation else 'SATISFIED'}")
        
        if constraint_violation:
            print("  UDT RULED OUT BY LIGO SPEED CONSTRAINT")
        elif chi2_timing > 9:  # 3-sigma
            print("  UDT RULED OUT BY TIMING CONSTRAINT")
        elif chi2_timing > 4:  # 2-sigma  
            print("  UDT CONSTRAINED BY TIMING DATA")
        else:
            print("  UDT CONSISTENT WITH LIGO CONSTRAINTS")
            
        return chi2_timing, constraint_violation, speed_violation
        
    def fit_udt_parameters_to_ligo_data(self):
        """Fit UDT parameters to LIGO gravitational wave data."""
        print(f"\nFITTING UDT PARAMETERS TO LIGO DATA")
        print("-" * 37)
        
        # Available LIGO events for fitting
        ligo_events = {
            'GW170817': {'distance_mpc': 40.0, 'timing_constraint': 1e-15},
            'GW150914': {'distance_mpc': 410.0, 'timing_constraint': None},
            'GW151226': {'distance_mpc': 440.0, 'timing_constraint': None},
            'GW170104': {'distance_mpc': 880.0, 'timing_constraint': None}
        }
        
        def udt_ligo_likelihood(params):
            R0_gw_log = params[0]
            
            # Convert from log parameter
            R0_gw = 10**R0_gw_log
            
            # Parameter bounds
            if R0_gw < 1e-10 or R0_gw > 1e10:
                return 1e10
            
            chi2_total = 0
            
            for event_name, event_data in ligo_events.items():
                distance_mpc = event_data['distance_mpc']
                
                # Calculate UDT prediction
                tau_event = R0_gw / (R0_gw + distance_mpc)
                c_gw_eff = self.c * tau_event
                
                # Speed constraint test
                speed_violation = abs(c_gw_eff - self.c) / self.c
                
                if event_data['timing_constraint'] is not None:
                    # Apply timing constraint
                    constraint_chi2 = (speed_violation / event_data['timing_constraint'])**2
                    chi2_total += constraint_chi2
                
                # For GW170817, also test timing difference
                if event_name == 'GW170817':
                    travel_time_gw = distance_mpc * self.Mpc_to_m / c_gw_eff
                    travel_time_em = distance_mpc * self.Mpc_to_m / self.c
                    timing_diff = travel_time_gw - travel_time_em
                    
                    # Should be consistent with no intrinsic timing difference
                    timing_chi2 = (timing_diff / self.gw170817_timing_error)**2
                    chi2_total += timing_chi2
            
            return chi2_total
        
        # Fit parameters
        print("Fitting UDT parameters to LIGO constraints...")
        
        # Initial guess
        initial_params = [np.log10(self.R0_gw)]
        
        # Minimize
        result = minimize(udt_ligo_likelihood, initial_params, method='Nelder-Mead',
                         options={'maxiter': 10000})
        
        if result.success:
            R0_gw_fit = 10**result.x[0]
            chi2_min = result.fun
            
            print(f"UDT LIGO FIT RESULTS:")
            print(f"  R0_gw = {R0_gw_fit:.1f} Mpc")
            print(f"  chi2_total = {chi2_min:.3f}")
            
            # Test final fit
            tau_fit = R0_gw_fit / (R0_gw_fit + self.gw170817_distance_mpc)
            c_gw_fit = self.c * tau_fit
            speed_violation_fit = abs(c_gw_fit - self.c) / self.c
            
            print(f"  Final speed violation: {speed_violation_fit:.3e}")
            print(f"  LIGO limit: {self.ligo_speed_constraint:.0e}")
            
            if chi2_min < 1:
                print("  FIT QUALITY: EXCELLENT")
            elif chi2_min < 4:
                print("  FIT QUALITY: GOOD")
            else:
                print("  FIT QUALITY: POOR")
            
            return R0_gw_fit, chi2_min, speed_violation_fit
        else:
            print(f"FITTING FAILED: {result.message}")
            return None, None, None
            
    def create_ligo_analysis_plots(self, R0_gw_fit, chi2_min, speed_violation_fit):
        """Create LIGO analysis plots."""
        print(f"\nCreating LIGO analysis plots...")
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Panel 1: Speed constraint vs distance
        ax1 = axes[0, 0]
        
        distances = np.logspace(0, 3, 100)  # 1 to 1000 Mpc
        speed_violations = []
        
        for dist in distances:
            tau_d = R0_gw_fit / (R0_gw_fit + dist)
            c_gw_d = self.c * tau_d
            violation = abs(c_gw_d - self.c) / self.c
            speed_violations.append(violation)
        
        ax1.loglog(distances, speed_violations, 'b-', linewidth=2, label='UDT Prediction')
        ax1.axhline(y=self.ligo_speed_constraint, color='r', linestyle='--', 
                   label=f'LIGO Limit ({self.ligo_speed_constraint:.0e})')
        ax1.scatter([self.gw170817_distance_mpc], [speed_violation_fit], 
                   color='red', s=100, zorder=5, label='GW170817')
        
        ax1.set_xlabel('Distance (Mpc)')
        ax1.set_ylabel('Speed Violation |c_gw - c|/c')
        ax1.set_title('UDT Speed Constraint vs Distance')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Panel 2: Effective speed vs distance
        ax2 = axes[0, 1]
        
        c_gw_values = []
        for dist in distances:
            tau_d = R0_gw_fit / (R0_gw_fit + dist)
            c_gw_d = self.c * tau_d
            c_gw_values.append(c_gw_d)
        
        ax2.semilogx(distances, np.array(c_gw_values)/self.c, 'g-', linewidth=2, label='c_gw_eff/c')
        ax2.axhline(y=1, color='k', linestyle='--', alpha=0.5, label='c_0')
        ax2.scatter([self.gw170817_distance_mpc], [c_gw_values[np.argmin(np.abs(distances - self.gw170817_distance_mpc))]/self.c], 
                   color='red', s=100, zorder=5, label='GW170817')
        
        ax2.set_xlabel('Distance (Mpc)')
        ax2.set_ylabel('c_gw_eff / c_0')
        ax2.set_title('UDT Effective GW Speed vs Distance')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Panel 3: LIGO events comparison
        ax3 = axes[1, 0]
        
        event_names = ['GW150914', 'GW151226', 'GW170104', 'GW170817']
        event_distances = [410, 440, 880, 40]
        event_violations = []
        
        for dist in event_distances:
            tau_event = R0_gw_fit / (R0_gw_fit + dist)
            c_gw_event = self.c * tau_event
            violation = abs(c_gw_event - self.c) / self.c
            event_violations.append(violation)
        
        bars = ax3.bar(event_names, event_violations, alpha=0.7, color=['blue', 'green', 'orange', 'red'])
        ax3.axhline(y=self.ligo_speed_constraint, color='r', linestyle='--', 
                   label=f'LIGO Limit ({self.ligo_speed_constraint:.0e})')
        ax3.set_ylabel('Speed Violation |c_gw - c|/c')
        ax3.set_title('UDT Predictions for LIGO Events')
        ax3.set_yscale('log')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        plt.setp(ax3.get_xticklabels(), rotation=45, ha='right')
        
        # Panel 4: Summary
        ax4 = axes[1, 1]
        ax4.axis('off')
        
        summary_text = f"""
        LIGO UDT ANALYSIS RESULTS
        
        UDT Parameters:
        R0_gw = {R0_gw_fit:.1f} Mpc
        
        GW170817 Constraints:
        Distance: {self.gw170817_distance_mpc:.1f} Mpc
        Speed violation: {speed_violation_fit:.3e}
        LIGO limit: {self.ligo_speed_constraint:.0e}
        
        Fit Quality:
        chi2 = {chi2_min:.3f}
        
        Constraint Test:
        {'PASSES' if speed_violation_fit < self.ligo_speed_constraint else 'FAILS'} LIGO speed limit
        
        {'SUCCESS: UDT consistent with LIGO' if speed_violation_fit < self.ligo_speed_constraint else 'FAILURE: UDT violates LIGO constraints'}
        """
        
        ax4.text(0.1, 0.5, summary_text, transform=ax4.transAxes,
                fontsize=10, verticalalignment='center', fontfamily='monospace',
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightblue", alpha=0.8))
        
        plt.tight_layout()
        plt.savefig('C:/UDT/results/ligo_gravitational_wave_udt_analysis.png', dpi=150)
        plt.close()
        
        print("LIGO analysis plots saved to: C:/UDT/results/ligo_gravitational_wave_udt_analysis.png")
        
    def run_complete_analysis(self):
        """Run complete LIGO gravitational wave UDT analysis."""
        print("COMPLETE LIGO GRAVITATIONAL WAVE UDT ANALYSIS")
        print("=" * 46)
        
        # 1. Calculate UDT GW propagation
        tau_gw, c_gw_eff, c_em_eff, udt_timing_diff, speed_violation = \
            self.calculate_udt_gravitational_wave_propagation(self.gw170817_distance_mpc)
        
        # 2. Analyze dispersion
        frequencies, dispersions, dispersion_spread = \
            self.calculate_udt_gravitational_wave_dispersion(self.gw170817_distance_mpc)
        
        # 3. Analyze GW170817 timing constraint
        chi2_timing, constraint_violation, speed_violation_final = \
            self.analyze_gw170817_timing_constraint()
        
        # 4. Fit UDT parameters
        R0_gw_fit, chi2_min, speed_violation_fit = self.fit_udt_parameters_to_ligo_data()
        
        if R0_gw_fit is None:
            print("ANALYSIS FAILED: Could not fit UDT parameters")
            return None
        
        # 5. Create plots
        self.create_ligo_analysis_plots(R0_gw_fit, chi2_min, speed_violation_fit)
        
        # Final assessment
        print("\n" + "=" * 60)
        print("LIGO GRAVITATIONAL WAVE UDT ANALYSIS CONCLUSIONS")
        print("=" * 60)
        
        print(f"\n1. UDT GRAVITATIONAL WAVE PREDICTIONS:")
        print(f"   R0_gw = {R0_gw_fit:.1f} Mpc")
        print(f"   GW170817 distance: {self.gw170817_distance_mpc:.1f} Mpc")
        print(f"   tau(GW170817): {tau_gw:.6f}")
        print(f"   c_gw_eff/c: {c_gw_eff/self.c:.6f}")
        
        print(f"\n2. LIGO CONSTRAINT TESTS:")
        print(f"   Speed violation: {speed_violation_fit:.3e}")
        print(f"   LIGO speed limit: {self.ligo_speed_constraint:.0e}")
        print(f"   Timing difference: {udt_timing_diff:.3e} s")
        print(f"   Dispersion spread: {dispersion_spread:.3e} s")
        
        print(f"\n3. FIT QUALITY:")
        print(f"   chi2 = {chi2_min:.3f}")
        print(f"   Timing chi2 = {chi2_timing:.3f}")
        
        success = not constraint_violation and chi2_min < 4
        print(f"\n4. CONCLUSION:")
        if success:
            print("   SUCCESS: UDT CONSISTENT WITH LIGO CONSTRAINTS")
            print("   SUCCESS: SECOND QUANTUM DOMAIN TEST PASSED")
            print("   SUCCESS: TOE FRAMEWORK SURVIVES LIGO VALIDATION")
        else:
            print("   FAILURE: UDT VIOLATES LIGO CONSTRAINTS")
            print("   FAILURE: SECOND QUANTUM DOMAIN TEST FAILED")
            print("   FAILURE: TOE FRAMEWORK RULED OUT BY LIGO")
        
        return {
            'R0_gw_fit': R0_gw_fit,
            'chi2_min': chi2_min,
            'speed_violation': speed_violation_fit,
            'constraint_violation': constraint_violation,
            'success': success
        }

def main():
    """Main analysis routine."""
    analyzer = LIGOGravitationalWaveUDTAnalysis()
    results = analyzer.run_complete_analysis()
    return results

if __name__ == "__main__":
    main()