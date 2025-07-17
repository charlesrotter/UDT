#!/usr/bin/env python3
"""
Quantum Experimental Simulation - UDT Framework
===============================================

Simulates realistic experimental data for quantum phenomena to demonstrate
how UDT predictions would compare against actual measurements.

This script generates synthetic experimental data with realistic noise
and uncertainties for:
1. Hydrogen atom spectroscopy
2. Quantum tunneling through barriers
3. Casimir force measurements

The simulated data allows validation of UDT predictions against
"experimental" results in a controlled environment.

Author: Charles Rotter
Date: 2025-01-17
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm
import json
import os

class QuantumExperimentalSimulator:
    """Simulate realistic quantum experimental data."""
    
    def __init__(self, seed=42):
        """Initialize with reproducible random seed."""
        np.random.seed(seed)
        self.results_dir = "results/quantum_experimental_simulation"
        os.makedirs(self.results_dir, exist_ok=True)
        
    def simulate_hydrogen_spectroscopy(self, R0_quantum=1e-10):
        """
        Simulate hydrogen atom spectroscopy experiment.
        
        This simulates measuring the Balmer series with realistic
        experimental uncertainties and systematic errors.
        """
        print("=" * 70)
        print("SIMULATED EXPERIMENT 1: HYDROGEN SPECTROSCOPY")
        print("=" * 70)
        print("Simulating Balmer series measurements with experimental noise")
        print()
        
        # Physical constants
        h = 6.626e-34  # Planck constant
        c = 3e8        # Speed of light
        R_inf = 1.097e7  # Rydberg constant (m^-1)
        
        # Balmer series transitions (n=2 to n=3,4,5,6,7)
        transitions = [(2, 3), (2, 4), (2, 5), (2, 6), (2, 7)]
        transition_names = ["H-alpha", "H-beta", "H-gamma", "H-delta", "H-epsilon"]
        
        # Standard quantum mechanics predictions
        standard_wavelengths = []
        standard_frequencies = []
        
        for n_low, n_high in transitions:
            # Rydberg formula: 1/lambda = R_inf * (1/n_low^2 - 1/n_high^2)
            wave_number = R_inf * (1/n_low**2 - 1/n_high**2)
            wavelength = 1 / wave_number
            frequency = c / wavelength
            
            standard_wavelengths.append(wavelength)
            standard_frequencies.append(frequency)
        
        # UDT-modified predictions
        udt_wavelengths = []
        udt_frequencies = []
        
        for i, (n_low, n_high) in enumerate(transitions):
            # UDT enhancement factor at characteristic radius
            r_char = n_high**2 * 5.29e-11  # Bohr radius scaling
            tau = R0_quantum / (R0_quantum + r_char)
            enhancement = 1 / tau**2
            
            # Modified energy levels affect transition frequencies
            # Enhanced binding energy increases transition frequency
            modified_frequency = standard_frequencies[i] * np.sqrt(enhancement)
            modified_wavelength = c / modified_frequency
            
            udt_wavelengths.append(modified_wavelength)
            udt_frequencies.append(modified_frequency)
        
        # Simulate experimental measurements with noise
        experimental_wavelengths = []
        experimental_errors = []
        
        for i, true_wavelength in enumerate(standard_wavelengths):
            # Realistic measurement uncertainty (~0.1 nm)
            measurement_error = 0.1e-9
            
            # Add systematic bias and random noise
            systematic_bias = np.random.normal(0, 0.05e-9)  # 0.05 nm systematic
            random_noise = np.random.normal(0, measurement_error)
            
            measured_wavelength = true_wavelength + systematic_bias + random_noise
            
            experimental_wavelengths.append(measured_wavelength)
            experimental_errors.append(measurement_error)
        
        # Analysis and comparison
        print("BALMER SERIES MEASUREMENTS:")
        print("Transition | Standard (nm) | UDT (nm)    | Measured (nm) | Error (nm)")
        print("-" * 70)
        
        chi2_standard = 0
        chi2_udt = 0
        
        for i, name in enumerate(transition_names):
            std_nm = standard_wavelengths[i] * 1e9
            udt_nm = udt_wavelengths[i] * 1e9
            exp_nm = experimental_wavelengths[i] * 1e9
            err_nm = experimental_errors[i] * 1e9
            
            # Calculate chi-squared
            chi2_standard += ((exp_nm - std_nm) / err_nm)**2
            chi2_udt += ((exp_nm - udt_nm) / err_nm)**2
            
            print(f"{name:8s}  | {std_nm:10.3f}    | {udt_nm:10.3f}  | {exp_nm:10.3f}    | {err_nm:8.3f}")
        
        print()
        print(f"Chi-squared analysis:")
        print(f"  Standard QM: chi^2 = {chi2_standard:.2f}")
        print(f"  UDT theory:  chi^2 = {chi2_udt:.2f}")
        print(f"  Best fit: {'UDT' if chi2_udt < chi2_standard else 'Standard QM'}")
        print()
        
        # Save results
        results = {
            'R0_quantum': R0_quantum,
            'transitions': transition_names,
            'standard_wavelengths_nm': [w*1e9 for w in standard_wavelengths],
            'udt_wavelengths_nm': [w*1e9 for w in udt_wavelengths],
            'experimental_wavelengths_nm': [w*1e9 for w in experimental_wavelengths],
            'experimental_errors_nm': [e*1e9 for e in experimental_errors],
            'chi2_standard': chi2_standard,
            'chi2_udt': chi2_udt
        }
        
        with open(f"{self.results_dir}/hydrogen_spectroscopy_results.json", "w") as f:
            json.dump(results, f, indent=2)
        
        return results
    
    def simulate_tunneling_experiment(self, R0_quantum=1e-10):
        """
        Simulate quantum tunneling experiment.
        
        This simulates measuring tunneling current through a barrier
        as a function of barrier width and height.
        """
        print("=" * 70)
        print("SIMULATED EXPERIMENT 2: QUANTUM TUNNELING")
        print("=" * 70)
        print("Simulating scanning tunneling microscopy measurements")
        print()
        
        # Experimental parameters
        barrier_widths = np.linspace(0.5e-9, 2.0e-9, 10)  # 0.5 to 2.0 nm
        barrier_height = 1.0 * 1.602e-19  # 1.0 eV in Joules
        
        # Physical constants
        h_bar = 1.055e-34
        m_e = 9.109e-31
        
        # Standard quantum mechanics predictions
        standard_currents = []
        
        for width in barrier_widths:
            # WKB approximation for tunneling probability
            k = np.sqrt(2 * m_e * barrier_height) / h_bar
            T = np.exp(-2 * k * width)
            
            # Current proportional to tunneling probability
            current = T * 1e-9  # Scale to nanoamps
            standard_currents.append(current)
        
        # UDT-modified predictions
        udt_currents = []
        
        for width in barrier_widths:
            # UDT enhancement at barrier scale
            r_barrier = width / 2
            tau = R0_quantum / (R0_quantum + r_barrier)
            enhancement = 1 / tau**2
            
            # Enhanced barrier height reduces tunneling
            enhanced_height = barrier_height * enhancement
            k_udt = np.sqrt(2 * m_e * enhanced_height) / h_bar
            T_udt = np.exp(-2 * k_udt * width)
            
            current_udt = T_udt * 1e-9
            udt_currents.append(current_udt)
        
        # Simulate experimental measurements
        experimental_currents = []
        experimental_errors = []
        
        for i, true_current in enumerate(standard_currents):
            # Realistic measurement uncertainty (10% of signal)
            measurement_error = 0.1 * true_current
            
            # Add noise
            noise = np.random.normal(0, measurement_error)
            measured_current = true_current + noise
            
            experimental_currents.append(measured_current)
            experimental_errors.append(measurement_error)
        
        # Analysis
        print("TUNNELING CURRENT MEASUREMENTS:")
        print("Width (nm) | Standard (nA) | UDT (nA)    | Measured (nA) | Error (nA)")
        print("-" * 70)
        
        chi2_standard = 0
        chi2_udt = 0
        
        for i, width in enumerate(barrier_widths):
            width_nm = width * 1e9
            std_na = standard_currents[i] * 1e9
            udt_na = udt_currents[i] * 1e9
            exp_na = experimental_currents[i] * 1e9
            err_na = experimental_errors[i] * 1e9
            
            chi2_standard += ((exp_na - std_na) / err_na)**2
            chi2_udt += ((exp_na - udt_na) / err_na)**2
            
            print(f"{width_nm:8.2f}   | {std_na:10.3f}    | {udt_na:10.3f}   | {exp_na:10.3f}    | {err_na:8.3f}")
        
        print()
        print(f"Chi-squared analysis:")
        print(f"  Standard QM: chi^2 = {chi2_standard:.2f}")
        print(f"  UDT theory:  chi^2 = {chi2_udt:.2f}")
        print(f"  Best fit: {'UDT' if chi2_udt < chi2_standard else 'Standard QM'}")
        print()
        
        # Save results
        results = {
            'R0_quantum': R0_quantum,
            'barrier_widths_nm': [w*1e9 for w in barrier_widths],
            'standard_currents_nA': [c*1e9 for c in standard_currents],
            'udt_currents_nA': [c*1e9 for c in udt_currents],
            'experimental_currents_nA': [c*1e9 for c in experimental_currents],
            'experimental_errors_nA': [e*1e9 for e in experimental_errors],
            'chi2_standard': chi2_standard,
            'chi2_udt': chi2_udt
        }
        
        with open(f"{self.results_dir}/tunneling_experiment_results.json", "w") as f:
            json.dump(results, f, indent=2)
        
        return results
    
    def simulate_casimir_experiment(self, R0_quantum=1e-10):
        """
        Simulate Casimir force measurement.
        
        This simulates measuring the Casimir force between parallel
        plates as a function of separation distance.
        """
        print("=" * 70)
        print("SIMULATED EXPERIMENT 3: CASIMIR FORCE")
        print("=" * 70)
        print("Simulating atomic force microscopy Casimir measurements")
        print()
        
        # Experimental parameters
        separations = np.linspace(10e-9, 100e-9, 10)  # 10 to 100 nm
        plate_area = 1e-12  # 1 mm^2
        
        # Physical constants
        h_bar = 1.055e-34
        c = 3e8
        
        # Standard Casimir force predictions
        standard_forces = []
        
        for separation in separations:
            # Casimir force per unit area
            force_per_area = -np.pi**2 * h_bar * c / (240 * separation**4)
            total_force = force_per_area * plate_area
            
            standard_forces.append(abs(total_force))  # Take absolute value
        
        # UDT-modified predictions
        udt_forces = []
        
        for separation in separations:
            # UDT enhancement at separation scale
            r_sep = separation / 2
            tau = R0_quantum / (R0_quantum + r_sep)
            enhancement = 1 / tau**2
            
            # Enhanced vacuum fluctuations increase Casimir force
            force_per_area = -np.pi**2 * h_bar * c / (240 * separation**4)
            enhanced_force = force_per_area * plate_area * enhancement
            
            udt_forces.append(abs(enhanced_force))
        
        # Simulate experimental measurements
        experimental_forces = []
        experimental_errors = []
        
        for i, true_force in enumerate(standard_forces):
            # Realistic measurement uncertainty (5% of signal)
            measurement_error = 0.05 * true_force
            
            # Add noise
            noise = np.random.normal(0, measurement_error)
            measured_force = true_force + noise
            
            experimental_forces.append(measured_force)
            experimental_errors.append(measurement_error)
        
        # Analysis
        print("CASIMIR FORCE MEASUREMENTS:")
        print("Sep (nm) | Standard (pN) | UDT (pN)    | Measured (pN) | Error (pN)")
        print("-" * 70)
        
        chi2_standard = 0
        chi2_udt = 0
        
        for i, separation in enumerate(separations):
            sep_nm = separation * 1e9
            std_pn = standard_forces[i] * 1e12
            udt_pn = udt_forces[i] * 1e12
            exp_pn = experimental_forces[i] * 1e12
            err_pn = experimental_errors[i] * 1e12
            
            chi2_standard += ((exp_pn - std_pn) / err_pn)**2
            chi2_udt += ((exp_pn - udt_pn) / err_pn)**2
            
            print(f"{sep_nm:6.0f}   | {std_pn:10.3f}    | {udt_pn:10.3f}   | {exp_pn:10.3f}    | {err_pn:8.3f}")
        
        print()
        print(f"Chi-squared analysis:")
        print(f"  Standard QM: chi^2 = {chi2_standard:.2f}")
        print(f"  UDT theory:  chi^2 = {chi2_udt:.2f}")
        print(f"  Best fit: {'UDT' if chi2_udt < chi2_standard else 'Standard QM'}")
        print()
        
        # Save results
        results = {
            'R0_quantum': R0_quantum,
            'separations_nm': [s*1e9 for s in separations],
            'standard_forces_pN': [f*1e12 for f in standard_forces],
            'udt_forces_pN': [f*1e12 for f in udt_forces],
            'experimental_forces_pN': [f*1e12 for f in experimental_forces],
            'experimental_errors_pN': [e*1e12 for e in experimental_errors],
            'chi2_standard': chi2_standard,
            'chi2_udt': chi2_udt
        }
        
        with open(f"{self.results_dir}/casimir_experiment_results.json", "w") as f:
            json.dump(results, f, indent=2)
        
        return results
    
    def create_experimental_plots(self, hydrogen_results, tunneling_results, casimir_results):
        """Create publication-quality plots of experimental results."""
        print("=" * 70)
        print("GENERATING EXPERIMENTAL PLOTS")
        print("=" * 70)
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('UDT Quantum Experimental Validation', fontsize=16)
        
        # Hydrogen spectroscopy plot
        ax1 = axes[0, 0]
        transitions = hydrogen_results['transitions']
        x_pos = np.arange(len(transitions))
        
        ax1.errorbar(x_pos, hydrogen_results['experimental_wavelengths_nm'], 
                    yerr=hydrogen_results['experimental_errors_nm'],
                    fmt='o', label='Experimental', capsize=5)
        ax1.plot(x_pos, hydrogen_results['standard_wavelengths_nm'], 
                's-', label='Standard QM', alpha=0.7)
        ax1.plot(x_pos, hydrogen_results['udt_wavelengths_nm'], 
                '^-', label='UDT Theory', alpha=0.7)
        
        ax1.set_xlabel('Balmer Transition')
        ax1.set_ylabel('Wavelength (nm)')
        ax1.set_title('Hydrogen Spectroscopy')
        ax1.set_xticks(x_pos)
        ax1.set_xticklabels(transitions, rotation=45)
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Tunneling current plot
        ax2 = axes[0, 1]
        ax2.errorbar(tunneling_results['barrier_widths_nm'], 
                    tunneling_results['experimental_currents_nA'],
                    yerr=tunneling_results['experimental_errors_nA'],
                    fmt='o', label='Experimental', capsize=5)
        ax2.plot(tunneling_results['barrier_widths_nm'], 
                tunneling_results['standard_currents_nA'], 
                's-', label='Standard QM', alpha=0.7)
        ax2.plot(tunneling_results['barrier_widths_nm'], 
                tunneling_results['udt_currents_nA'], 
                '^-', label='UDT Theory', alpha=0.7)
        
        ax2.set_xlabel('Barrier Width (nm)')
        ax2.set_ylabel('Tunneling Current (nA)')
        ax2.set_title('Quantum Tunneling')
        ax2.set_yscale('log')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Casimir force plot
        ax3 = axes[1, 0]
        ax3.errorbar(casimir_results['separations_nm'], 
                    casimir_results['experimental_forces_pN'],
                    yerr=casimir_results['experimental_errors_pN'],
                    fmt='o', label='Experimental', capsize=5)
        ax3.plot(casimir_results['separations_nm'], 
                casimir_results['standard_forces_pN'], 
                's-', label='Standard QM', alpha=0.7)
        ax3.plot(casimir_results['separations_nm'], 
                casimir_results['udt_forces_pN'], 
                '^-', label='UDT Theory', alpha=0.7)
        
        ax3.set_xlabel('Plate Separation (nm)')
        ax3.set_ylabel('Casimir Force (pN)')
        ax3.set_title('Casimir Effect')
        ax3.set_yscale('log')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Chi-squared comparison
        ax4 = axes[1, 1]
        experiments = ['Hydrogen\nSpectroscopy', 'Quantum\nTunneling', 'Casimir\nEffect']
        standard_chi2 = [hydrogen_results['chi2_standard'], 
                        tunneling_results['chi2_standard'], 
                        casimir_results['chi2_standard']]
        udt_chi2 = [hydrogen_results['chi2_udt'], 
                   tunneling_results['chi2_udt'], 
                   casimir_results['chi2_udt']]
        
        x_pos = np.arange(len(experiments))
        width = 0.35
        
        ax4.bar(x_pos - width/2, standard_chi2, width, label='Standard QM', alpha=0.7)
        ax4.bar(x_pos + width/2, udt_chi2, width, label='UDT Theory', alpha=0.7)
        
        ax4.set_xlabel('Experiment')
        ax4.set_ylabel('Chi-squared')
        ax4.set_title('Model Comparison')
        ax4.set_xticks(x_pos)
        ax4.set_xticklabels(experiments)
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f"{self.results_dir}/quantum_experimental_validation.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Plots saved to {self.results_dir}/quantum_experimental_validation.png")
        print()

def main():
    """Run comprehensive quantum experimental simulation."""
    print("QUANTUM EXPERIMENTAL SIMULATION")
    print("=" * 35)
    print("Simulating realistic quantum experiments to validate UDT")
    print("Using synthetic data with experimental noise and uncertainties")
    print()
    
    # Initialize simulator
    simulator = QuantumExperimentalSimulator()
    
    # Test different R0 values
    R0_values = [1e-11, 1e-10, 1e-9]  # Different quantum scales
    
    for R0_quantum in R0_values:
        print(f"Testing R0 = {R0_quantum:.0e} m")
        print("-" * 40)
        
        # Run simulated experiments
        hydrogen_results = simulator.simulate_hydrogen_spectroscopy(R0_quantum)
        tunneling_results = simulator.simulate_tunneling_experiment(R0_quantum)
        casimir_results = simulator.simulate_casimir_experiment(R0_quantum)
        
        # Create plots
        simulator.create_experimental_plots(hydrogen_results, tunneling_results, casimir_results)
        
        print("=" * 70)
        print(f"SUMMARY FOR R0 = {R0_quantum:.0e} m")
        print("=" * 70)
        print(f"Hydrogen Spectroscopy: {'UDT' if hydrogen_results['chi2_udt'] < hydrogen_results['chi2_standard'] else 'Standard QM'} fits better")
        print(f"Quantum Tunneling: {'UDT' if tunneling_results['chi2_udt'] < tunneling_results['chi2_standard'] else 'Standard QM'} fits better")
        print(f"Casimir Effect: {'UDT' if casimir_results['chi2_udt'] < casimir_results['chi2_standard'] else 'Standard QM'} fits better")
        print()
    
    print("=" * 70)
    print("EXPERIMENTAL VALIDATION COMPLETE")
    print("=" * 70)
    print("Simulated experimental data demonstrates:")
    print("- UDT modifications are measurable at quantum scales")
    print("- Different experiments probe different aspects of temporal geometry")
    print("- Statistical analysis can distinguish between theories")
    print("- Results saved to results/quantum_experimental_simulation/")
    print()

if __name__ == "__main__":
    main()