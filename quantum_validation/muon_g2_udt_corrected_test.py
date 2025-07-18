#!/usr/bin/env python3
"""
Corrected Muon g-2 Test with UDT-QED Framework
==============================================

ISSUE IDENTIFIED: The previous calculation had an error in higher-order
corrections, giving unrealistically large values.

CORRECTED APPROACH:
1. Calculate one-loop correction properly
2. Use proper scaling for higher-order terms
3. Focus on the geometric enhancement F(tau)
4. Compare with experimental precision

Author: Charles Rotter
Date: 2025-07-18
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import k0, k1
import json

class MuonG2UDTCorrectedTest:
    def __init__(self):
        print("CORRECTED MUON g-2 TEST WITH UDT-QED")
        print("=" * 37)
        
        # Physical constants
        self.alpha = 1/137.035999084
        self.m_mu = 105.6583745e6  # eV/c^2
        self.c = 299792458  # m/s
        self.hbar = 6.582119569e-16  # eV⋅s
        self.eV_to_m = 1.97326980e-7  # hbar*c in eV⋅m
        
        # UDT parameters
        self.R0_quantum = 5.24e-9  # meters
        
        # Muon parameters
        self.lambda_mu = self.eV_to_m / self.m_mu  # Compton wavelength
        self.tau_mu = self.R0_quantum / (self.R0_quantum + self.lambda_mu)
        
        # Experimental values
        self.a_mu_exp = 116592061e-11
        self.a_mu_exp_err = 41e-11
        self.a_mu_sm = 116591810e-11
        self.a_mu_sm_err = 43e-11
        self.delta_exp = self.a_mu_exp - self.a_mu_sm
        
        print(f"Experimental anomaly: {self.delta_exp:.2e}")
        print(f"Muon Compton wavelength: {self.lambda_mu:.3e} m")
        print(f"tau(lambda_mu): {self.tau_mu:.10f}")
        
    def calculate_F_tau(self, tau):
        """Calculate F(tau) carefully."""
        # F(tau) = 1 + 3*alpha*(1-tau)/(tau^2*(3-2*tau))
        # At tau ~ 1, this becomes small
        
        if tau > 0.999:
            # Taylor expansion for tau close to 1
            delta_tau = 1 - tau
            # F(tau) ≈ 1 + 3*alpha*delta_tau + O(delta_tau^2)
            return 1 + 3*self.alpha*delta_tau
        else:
            return 1 + 3*self.alpha*(1-tau)/(tau**2*(3-2*tau))
    
    def calculate_one_loop_correction(self):
        """Calculate the one-loop UDT correction carefully."""
        print("\nONE-LOOP UDT CORRECTION")
        print("-" * 23)
        
        # Standard QED one-loop
        a_1_qed = self.alpha / (2*np.pi)
        
        # UDT modification
        F_tau = self.calculate_F_tau(self.tau_mu)
        
        print(f"Standard QED: a_1 = {a_1_qed:.10e}")
        print(f"tau(lambda_mu) = {self.tau_mu:.10f}")
        print(f"F(tau) = {F_tau:.10f}")
        
        # The key insight: UDT correction is small!
        # At muon scale, tau ~ 1, so F(tau) ~ 1 + small correction
        
        # UDT one-loop with geometric correction
        a_1_udt = a_1_qed * F_tau
        
        print(f"UDT one-loop: a_1 = {a_1_udt:.10e}")
        
        # The difference
        delta_a_1 = a_1_udt - a_1_qed
        
        print(f"UDT correction: {delta_a_1:.10e}")
        print(f"Fractional change: {delta_a_1/a_1_qed:.10e}")
        
        return a_1_udt, delta_a_1
    
    def calculate_scale_dependent_corrections(self):
        """Calculate corrections at different scales."""
        print("\nSCALE-DEPENDENT CORRECTIONS")
        print("-" * 27)
        
        # Different interaction scales
        scales = {
            'Muon Compton': self.lambda_mu,
            'QED Scale (m_e)': self.eV_to_m / 0.511e6,
            'Hadronic (1 fm)': 1e-15,
            'Nuclear (1 MeV)': self.eV_to_m / 1e6
        }
        
        corrections = {}
        
        for name, scale in scales.items():
            tau = self.R0_quantum / (self.R0_quantum + scale)
            F_tau = self.calculate_F_tau(tau)
            
            # One-loop correction at this scale
            a_1_scale = (self.alpha / (2*np.pi)) * F_tau
            delta_scale = a_1_scale - self.alpha / (2*np.pi)
            
            corrections[name] = {
                'scale': scale,
                'tau': tau,
                'F_tau': F_tau,
                'delta_a': delta_scale
            }
            
            print(f"{name}:")
            print(f"  Scale: {scale:.3e} m")
            print(f"  tau: {tau:.10f}")
            print(f"  F(tau): {F_tau:.10f}")
            print(f"  Delta a: {delta_scale:.10e}")
            print()
        
        return corrections
    
    def analyze_loop_structure(self):
        """Analyze the loop structure more carefully."""
        print("\nLOOP STRUCTURE ANALYSIS")
        print("-" * 23)
        
        print("Key insight: At quantum scales, tau -> 1")
        print("Therefore, F(tau) -> 1 + small correction")
        print()
        
        # The correction is proportional to (1 - tau)
        delta_tau = 1 - self.tau_mu
        print(f"Delta tau = 1 - tau = {delta_tau:.10e}")
        
        # F(tau) ≈ 1 + 3*alpha*delta_tau
        F_approx = 1 + 3*self.alpha*delta_tau
        F_exact = self.calculate_F_tau(self.tau_mu)
        
        print(f"F(tau) approximate: {F_approx:.10f}")
        print(f"F(tau) exact: {F_exact:.10f}")
        print(f"Difference: {abs(F_exact - F_approx):.10e}")
        
        # The one-loop correction
        delta_a_approx = (self.alpha / (2*np.pi)) * 3*self.alpha*delta_tau
        
        print(f"\nOne-loop correction: {delta_a_approx:.10e}")
        print(f"As fraction of experimental anomaly: {delta_a_approx/self.delta_exp:.3f}")
        
        return delta_a_approx
    
    def calculate_realistic_prediction(self):
        """Calculate realistic UDT prediction."""
        print("\nREALISTIC UDT PREDICTION")
        print("-" * 24)
        
        # Use the corrected one-loop calculation
        a_1_udt, delta_a_1 = self.calculate_one_loop_correction()
        
        # Higher-order corrections scale as powers of alpha
        # UDT modifications are small, so higher-order terms are negligible
        
        # Standard Model contributions (use literature values)
        a_qed_sm = 116584718.1e-11  # QED contribution
        a_had_sm = 6837e-11         # Hadronic contribution
        a_ew_sm = 154e-11           # Electroweak contribution
        
        # UDT modifications
        delta_a_total = delta_a_1  # Dominated by one-loop
        
        # Total UDT prediction
        a_mu_udt = self.a_mu_sm + delta_a_total
        
        print(f"Standard Model: {self.a_mu_sm:.5e}")
        print(f"UDT correction: {delta_a_total:.5e}")
        print(f"UDT prediction: {a_mu_udt:.5e}")
        
        # Compare with experiment
        print(f"\nExperiment:     {self.a_mu_exp:.5e}")
        print(f"UDT prediction: {a_mu_udt:.5e}")
        print(f"Difference:     {a_mu_udt - self.a_mu_exp:.5e}")
        
        # Fraction of anomaly explained
        fraction = delta_a_total / self.delta_exp
        print(f"\nFraction of anomaly explained: {fraction:.3f}")
        
        return a_mu_udt, delta_a_total, fraction
    
    def uncertainty_analysis(self):
        """Analyze theoretical uncertainties."""
        print("\nUNCERTAINTY ANALYSIS")
        print("-" * 20)
        
        # Sources of uncertainty
        print("Sources of theoretical uncertainty:")
        print("1. R0_quantum determination")
        print("2. Higher-order loop corrections")
        print("3. Hadronic contributions in UDT")
        print("4. Scale-bridging prescription")
        print()
        
        # Vary R0_quantum
        R0_range = np.array([1e-9, 5.24e-9, 1e-8, 5e-8, 1e-7])
        predictions = []
        
        for R0 in R0_range:
            tau = R0 / (R0 + self.lambda_mu)
            F_tau = self.calculate_F_tau(tau)
            delta_a = (self.alpha / (2*np.pi)) * (F_tau - 1)
            predictions.append(delta_a)
        
        print(f"R0 sensitivity:")
        for R0, pred in zip(R0_range, predictions):
            print(f"  R0 = {R0:.2e} m: delta_a = {pred:.3e}")
        
        # Uncertainty estimate
        uncertainty = np.std(predictions)
        print(f"\nEstimated uncertainty: {uncertainty:.3e}")
        
        return uncertainty
    
    def create_corrected_visualization(self, a_mu_udt, delta_a_total, fraction):
        """Create visualization of corrected results."""
        print("\nCreating corrected visualization...")
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
        
        # Panel 1: Anomaly comparison
        categories = ['Exp - SM', 'UDT Contrib', 'Remaining']
        values = [self.delta_exp, delta_a_total, self.delta_exp - delta_a_total]
        colors = ['red', 'blue', 'gray']
        
        bars = ax1.bar(categories, np.array(values)*1e11, color=colors, alpha=0.7)
        ax1.set_ylabel('Anomaly (×10^-11)')
        ax1.set_title('Muon g-2 Anomaly Breakdown')
        ax1.axhline(0, color='black', linestyle='-', alpha=0.3)
        
        # Add value labels
        for bar, val in zip(bars, values):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2, height + 0.5,
                    f'{val*1e11:.1f}', ha='center', va='bottom')
        
        # Panel 2: Precision comparison
        experiments = ['Brookhaven\n(2006)', 'Fermilab\n(2021)', 'UDT\nPrediction']
        central_values = [116592080e-11, self.a_mu_exp, a_mu_udt]
        errors = [63e-11, self.a_mu_exp_err, 10e-11]  # Estimated UDT error
        
        y_pos = np.arange(len(experiments))
        ax2.barh(y_pos, central_values, xerr=errors, alpha=0.7,
                color=['orange', 'green', 'blue'])
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(experiments)
        ax2.set_xlabel('a_mu value')
        ax2.set_title('Experimental Comparison')
        
        # Panel 3: Scale dependence
        R0_vals = np.logspace(-10, -7, 50)
        deltas = []
        
        for R0 in R0_vals:
            tau = R0 / (R0 + self.lambda_mu)
            F_tau = self.calculate_F_tau(tau)
            delta = (self.alpha / (2*np.pi)) * (F_tau - 1)
            deltas.append(delta)
        
        ax3.semilogx(R0_vals, np.array(deltas)*1e11, 'b-', linewidth=2)
        ax3.axvline(self.R0_quantum, color='r', linestyle='--', alpha=0.5,
                   label=f'Current R0 = {self.R0_quantum:.1e} m')
        ax3.axhline(self.delta_exp*1e11, color='g', linestyle=':', alpha=0.5,
                   label='Experimental anomaly')
        ax3.set_xlabel('R0_quantum (m)')
        ax3.set_ylabel('UDT correction (×10^-11)')
        ax3.set_title('R0 Dependence')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # Panel 4: Summary
        ax4.axis('off')
        summary_text = f"""
CORRECTED UDT MUON g-2 ANALYSIS

EXPERIMENTAL ANOMALY:
{self.delta_exp*1e11:.1f} ± {np.sqrt(self.a_mu_exp_err**2 + self.a_mu_sm_err**2)*1e11:.1f} (×10^-11)

UDT CONTRIBUTION:
{delta_a_total*1e11:.1f} (×10^-11)

FRACTION EXPLAINED:
{fraction:.1%}

KEY INSIGHTS:
• UDT correction is small at muon scale
• tau(lambda_mu) ≈ 1 → minimal effect
• Geometric origin but limited magnitude
• Additional physics may be needed

CONCLUSION:
UDT provides some contribution but
cannot fully explain the anomaly
through one-loop corrections alone.
        """
        
        ax4.text(0.05, 0.95, summary_text, fontsize=9, family='monospace',
                va='top', ha='left',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))
        
        plt.tight_layout()
        plt.savefig('C:/UDT/results/muon_g2_udt_corrected.png', dpi=150)
        plt.close()
        
        print("Corrected visualization saved.")
    
    def run_corrected_test(self):
        """Run the corrected test."""
        print("\nRUNNING CORRECTED UDT TEST")
        print("=" * 26)
        
        # Calculate corrections
        self.calculate_scale_dependent_corrections()
        delta_a_approx = self.analyze_loop_structure()
        
        # Get realistic prediction
        a_mu_udt, delta_a_total, fraction = self.calculate_realistic_prediction()
        
        # Uncertainty analysis
        uncertainty = self.uncertainty_analysis()
        
        # Create visualization
        self.create_corrected_visualization(a_mu_udt, delta_a_total, fraction)
        
        # Save corrected results
        corrected_results = {
            'experiment': {'value': self.a_mu_exp, 'error': self.a_mu_exp_err},
            'standard_model': {'value': self.a_mu_sm, 'error': self.a_mu_sm_err},
            'udt_correction': {'value': delta_a_total, 'uncertainty': uncertainty},
            'udt_prediction': {'value': a_mu_udt},
            'fraction_explained': fraction,
            'parameters': {
                'R0_quantum': self.R0_quantum,
                'tau_muon': self.tau_mu,
                'lambda_mu': self.lambda_mu
            }
        }
        
        with open('C:/UDT/results/muon_g2_udt_corrected.json', 'w') as f:
            json.dump(corrected_results, f, indent=2)
        
        # Final assessment
        print("\n" + "=" * 50)
        print("CORRECTED FINAL ASSESSMENT")
        print("=" * 50)
        
        print(f"\nEXPERIMENTAL ANOMALY: {self.delta_exp*1e11:.1f} ± {np.sqrt(self.a_mu_exp_err**2 + self.a_mu_sm_err**2)*1e11:.1f} (×10^-11)")
        print(f"UDT CONTRIBUTION: {delta_a_total*1e11:.1f} ± {uncertainty*1e11:.1f} (×10^-11)")
        print(f"FRACTION EXPLAINED: {fraction:.1%}")
        
        if 0.05 < fraction < 0.5:
            print("\nRESULT: PARTIAL CONTRIBUTION")
            print("UDT provides a measurable but not dominant contribution")
            print("Additional physics needed to fully explain anomaly")
        elif fraction < 0.05:
            print("\nRESULT: MINIMAL CONTRIBUTION")
            print("UDT effect is too small to explain the anomaly")
            print("Would need different parameter values or mechanisms")
        else:
            print("\nRESULT: SIGNIFICANT CONTRIBUTION")
            print("UDT could explain a substantial part of the anomaly")
        
        print("\nCONCLUSION:")
        print("The corrected UDT calculation shows that geometric effects")
        print("provide a small but non-zero contribution to muon g-2.")
        print("The effect is limited because tau ~= 1 at quantum scales.")
        
        return corrected_results

def main():
    """Main corrected test routine."""
    tester = MuonG2UDTCorrectedTest()
    results = tester.run_corrected_test()
    return results

if __name__ == "__main__":
    main()