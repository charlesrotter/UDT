#!/usr/bin/env python3
"""
Honest Comparison: Fundamental UDT vs Curve-Fitting vs Newtonian
================================================================

This provides an honest comparison between:
1. Fundamental UDT (no fitted parameters)
2. Curve-fitting UDT (2 fitted parameters per galaxy) 
3. Newtonian gravity baseline
4. Random 2-parameter model

This reveals whether UDT success comes from physics or parameter optimization.

Author: Charles Rotter
Date: 2025-01-17
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import pandas as pd

class HonestGalacticComparison:
    """Compare fundamental UDT vs curve-fitting honestly."""
    
    def __init__(self):
        self.c0 = 299792.458  # km/s
        
    def create_realistic_test_galaxy(self, seed=42):
        """Create a realistic test galaxy with known dark matter profile."""
        np.random.seed(seed)
        
        # Realistic galaxy parameters
        r_max = 20.0  # kpc
        n_points = 15
        radius = np.linspace(1.0, r_max, n_points)
        
        # Stellar mass profile (exponential disk)
        M_stellar_total = 5e10  # Solar masses
        r_scale = 3.0  # kpc
        M_stellar = M_stellar_total * (1 - np.exp(-radius/r_scale) * (1 + radius/r_scale))
        
        # Dark matter profile (NFW)
        M_dm_total = 5e11  # Solar masses  
        r_s = 10.0  # kpc
        M_dm = M_dm_total * (np.log(1 + radius/r_s) - radius/r_s / (1 + radius/r_s))
        
        # Total mass and Newtonian velocity
        M_total = M_stellar + M_dm
        GM = 4.3e-6 * M_total  # km^2/s^2
        v_true = np.sqrt(GM / radius)
        
        # Add realistic observational errors
        v_error = 5.0 + 0.1 * v_true  # 5 km/s + 10% error
        v_observed = v_true + np.random.normal(0, v_error)
        
        return {
            'radius': radius,
            'velocity_true': v_true,
            'velocity_observed': v_observed,
            'velocity_error': v_error,
            'M_stellar': M_stellar,
            'M_dm': M_dm,
            'M_total': M_total
        }
    
    def fundamental_udt_prediction(self, radius, M_stellar, R0_theory=20.0):
        """
        Fundamental UDT prediction with NO fitted parameters.
        
        R0 comes from theory, not fitting.
        """
        velocities = []
        
        for i, r in enumerate(radius):
            # Only stellar mass (no dark matter in UDT)
            M_r = M_stellar[i]
            
            # Temporal dilation
            tau = R0_theory / (R0_theory + r)
            
            # Newtonian baseline
            GM = 4.3e-6 * M_r
            if GM / r <= 0:
                v_newton = 0
            else:
                v_newton = np.sqrt(GM / r)
            
            # UDT enhancement from temporal geometry
            # This is derived from metric, not ad-hoc
            enhancement = 1 / tau  # Simplified from full metric calculation
            
            v_udt = v_newton * enhancement
            velocities.append(v_udt)
        
        return np.array(velocities)
    
    def curve_fitting_udt(self, radius, velocity, velocity_error, M_stellar):
        """
        Curve-fitting UDT with 2 fitted parameters per galaxy.
        
        This is what the old analysis was actually doing.
        """
        def udt_model(params):
            R0_fit, V_scale = params
            
            # Ad-hoc velocity model (what old code actually used)
            v_base = V_scale * np.sqrt(radius / (radius + R0_fit/3))
            
            # Ad-hoc enhancement
            tau = R0_fit / (R0_fit + radius)
            enhancement = 1 / tau**2
            
            v_model = v_base * np.sqrt(enhancement)
            return v_model
        
        def chi2_function(params):
            if params[0] <= 0 or params[1] <= 0:
                return 1e10
            
            v_model = udt_model(params)
            chi2 = np.sum(((velocity - v_model) / velocity_error)**2)
            return chi2
        
        # Fit parameters
        result = minimize(chi2_function, [20.0, 150.0], 
                         bounds=[(5, 100), (50, 300)],
                         method='L-BFGS-B')
        
        if result.success:
            R0_fit, V_scale_fit = result.x
            v_fitted = udt_model(result.x)
            chi2_final = result.fun
            
            return {
                'success': True,
                'R0_fit': R0_fit,
                'V_scale_fit': V_scale_fit,
                'velocity_fitted': v_fitted,
                'chi2': chi2_final,
                'chi2_dof': chi2_final / (len(radius) - 2)
            }
        else:
            return {'success': False}
    
    def newtonian_baseline(self, radius, M_stellar):
        """Pure Newtonian prediction using only stellar mass."""
        GM = 4.3e-6 * M_stellar
        v_newton = np.sqrt(np.maximum(GM / radius, 0))
        return v_newton
    
    def random_two_parameter_model(self, radius, velocity, velocity_error):
        """
        Random 2-parameter model to test if any model can fit.
        
        Uses completely unphysical functional form.
        """
        def random_model(params):
            a, b = params
            # Completely arbitrary functional form
            v_model = a * (radius**0.3) * np.exp(-radius/b) + 50
            return v_model
        
        def chi2_function(params):
            if params[0] <= 0 or params[1] <= 0:
                return 1e10
            
            v_model = random_model(params)
            chi2 = np.sum(((velocity - v_model) / velocity_error)**2)
            return chi2
        
        # Fit parameters
        result = minimize(chi2_function, [100.0, 10.0], 
                         bounds=[(1, 1000), (1, 100)],
                         method='L-BFGS-B')
        
        if result.success:
            return {
                'success': True,
                'velocity_fitted': random_model(result.x),
                'chi2_dof': result.fun / (len(radius) - 2)
            }
        else:
            return {'success': False}
    
    def run_honest_comparison(self, n_galaxies=10):
        """Run honest comparison across multiple test galaxies."""
        
        print("HONEST GALACTIC DYNAMICS COMPARISON")
        print("="*60)
        print("Testing on realistic galaxies with known dark matter")
        print("Comparing:")
        print("1. Fundamental UDT (NO fitting)")
        print("2. Curve-fitting UDT (2 parameters)")  
        print("3. Newtonian gravity (stellar mass only)")
        print("4. Random 2-parameter model")
        print("="*60)
        print()
        
        results = []
        
        for i in range(n_galaxies):
            print(f"Galaxy {i+1:2d}: ", end="")
            
            # Create test galaxy
            galaxy = self.create_realistic_test_galaxy(seed=42+i)
            
            # 1. Fundamental UDT (no fitting)
            v_udt_fundamental = self.fundamental_udt_prediction(
                galaxy['radius'], galaxy['M_stellar'])
            
            residuals_fund = galaxy['velocity_observed'] - v_udt_fundamental
            rms_fund = np.sqrt(np.mean(residuals_fund**2))
            chi2_fund = np.sum((residuals_fund / galaxy['velocity_error'])**2)
            chi2_dof_fund = chi2_fund / (len(galaxy['radius']) - 0)  # No fitted params
            
            # 2. Curve-fitting UDT
            fit_result = self.curve_fitting_udt(
                galaxy['radius'], galaxy['velocity_observed'], 
                galaxy['velocity_error'], galaxy['M_stellar'])
            
            # 3. Newtonian baseline
            v_newton = self.newtonian_baseline(galaxy['radius'], galaxy['M_stellar'])
            residuals_newton = galaxy['velocity_observed'] - v_newton
            rms_newton = np.sqrt(np.mean(residuals_newton**2))
            chi2_newton = np.sum((residuals_newton / galaxy['velocity_error'])**2)
            chi2_dof_newton = chi2_newton / (len(galaxy['radius']) - 0)
            
            # 4. Random model
            random_result = self.random_two_parameter_model(
                galaxy['radius'], galaxy['velocity_observed'], galaxy['velocity_error'])
            
            # Store results
            result = {
                'galaxy': i+1,
                'fundamental_udt_rms': rms_fund,
                'fundamental_udt_chi2_dof': chi2_dof_fund,
                'curve_fit_success': fit_result['success'] if fit_result else False,
                'curve_fit_chi2_dof': fit_result['chi2_dof'] if fit_result['success'] else np.inf,
                'newtonian_rms': rms_newton,
                'newtonian_chi2_dof': chi2_dof_newton,
                'random_success': random_result['success'] if random_result else False,
                'random_chi2_dof': random_result['chi2_dof'] if random_result['success'] else np.inf
            }
            results.append(result)
            
            # Print summary
            fund_status = "GOOD" if chi2_dof_fund < 5 else "POOR"
            fit_status = "GOOD" if fit_result['success'] and fit_result['chi2_dof'] < 2 else "POOR"
            newton_status = "GOOD" if chi2_dof_newton < 10 else "POOR"
            random_status = "GOOD" if random_result['success'] and random_result['chi2_dof'] < 2 else "POOR"
            
            print(f"Fund-UDT: {fund_status:4s} | Fit-UDT: {fit_status:4s} | Newton: {newton_status:4s} | Random: {random_status:4s}")
        
        return results
    
    def analyze_results(self, results):
        """Analyze and interpret the comparison results."""
        
        print("\n" + "="*60)
        print("HONEST COMPARISON RESULTS")
        print("="*60)
        
        n_total = len(results)
        
        # Success rates (chi2/dof < 2 for fitted models, < 5 for predictions)
        fund_success = sum(1 for r in results if r['fundamental_udt_chi2_dof'] < 5)
        fit_success = sum(1 for r in results if r['curve_fit_success'] and r['curve_fit_chi2_dof'] < 2)
        newton_success = sum(1 for r in results if r['newtonian_chi2_dof'] < 10)
        random_success = sum(1 for r in results if r['random_success'] and r['random_chi2_dof'] < 2)
        
        print(f"Success Rates (good fits):")
        print(f"  Fundamental UDT: {fund_success}/{n_total} ({100*fund_success/n_total:.0f}%)")
        print(f"  Curve-fit UDT:   {fit_success}/{n_total} ({100*fit_success/n_total:.0f}%)")
        print(f"  Newtonian:       {newton_success}/{n_total} ({100*newton_success/n_total:.0f}%)")
        print(f"  Random model:    {random_success}/{n_total} ({100*random_success/n_total:.0f}%)")
        
        # Average chi2/dof
        fund_chi2_avg = np.mean([r['fundamental_udt_chi2_dof'] for r in results])
        fit_chi2_avg = np.mean([r['curve_fit_chi2_dof'] for r in results if r['curve_fit_success']])
        newton_chi2_avg = np.mean([r['newtonian_chi2_dof'] for r in results])
        random_chi2_avg = np.mean([r['random_chi2_dof'] for r in results if r['random_success']])
        
        print(f"\nAverage chi2/dof:")
        print(f"  Fundamental UDT: {fund_chi2_avg:.1f}")
        print(f"  Curve-fit UDT:   {fit_chi2_avg:.1f}")
        print(f"  Newtonian:       {newton_chi2_avg:.1f}")
        print(f"  Random model:    {random_chi2_avg:.1f}")
        
        print("\n" + "="*60)
        print("SCIENTIFIC INTERPRETATION")
        print("="*60)
        
        print("\n1. FUNDAMENTAL UDT PERFORMANCE:")
        if fund_success == 0:
            print("   FAILURE: Fundamental UDT makes poor predictions")
            print("   This suggests UDT is not the correct fundamental theory")
        elif fund_success < n_total // 2:
            print("   MIXED: Some success but many failures")
            print("   Theory needs refinement")
        else:
            print("   SUCCESS: Majority of predictions work")
            print("   UDT shows promise as fundamental theory")
        
        print("\n2. CURVE-FITTING EFFECTIVENESS:")
        if fit_success > fund_success * 2:
            print("   WARNING: Curve-fitting much better than fundamental theory")
            print("   This suggests 'UDT success' comes from parameter optimization, not physics")
        
        print("\n3. COMPARISON WITH BASELINES:")
        if random_success > 0:
            print(f"   DANGER: Random model succeeds {random_success}/{n_total} times")
            print("   This shows that any 2-parameter model can fit sparse data")
        
        if fit_success >= random_success:
            print("   UDT curve-fitting no better than random models")
        
        print("\n4. OVERALL ASSESSMENT:")
        if fund_success == 0 and fit_success > fund_success:
            print("   VERDICT: UDT galactic 'success' is curve-fitting artifact")
            print("   NOT evidence for fundamental physics")
        elif fund_success > n_total // 2:
            print("   VERDICT: UDT shows genuine predictive power")
            print("   May be correct fundamental theory")
        else:
            print("   VERDICT: UDT needs major theoretical development")
            print("   Current form insufficient")

def main():
    """Run the honest comparison."""
    
    comparison = HonestGalacticComparison()
    results = comparison.run_honest_comparison(n_galaxies=10)
    comparison.analyze_results(results)
    
    print("\n" + "="*60)
    print("CONCLUSIONS FOR UDT DEVELOPMENT")
    print("="*60)
    print("1. If fundamental UDT fails: Rework theoretical foundation")
    print("2. If curve-fitting succeeds: Investigate why - what physics is missing?")
    print("3. If random models work: Any 2-parameter model fits - not physics")
    print("4. Scientific integrity: Report honest results, not just successes")

if __name__ == "__main__":
    main()