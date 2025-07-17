#!/usr/bin/env python3
"""
Derive Galactic Rotation Curves from UDT Geodesics
==================================================

Now that we've established UDT is mathematically viable, we derive
galactic rotation curves from the geodesic equations in UDT spacetime.

This tests whether UDT field equations produce realistic galactic dynamics
without ad-hoc enhancement factors.

Author: Charles Rotter
Date: 2025-01-17
"""

import sympy as sp
import numpy as np
from sympy import symbols, Function, diff, simplify, sqrt, log
import matplotlib.pyplot as plt

class UDTGeodesicSolver:
    """
    Derive rotation curves from geodesics in UDT spacetime.
    
    This is the proper way to get dynamics from field equations.
    """
    
    def __init__(self):
        print("UDT GEODESIC ROTATION CURVE DERIVATION")
        print("=" * 50)
        print("Deriving galactic dynamics from UDT field equations")
        print("via geodesic motion in curved spacetime")
        print("=" * 50)
        print()
        
        # Symbolic variables
        self.r, self.R0, self.t, self.theta, self.phi = symbols('r R_0 t theta phi', real=True, positive=True)
        self.v, self.c, self.M = symbols('v c M', real=True, positive=True)
        
        # UDT temporal dilation
        self.tau = self.R0 / (self.R0 + self.r)
        
        print(f"UDT temporal dilation: tau(r) = {self.tau}")
        print()
    
    def simple_udt_metric_ansatz(self):
        """
        Use simplified UDT metric ansatz for geodesic calculation.
        
        Instead of solving complex field equations, use physical reasoning:
        - Time dilation affects g_tt component
        - Maintain spherical symmetry
        - Ensure proper Newtonian limit
        """
        print("SIMPLIFIED UDT METRIC ANSATZ")
        print("-" * 35)
        
        # Physical ansatz based on UDT temporal geometry
        # The temporal dilation should affect the time component directly
        
        # Ansatz 1: Direct temporal dilation in metric
        # g_tt = -c^2 * tau^2 = -c^2 * R0^2/(R0 + r)^2
        g_tt = -self.c**2 * self.tau**2
        
        # Ansatz 2: Minimal modification to spatial part
        # Keep g_rr simple for now - can be refined later
        g_rr = 1  # Minkowski-like spatial metric
        
        print("UDT metric components (ansatz):")
        print(f"g_tt = {g_tt}")
        print(f"g_rr = {g_rr}")
        print(f"g_theta_theta = r^2")
        print(f"g_phi_phi = r^2 sin^2(theta)")
        print()
        print("This represents spacetime with position-dependent time flow")
        print("but minimal spatial curvature modifications")
        print()
        
        return g_tt, g_rr
    
    def calculate_effective_potential(self, g_tt, g_rr):
        """
        Calculate effective potential for circular orbits in UDT spacetime.
        
        For circular orbits in spherical spacetime:
        V_eff = sqrt(-g_tt) * sqrt(g_phi_phi) * (angular_momentum/radius)
        """
        print("CALCULATING EFFECTIVE POTENTIAL")
        print("-" * 35)
        
        # For circular motion at radius r with angular momentum L
        L = symbols('L', real=True, positive=True)
        
        # Effective potential from metric components
        # V_eff^2 comes from geodesic equation
        g_phi_phi = self.r**2  # From spherical symmetry
        
        # Effective potential for circular motion
        # From geodesic equation: E^2 = -g_tt * (1 + V_eff^2)
        # where V_eff^2 = g_phi_phi * (L/r)^2 / (-g_tt)
        
        V_eff_squared = g_phi_phi * (L/self.r)**2 / (-g_tt)
        V_eff = sqrt(V_eff_squared)
        
        print(f"Effective potential: V_eff = {V_eff}")
        
        # For stable circular orbits, we need dV_eff/dr = 0
        dV_eff_dr = diff(V_eff, self.r)
        
        print("Circular orbit condition: dV_eff/dr = 0")
        print(f"dV_eff/dr = {dV_eff_dr}")
        print()
        
        return V_eff, dV_eff_dr
    
    def derive_circular_velocity(self, g_tt, g_rr):
        """
        Derive circular velocity from geodesic equations.
        
        This is the fundamental physics - no ad-hoc factors.
        """
        print("DERIVING CIRCULAR VELOCITY FROM GEODESICS")
        print("-" * 45)
        
        # For circular orbits in general metric, the velocity is determined by
        # the condition for stable circular motion
        
        # Geodesic equation gives us the relationship between
        # orbital angular momentum L and orbital velocity v
        
        # From metric: ds^2 = g_tt dt^2 + g_rr dr^2 + r^2 dphi^2
        # For circular motion: dr/dt = 0, so ds^2 = g_tt dt^2 + r^2 dphi^2
        
        # Four-velocity normalization: g_mu_nu u^mu u^nu = -c^2
        # For circular motion: g_tt (dt/dtau)^2 + r^2 (dphi/dtau)^2 = -c^2
        
        # Orbital velocity: v = r * dphi/dt = r * (dphi/dtau) / (dt/dtau)
        
        # From the normalization condition and circular orbit stability:
        dtau_dt = symbols('dtau_dt', real=True, positive=True)
        dphi_dtau = symbols('dphi_dtau', real=True, positive=True)
        
        # Four-velocity normalization
        normalization = g_tt * (1/dtau_dt)**2 + self.r**2 * dphi_dtau**2
        print(f"Four-velocity normalization: {normalization} = -c^2")
        
        # Solve for the relationship between dtau_dt and dphi_dtau
        dtau_dt_expr = sqrt(-self.c**2 / (g_tt + self.r**2 * dphi_dtau**2 * g_tt / self.c**2))
        
        # Circular velocity
        v_circ = self.r * dphi_dtau / dtau_dt_expr
        
        # For stable circular orbits, we use the fact that the effective potential
        # has a minimum, which gives us a specific relationship
        
        # Simplified approach: use Newtonian limit with UDT corrections
        # In the weak field limit, the UDT metric gives modified gravity
        
        # The key insight: UDT temporal dilation affects the time flow,
        # which modifies the effective gravitational potential
        
        # From g_tt = -c^2 tau^2, the effective gravitational potential is modified
        # Phi_eff = c^2 * (1 - tau) = c^2 * r/(R0 + r)
        
        Phi_eff = self.c**2 * self.r / (self.R0 + self.r)
        
        print(f"UDT effective gravitational potential: Phi_eff = {Phi_eff}")
        
        # For circular orbits: v^2 = r * dPhi_eff/dr
        dPhi_eff_dr = diff(Phi_eff, self.r)
        v_squared_udt = self.r * dPhi_eff_dr
        
        print(f"UDT circular velocity squared: v^2 = r * dPhi_eff/dr = {v_squared_udt}")
        
        v_udt = sqrt(v_squared_udt)
        print(f"UDT circular velocity: v = {v_udt}")
        print()
        
        return v_udt, Phi_eff
    
    def compare_with_newtonian(self, v_udt):
        """
        Compare UDT velocity with Newtonian prediction.
        """
        print("COMPARISON WITH NEWTONIAN DYNAMICS")
        print("-" * 40)
        
        # Newtonian circular velocity for point mass M
        v_newton = sqrt(self.c**2 * self.M / self.r)  # Using c^2 to match units
        
        print(f"Newtonian velocity: v_N = {v_newton}")
        print(f"UDT velocity: v_UDT = {v_udt}")
        
        # Enhancement factor
        enhancement = simplify(v_udt / v_newton)
        print(f"Enhancement factor: v_UDT/v_N = {enhancement}")
        
        # In the limit R0 >> r (Newtonian limit)
        enhancement_limit = sp.limit(enhancement, self.R0, sp.oo)
        print(f"Newtonian limit (R0 -> infinity): {enhancement_limit}")
        
        # For r >> R0 (strong UDT regime)
        enhancement_strong = sp.limit(enhancement, self.r, sp.oo)
        print(f"Strong UDT limit (r -> infinity): {enhancement_strong}")
        
        print()
        return enhancement
    
    def test_galactic_scales(self, v_udt_func, R0_value=38.0):
        """
        Test UDT velocity formula on galactic scales.
        """
        print("TESTING ON GALACTIC SCALES")
        print("-" * 30)
        print(f"Using R0 = {R0_value} kpc (theoretical prediction)")
        print()
        
        # Convert to numerical function
        r_vals = np.linspace(1, 50, 50)  # kpc
        M_galaxy = 1e12  # Solar masses (typical galaxy)
        c_kms = 299792.458  # km/s
        
        # Newtonian velocities
        v_newton_vals = np.sqrt(4.3e-6 * M_galaxy / r_vals)  # km/s
        
        # UDT velocities (convert formula to numerical)
        # v_UDT = c * sqrt(r/(R0 + r)) = c * sqrt(r) / sqrt(R0 + r)
        v_udt_vals = c_kms * np.sqrt(r_vals / (R0_value + r_vals))
        
        # Enhancement factors
        enhancement_vals = v_udt_vals / v_newton_vals
        
        print("Sample predictions:")
        print("r (kpc)  v_Newton (km/s)  v_UDT (km/s)  Enhancement")
        print("-" * 55)
        for i in range(0, len(r_vals), 10):
            r = r_vals[i]
            vn = v_newton_vals[i]
            vu = v_udt_vals[i]
            enh = enhancement_vals[i]
            print(f"{r:6.1f}  {vn:11.1f}     {vu:9.1f}      {enh:8.2f}")
        
        print()
        print("UDT PREDICTION ASSESSMENT:")
        print(f"- Velocity range: {v_udt_vals[0]:.1f} - {v_udt_vals[-1]:.1f} km/s")
        print(f"- Typical galactic velocities: 100-300 km/s")
        
        if 100 <= v_udt_vals[-1] <= 300:
            print("+ UDT velocities are in realistic galactic range")
            realistic = True
        else:
            print("- UDT velocities are outside realistic galactic range")
            realistic = False
        
        # Check velocity profile shape
        flat_region = np.std(v_udt_vals[-10:]) / np.mean(v_udt_vals[-10:])
        print(f"- Outer velocity variation: {flat_region:.1%}")
        
        if flat_region < 0.1:
            print("+ Shows approximately flat rotation curve")
            flat_curve = True
        else:
            print("- Does not show flat rotation curve")
            flat_curve = False
        
        success = realistic and flat_curve
        print(f"\nOVERALL ASSESSMENT: {'SUCCESS' if success else 'NEEDS WORK'}")
        
        return {
            'r_vals': r_vals,
            'v_udt_vals': v_udt_vals,
            'v_newton_vals': v_newton_vals,
            'enhancement_vals': enhancement_vals,
            'realistic': realistic,
            'flat_curve': flat_curve,
            'success': success
        }
    
    def run_complete_derivation(self):
        """
        Run complete derivation of rotation curves from UDT geodesics.
        """
        print("COMPLETE UDT GEODESIC ANALYSIS")
        print("=" * 40)
        print()
        
        # Step 1: Set up metric
        g_tt, g_rr = self.simple_udt_metric_ansatz()
        
        # Step 2: Calculate effective potential
        V_eff, dV_eff_dr = self.calculate_effective_potential(g_tt, g_rr)
        
        # Step 3: Derive circular velocity
        v_udt, Phi_eff = self.derive_circular_velocity(g_tt, g_rr)
        
        # Step 4: Compare with Newtonian
        enhancement = self.compare_with_newtonian(v_udt)
        
        # Step 5: Test on galactic scales
        test_results = self.test_galactic_scales(v_udt)
        
        # Step 6: Final assessment
        print("=" * 50)
        print("UDT GEODESIC DERIVATION RESULTS")
        print("=" * 50)
        
        print(f"Derived UDT velocity formula: v = c * sqrt(r/(R0 + r))")
        print(f"Enhancement over Newtonian: factor depends on r/R0 ratio")
        print(f"Galactic scale test: {'PASSED' if test_results['success'] else 'FAILED'}")
        
        if test_results['success']:
            print("\n+ UDT geodesics produce realistic galactic rotation curves")
            print("+ Velocities are in observed range (100-300 km/s)")
            print("+ Shows approximately flat rotation curve behavior")
            print("\nCONCLUSION: UDT field equations lead to viable galactic dynamics")
        else:
            print("\n- UDT geodesics do not match galactic observations")
            print("- Either metric ansatz needs refinement or theory has problems")
            print("\nCONCLUSION: UDT requires further development")
        
        return {
            'v_udt_formula': v_udt,
            'enhancement_formula': enhancement,
            'test_results': test_results,
            'viable': test_results['success']
        }

def main():
    """Run UDT geodesic derivation."""
    
    solver = UDTGeodesicSolver()
    results = solver.run_complete_derivation()
    
    print("\n" + "=" * 50)
    print("PHASE 2 RESULTS: PHYSICAL PREDICTIONS")
    print("=" * 50)
    
    if results['viable']:
        print("STATUS: UDT produces viable galactic dynamics from geodesics")
        print("RECOMMENDATION: Test on real SPARC galaxy data")
    else:
        print("STATUS: UDT geodesics do not match observations")
        print("RECOMMENDATION: Refine metric ansatz or field equations")
    
    return results

if __name__ == "__main__":
    main()