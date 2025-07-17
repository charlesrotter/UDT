#!/usr/bin/env python3
"""
Diagnose UDT Velocity Problem
=============================

The geodesic derivation produced velocities ~50,000 km/s instead of ~200 km/s.
This is off by a factor of ~250. We need to find the source of this problem.

Possible issues:
1. Wrong metric ansatz
2. Units/scaling problem
3. Fundamental problem with UDT approach
4. Missing physics in derivation

Author: Charles Rotter
Date: 2025-01-17
"""

import sympy as sp
import numpy as np
from sympy import symbols, Function, diff, simplify, sqrt, log

class UDTVelocityDiagnostic:
    """Diagnose the source of unrealistic UDT velocities."""
    
    def __init__(self):
        print("UDT VELOCITY PROBLEM DIAGNOSTIC")
        print("=" * 40)
        print("Investigating why UDT geodesics give ~50,000 km/s")
        print("instead of realistic ~200 km/s galactic velocities")
        print("=" * 40)
        print()
        
        # Symbolic variables
        self.r, self.R0, self.c, self.M = symbols('r R_0 c M', real=True, positive=True)
        self.G = symbols('G', real=True, positive=True)
        
        # UDT temporal dilation
        self.tau = self.R0 / (self.R0 + self.r)
    
    def check_units_and_scaling(self):
        """Check if the problem is units or scaling."""
        
        print("CHECKING UNITS AND SCALING")
        print("-" * 30)
        
        # Realistic values
        R0_kpc = 38.0
        r_kpc = 10.0
        c_kms = 299792.458  # km/s
        M_solar = 1e12  # solar masses
        
        print(f"Test values:")
        print(f"R0 = {R0_kpc} kpc")
        print(f"r = {r_kpc} kpc") 
        print(f"c = {c_kms} km/s")
        print(f"M = {M_solar:.0e} solar masses")
        print()
        
        # Our UDT formula: v = c * sqrt(r/(R0 + r))
        v_udt = c_kms * np.sqrt(r_kpc / (R0_kpc + r_kpc))
        print(f"UDT velocity: v = c * sqrt(r/(R0+r)) = {v_udt:.1f} km/s")
        
        # This is the problem! We're getting c * (small factor), but c is huge
        factor = np.sqrt(r_kpc / (R0_kpc + r_kpc))
        print(f"Reduction factor: sqrt(r/(R0+r)) = {factor:.3f}")
        print(f"But c = {c_kms:.0f} km/s is enormous!")
        print()
        
        # Newtonian velocity for comparison
        v_newton = np.sqrt(4.3e-6 * M_solar / r_kpc)
        print(f"Newtonian velocity: v_N = {v_newton:.1f} km/s")
        print(f"Ratio: v_UDT/v_N = {v_udt/v_newton:.1f}")
        print()
        
        print("DIAGNOSIS: The problem is fundamental!")
        print("- UDT formula scales with speed of light c")
        print("- Galactic velocities should scale with sqrt(GM/r)")
        print("- These have completely different magnitude scales")
        print()
        
        return v_udt, v_newton
    
    def check_metric_ansatz(self):
        """Check if the metric ansatz is wrong."""
        
        print("CHECKING METRIC ANSATZ")
        print("-" * 25)
        
        print("Our ansatz was:")
        print("g_tt = -c^2 * tau^2 = -c^2 * R0^2/(R0+r)^2")
        print("g_rr = 1")
        print()
        
        print("This gives effective potential:")
        print("Phi_eff = c^2 * r/(R0+r)")
        print("And velocity: v^2 = r * dPhi_eff/dr")
        print()
        
        # Alternative ansatz - maybe we need Newtonian base
        print("ALTERNATIVE ANSATZ:")
        print("What if we start with Newtonian gravity + UDT corrections?")
        print()
        
        # Newtonian metric: g_tt = -(1 + 2*Phi/c^2) where Phi = GM/r
        print("Standard Newtonian metric:")
        print("g_tt = -(1 + 2*GM/(c^2*r))")
        print()
        
        # UDT modification - maybe tau modifies the Newtonian potential?
        print("UDT modification idea 1: tau modifies Newtonian potential")
        print("g_tt = -(1 + 2*GM*tau/(c^2*r))")
        print("This would give:")
        Phi_udt1 = self.G * self.M * self.tau / self.r
        print(f"Phi_UDT1 = {Phi_udt1}")
        
        v_squared_1 = self.r * diff(Phi_udt1, self.r)
        print(f"v^2 = r * dPhi/dr = {simplify(v_squared_1)}")
        print()
        
        print("UDT modification idea 2: tau appears as overall factor")
        print("g_tt = -tau^2 * (1 + 2*GM/(c^2*r))")
        print("This keeps Newtonian physics but modifies time flow")
        print()
        
        return Phi_udt1, v_squared_1
    
    def test_newtonian_based_ansatz(self):
        """Test metric ansatz that starts from Newtonian gravity."""
        
        print("TESTING NEWTONIAN-BASED UDT ANSATZ")
        print("-" * 40)
        
        # Ansatz: UDT modifies the Newtonian gravitational potential
        # Instead of Phi = GM/r, we get Phi_UDT = GM*tau/r = GM*R0/[r*(R0+r)]
        
        Phi_newton = self.G * self.M / self.r
        Phi_udt_modified = self.G * self.M * self.tau / self.r
        
        print(f"Newtonian potential: Phi_N = {Phi_newton}")
        print(f"UDT modified potential: Phi_UDT = {Phi_udt_modified}")
        print()
        
        # Circular velocities
        v_squared_newton = self.r * diff(Phi_newton, self.r)
        v_squared_udt = self.r * diff(Phi_udt_modified, self.r)
        
        print(f"Newtonian: v_N^2 = r * dPhi_N/dr = {v_squared_newton}")
        print(f"UDT: v_UDT^2 = r * dPhi_UDT/dr = {simplify(v_squared_udt)}")
        print()
        
        # Simplify UDT velocity
        v_udt_simplified = sqrt(simplify(v_squared_udt))
        print(f"UDT velocity: v_UDT = {v_udt_simplified}")
        
        # Enhancement factor
        enhancement = simplify(sqrt(v_squared_udt / v_squared_newton))
        print(f"Enhancement: v_UDT/v_N = {enhancement}")
        print()
        
        # Test numerically
        print("NUMERICAL TEST:")
        R0_val = 38.0  # kpc
        r_val = 10.0   # kpc  
        G_val = 4.3e-6  # (km/s)^2 * kpc / solar_mass
        M_val = 1e12    # solar masses
        
        # Newtonian
        v_n_num = np.sqrt(G_val * M_val / r_val)
        print(f"v_Newtonian = {v_n_num:.1f} km/s")
        
        # UDT (modified potential)
        tau_val = R0_val / (R0_val + r_val)
        enhancement_num = np.sqrt(tau_val * (2*R0_val + r_val) / (R0_val + r_val))
        v_udt_num = v_n_num * enhancement_num
        
        print(f"tau = {tau_val:.3f}")
        print(f"Enhancement factor = {enhancement_num:.3f}")
        print(f"v_UDT = {v_udt_num:.1f} km/s")
        print()
        
        if 100 <= v_udt_num <= 500:
            print("+ This gives realistic galactic velocities!")
            realistic = True
        else:
            print("- Still not in realistic range")
            realistic = False
        
        return v_udt_simplified, enhancement, realistic
    
    def test_alternative_approaches(self):
        """Test other possible UDT approaches."""
        
        print("TESTING ALTERNATIVE UDT APPROACHES")
        print("-" * 40)
        
        print("Approach 1: UDT as modified gravity")
        print("Standard: F = GMm/r^2")
        print("UDT: F = GMm*tau/r^2 = GMm*R0/[r^2*(R0+r)]")
        print("This gives: v^2 = GM*tau/r = GM*R0/[r*(R0+r)]")
        
        v_approach1 = sqrt(self.G * self.M * self.tau / self.r)
        print(f"v_1 = {v_approach1}")
        print()
        
        print("Approach 2: UDT as enhanced inertia")
        print("Standard: mv^2/r = GMm/r^2 -> v^2 = GM/r")
        print("UDT: m*tau*v^2/r = GMm/r^2 -> v^2 = GM/(tau*r)")
        print("This gives: v^2 = GM*(R0+r)/(R0*r)")
        
        v_approach2 = sqrt(self.G * self.M * (self.R0 + self.r) / (self.R0 * self.r))
        print(f"v_2 = {v_approach2}")
        print()
        
        print("Approach 3: UDT as effective mass enhancement")
        print("Standard: v^2 = GM/r")
        print("UDT: v^2 = GM_eff/r where M_eff = M/tau^2")
        print("This gives: v^2 = GM*(R0+r)^2/(R0^2*r)")
        
        v_approach3 = sqrt(self.G * self.M * (self.R0 + self.r)**2 / (self.R0**2 * self.r))
        print(f"v_3 = {v_approach3}")
        print()
        
        # Test numerically
        print("NUMERICAL COMPARISON:")
        R0_val = 38.0
        r_val = 10.0  
        G_val = 4.3e-6
        M_val = 1e12
        
        v_n = np.sqrt(G_val * M_val / r_val)
        
        tau_val = R0_val / (R0_val + r_val)
        v1 = np.sqrt(G_val * M_val * tau_val / r_val)
        v2 = np.sqrt(G_val * M_val * (R0_val + r_val) / (R0_val * r_val))  
        v3 = np.sqrt(G_val * M_val * (R0_val + r_val)**2 / (R0_val**2 * r_val))
        
        print(f"Newtonian: {v_n:.1f} km/s")
        print(f"Approach 1: {v1:.1f} km/s (factor {v1/v_n:.2f})")
        print(f"Approach 2: {v2:.1f} km/s (factor {v2/v_n:.2f})")  
        print(f"Approach 3: {v3:.1f} km/s (factor {v3/v_n:.2f})")
        print()
        
        # Check which are realistic
        approaches = [v1, v2, v3]
        names = ["Modified gravity", "Enhanced inertia", "Effective mass"]
        
        for i, (v, name) in enumerate(zip(approaches, names)):
            if 100 <= v <= 500:
                print(f"+ {name}: REALISTIC velocity")
            else:
                print(f"- {name}: unrealistic velocity")
        
        return approaches, names
    
    def run_complete_diagnosis(self):
        """Run complete diagnostic of UDT velocity problem."""
        
        print("COMPLETE UDT VELOCITY DIAGNOSTIC")
        print("=" * 45)
        print()
        
        # Step 1: Check units/scaling
        v_udt_original, v_newton = self.check_units_and_scaling()
        
        # Step 2: Check metric ansatz
        Phi_alt, v_alt = self.check_metric_ansatz()
        
        # Step 3: Test Newtonian-based ansatz
        v_udt_new, enhancement, realistic_new = self.test_newtonian_based_ansatz()
        
        # Step 4: Test alternative approaches
        approaches, names = self.test_alternative_approaches()
        
        # Final assessment
        print("=" * 45)
        print("DIAGNOSTIC CONCLUSIONS")
        print("=" * 45)
        
        print("PROBLEM IDENTIFIED:")
        print("1. Original ansatz g_tt = -c^2*tau^2 is fundamentally wrong")
        print("2. It scales velocities with speed of light c (~300,000 km/s)")
        print("3. Galactic velocities should scale with sqrt(GM/r) (~200 km/s)")
        print()
        
        print("SOLUTION APPROACHES:")
        if realistic_new:
            print("+ Newtonian-based UDT ansatz gives realistic velocities")
            print("  - Modify gravitational potential: Phi = GM*tau/r")
            print("  - This preserves Newtonian scaling while adding UDT effects")
        
        print()
        print("RECOMMENDATION:")
        print("- Abandon pure geodesic approach with arbitrary metric")
        print("- Use UDT as modification to Newtonian gravity")
        print("- Test the most promising approach on real galactic data")
        
        return {
            'problem_identified': True,
            'solution_found': realistic_new,
            'recommended_approach': 'Modified Newtonian gravity with UDT tau factor'
        }

def main():
    """Run UDT velocity diagnostic."""
    
    diagnostic = UDTVelocityDiagnostic()
    results = diagnostic.run_complete_diagnosis()
    
    print("\n" + "=" * 45)
    print("DIAGNOSTIC SUMMARY")
    print("=" * 45)
    
    if results['solution_found']:
        print("STATUS: Problem diagnosed and solution identified")
        print("NEXT STEP: Test modified Newtonian approach on real data")
    else:
        print("STATUS: Fundamental problems with UDT galactic dynamics")
        print("NEXT STEP: Reconsider theoretical foundations")
    
    return results

if __name__ == "__main__":
    main()