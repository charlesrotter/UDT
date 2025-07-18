#!/usr/bin/env python3
"""
Debug UDT Distance Calculation
==============================

Find why UDT distances are all coming out as 0.0 Mpc

Author: Charles Rotter
Date: 2025-07-18
"""

import numpy as np

class UDTDistanceDebug:
    def __init__(self):
        # UDT parameters
        self.R0_local = 38.0  # kpc (galactic scale)
        self.R0_cosmo = 4754.3  # Mpc (cosmological scale) 
        self.r_horizon = 27.0 * 1000 * 1000  # 27 Gly in kpc
        self.alpha = 3.0  # Power law exponent
        
        print("UDT DISTANCE DEBUG")
        print("=" * 20)
        print(f"R0_local = {self.R0_local} kpc")
        print(f"R0_cosmo = {self.R0_cosmo} Mpc")
        print(f"r_horizon = {self.r_horizon/1e6:.0f} Gly = {self.r_horizon} kpc")
        print(f"alpha = {self.alpha}")
        print()
    
    def R0_function_udt_debug(self, r_Mpc):
        """Debug version of R0(r) function."""
        print(f"  R0_function_udt_debug called with r_Mpc = {r_Mpc}")
        
        # Convert r from Mpc to kpc for horizon calculation
        r_kpc = r_Mpc * 1000
        print(f"  r_kpc = {r_kpc}")
        
        # R0(r) = R0_local × (1 + r/r_horizon)^α
        R0_local_kpc = self.R0_local  # 38 kpc
        ratio = r_kpc / self.r_horizon
        print(f"  r_kpc/r_horizon = {ratio}")
        
        scale_factor = (1 + ratio)**self.alpha
        print(f"  scale_factor = (1 + {ratio})^{self.alpha} = {scale_factor}")
        
        R0_r_kpc = R0_local_kpc * scale_factor
        print(f"  R0_r_kpc = {R0_local_kpc} * {scale_factor} = {R0_r_kpc}")
        
        R0_r_Mpc = R0_r_kpc / 1000  # Convert back to Mpc
        print(f"  R0_r_Mpc = {R0_r_Mpc}")
        
        return R0_r_Mpc
    
    def udt_distance_from_redshift_debug(self, z):
        """Debug version of UDT distance calculation."""
        print(f"\\nDEBUGGING UDT DISTANCE FOR z = {z}")
        print("-" * 40)
        
        # Initial guess using cosmological R0
        d_guess = z * self.R0_cosmo  # Mpc
        print(f"Initial guess: d_guess = {z} * {self.R0_cosmo} = {d_guess} Mpc")
        
        # Iterate to self-consistency
        for i in range(10):
            print(f"\\nIteration {i+1}:")
            print(f"  Current d_guess = {d_guess}")
            
            R0_r = self.R0_function_udt_debug(d_guess)
            d_new = z * R0_r
            print(f"  d_new = {z} * {R0_r} = {d_new}")
            
            diff = abs(d_new - d_guess)
            rel_diff = diff / d_guess if d_guess > 0 else float('inf')
            print(f"  |d_new - d_guess| = {diff}")
            print(f"  Relative difference = {rel_diff}")
            
            if rel_diff < 1e-6:
                print(f"  CONVERGED after {i+1} iterations")
                break
            
            d_guess = d_new
            print(f"  Updated d_guess = {d_guess}")
        
        print(f"\\nFINAL RESULT: d_UDT = {d_new} Mpc")
        return d_new, R0_r
    
    def test_sample_redshifts(self):
        """Test UDT distance calculation with sample redshifts."""
        print("TESTING SAMPLE REDSHIFTS")
        print("=" * 25)
        
        # Test redshifts from the analysis
        test_redshifts = [0.000785, 0.002052, 0.001803, 0.015504]
        expected_lcdm = [3.4, 8.8, 7.7, 66.4]  # Mpc
        
        for i, z in enumerate(test_redshifts):
            print(f"\\nTEST {i+1}: Galaxy with z = {z}")
            print(f"Expected LCDM distance: {expected_lcdm[i]} Mpc")
            
            d_udt, R0_r = self.udt_distance_from_redshift_debug(z)
            
            print(f"UDT distance: {d_udt} Mpc")
            print(f"Effective R0: {R0_r} Mpc")
            
            if expected_lcdm[i] > 0:
                ratio = d_udt / expected_lcdm[i]
                print(f"UDT/LCDM ratio: {ratio:.3f}")
            
            if d_udt == 0.0:
                print("ERROR: UDT distance is zero!")
                break

def main():
    debug = UDTDistanceDebug()
    debug.test_sample_redshifts()

if __name__ == "__main__":
    main()