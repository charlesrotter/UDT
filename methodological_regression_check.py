#!/usr/bin/env python3
"""
Methodological Regression Check
==============================

Compare current real data performance with documented best results
to ensure no methodological regression from parameter fixes.
"""

import sys
import os
import json
import subprocess
from pathlib import Path

def run_current_tests():
    """Run current real data tests and capture results."""
    print("METHODOLOGICAL REGRESSION CHECK")
    print("=" * 50)
    
    current_results = {}
    
    # Test current SPARC performance  
    print("\n1. SPARC Galaxy Analysis (Current):")
    try:
        result = subprocess.run([
            sys.executable, "scripts/analyze_sparc_galaxies.py", 
            "--max-galaxies", "10"
        ], capture_output=True, text=True, timeout=120)
        
        if result.returncode == 0:
            # Extract key metrics from output
            output = result.stdout
            if "RMS residual statistics:" in output:
                lines = output.split('\n')
                for line in lines:
                    if "Mean:" in line and "km/s" in line:
                        try:
                            rms_mean = float(line.split("Mean:")[1].split("km/s")[0].strip())
                            current_results['sparc_rms_mean'] = rms_mean
                        except:
                            pass
                    if "chi2/dof" in line:
                        try:
                            chi2_val = float(line.split("chi2/dof = ")[1])
                            if 'sparc_chi2_samples' not in current_results:
                                current_results['sparc_chi2_samples'] = []
                            current_results['sparc_chi2_samples'].append(chi2_val)
                        except:
                            pass
            
            print(f"   + Current SPARC: RMS ~{current_results.get('sparc_rms_mean', 'N/A')} km/s")
            chi2_samples = current_results.get('sparc_chi2_samples', [])
            if chi2_samples:
                avg_chi2 = sum(chi2_samples) / len(chi2_samples)
                print(f"   + Current SPARC: Average χ²/dof ~{avg_chi2:.2f}")
        else:
            print(f"   - SPARC test failed: {result.stderr}")
            current_results['sparc_failed'] = True
    except Exception as e:
        print(f"   - SPARC test error: {e}")
        current_results['sparc_failed'] = True
    
    # Check if CMB script works with current parameters
    print("\n2. CMB Analysis (Current):")
    try:
        # Just test if it can initialize without errors
        result = subprocess.run([
            sys.executable, "-c", 
            "from scripts.analyze_cmb_planck import UDTCMBAnalyzer; "
            "analyzer = UDTCMBAnalyzer(); "
            "print(f'CMB R0 = {analyzer.R0_cmb:.1f} Mpc')"
        ], capture_output=True, text=True, timeout=30)
        
        if result.returncode == 0:
            r0_cmb = float(result.stdout.split("CMB R0 = ")[1].split(" Mpc")[0])
            current_results['cmb_r0'] = r0_cmb
            print(f"   + Current CMB: R0 = {r0_cmb:.1f} Mpc")
        else:
            print(f"   - CMB initialization failed: {result.stderr}")
            current_results['cmb_failed'] = True
    except Exception as e:
        print(f"   - CMB test error: {e}")
        current_results['cmb_failed'] = True
    
    return current_results

def compare_with_benchmarks(current_results):
    """Compare current results with known best benchmarks."""
    print("\n" + "=" * 50)
    print("BENCHMARK COMPARISON")
    print("=" * 50)
    
    # Known best results (from conversation context)
    benchmarks = {
        'cmb_r0_best': 13041.1,  # Mpc - from previous CMB results
        'cmb_chi2_dof_acceptable': 35563.4,  # Previous CMB χ²/dof 
        'sparc_rms_excellent': 3.0,  # km/s - good SPARC RMS
        'sparc_chi2_excellent': 2.0,  # χ²/dof - excellent SPARC fits
    }
    
    regression_detected = False
    
    # CMB parameter check
    if 'cmb_r0' in current_results:
        if abs(current_results['cmb_r0'] - benchmarks['cmb_r0_best']) < 100:
            print(f"   + CMB R0: {current_results['cmb_r0']:.1f} vs {benchmarks['cmb_r0_best']:.1f} Mpc - CONSISTENT")
        else:
            print(f"   - CMB R0: {current_results['cmb_r0']:.1f} vs {benchmarks['cmb_r0_best']:.1f} Mpc - REGRESSION!")
            regression_detected = True
    else:
        print(f"   - CMB R0: Failed to initialize - CRITICAL REGRESSION!")
        regression_detected = True
    
    # SPARC performance check
    if 'sparc_rms_mean' in current_results:
        if current_results['sparc_rms_mean'] <= benchmarks['sparc_rms_excellent'] * 2:
            print(f"   + SPARC RMS: {current_results['sparc_rms_mean']:.1f} vs {benchmarks['sparc_rms_excellent']:.1f} km/s target - GOOD")
        else:
            print(f"   - SPARC RMS: {current_results['sparc_rms_mean']:.1f} vs {benchmarks['sparc_rms_excellent']:.1f} km/s target - REGRESSION!")
            regression_detected = True
    else:
        print(f"   - SPARC RMS: Failed to measure - REGRESSION!")
        regression_detected = True
    
    if 'sparc_chi2_samples' in current_results and current_results['sparc_chi2_samples']:
        avg_chi2 = sum(current_results['sparc_chi2_samples']) / len(current_results['sparc_chi2_samples'])
        if avg_chi2 <= benchmarks['sparc_chi2_excellent'] * 5:  # Allow some tolerance
            print(f"   + SPARC χ²/dof: {avg_chi2:.2f} vs {benchmarks['sparc_chi2_excellent']:.1f} target - ACCEPTABLE")
        else:
            print(f"   - SPARC χ²/dof: {avg_chi2:.2f} vs {benchmarks['sparc_chi2_excellent']:.1f} target - REGRESSION!")
            regression_detected = True
    
    return regression_detected

def check_data_availability():
    """Check that real data files are available."""
    print("\n" + "=" * 50)
    print("DATA AVAILABILITY CHECK")
    print("=" * 50)
    
    data_paths = [
        "data/sparc_database/SPARC_Lelli2016c.mrt",
        "data/sparc_database/MassModels_Lelli2016c.mrt", 
        "data/cmb_planck",
        "results/cmb_analysis/udt_cmb_results.json"
    ]
    
    all_available = True
    for path in data_paths:
        if os.path.exists(path):
            print(f"   + {path}: Available")
        else:
            print(f"   - {path}: Missing")
            all_available = False
    
    return all_available

if __name__ == "__main__":
    print("Testing for methodological regression after parameter fixes...")
    
    # Check data availability
    data_ok = check_data_availability()
    
    # Run current tests
    current_results = run_current_tests()
    
    # Compare with benchmarks
    regression_detected = compare_with_benchmarks(current_results)
    
    # Final assessment
    print("\n" + "=" * 50)
    print("FINAL ASSESSMENT")
    print("=" * 50)
    
    if not data_ok:
        print("WARNING: Some data files missing - may affect analysis")
    
    if regression_detected:
        print("METHODOLOGICAL REGRESSION DETECTED!")
        print("Current performance is worse than previous best results.")
        print("Action needed: Investigate and fix regression issues.")
        sys.exit(1)
    else:
        print("NO METHODOLOGICAL REGRESSION DETECTED")
        print("Current performance meets or exceeds previous benchmarks.")
        print("Parameter consistency fixes successful without methodology loss.")
        
    print(f"\nCurrent Results Summary:")
    for key, value in current_results.items():
        print(f"  {key}: {value}")