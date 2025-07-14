#!/usr/bin/env python3
"""
Complete UDT Validation Pipeline
Run this script to perform full validation against SPARC database

Usage:
    python run_udt_validation.py
    python run_udt_validation.py --quick-test  # For subset validation
    python run_udt_validation.py --generate-plots  # Generate all plots
"""

import argparse
import sys
import os
import time
from pathlib import Path
import pandas as pd
import numpy as np

# Add src to path if running from repository root
sys.path.append('src')

try:
    from udt_core_validated import UDTFramework, SPARCAnalysis
    from sparc_data_parser import SPARCDatabaseParser, SPARCValidationSuite
    from udt_plotting_tools import create_publication_plots
except ImportError as e:
    print(f"❌ Import Error: {e}")
    print("Make sure you're running from the repository root and all files are present")
    sys.exit(1)

def quick_validation_test():
    """Quick validation test with known sample galaxies"""
    
    print("🚀 QUICK UDT VALIDATION TEST")
    print("=" * 50)
    
    # Initialize UDT framework
    udt = UDTFramework()
    sparc = SPARCAnalysis(udt)
    
    # Test galaxies from our validation
    test_galaxies = [
        {'name': 'D631-7', 'type': 10, 'distance': 7.72, 
         'luminosity': 0.196, 'vflat': 57.7, 'vflat_err': 2.7, 'quality': 1},
        {'name': 'DDO064', 'type': 10, 'distance': 6.80, 
         'luminosity': 0.157, 'vflat': 46.1, 'vflat_err': 3.9, 'quality': 1},
        {'name': 'DDO154', 'type': 10, 'distance': 4.04, 
         'luminosity': 0.053, 'vflat': 47.0, 'vflat_err': 1.0, 'quality': 2}
    ]
    
    print(f"Framework parameters:")
    print(f"  α = {udt.alpha} (geometric constant)")
    print(f"  β = {udt.beta_galactic} (galactic optimized)")
    print(f"  R₀ = {udt.R0_galactic} kpc (characteristic radius)")
    print()
    
    # Create DataFrame for validation
    test_df = pd.DataFrame(test_galaxies)
    
    # Run validation
    validator = SPARCValidationSuite(udt)
    results_df = validator.full_sample_validation(test_df, save_results=False)
    
    # Display results
    print("Individual galaxy results:")
    print("-" * 60)
    for _, row in results_df.iterrows():
        status_icon = "✅" if row['fit_quality'] == 'Excellent' else "⚠️" if row['fit_quality'] == 'Good' else "❌"
        print(f"{row['name']:8s}: {row['v_predicted']:5.1f} vs {row['v_observed']:5.1f} km/s "
              f"({row['ratio_pred_obs']*100:3.0f}%) χ²={row['chi_squared']:4.1f} {status_icon}")
    
    # Summary statistics
    total_chi2 = results_df['chi_squared'].sum()
    mean_chi2 = results_df['chi_squared'].mean()
    good_fits = (results_df['chi_squared'] < 10).sum()
    
    print(f"\nSummary:")
    print(f"  Total χ² = {total_chi2:.1f}")
    print(f"  Mean χ²  = {mean_chi2:.1f}")
    print(f"  Good fits (χ² < 10) = {good_fits}/{len(results_df)} ({good_fits/len(results_df)*100:.0f}%)")
    
    if mean_chi2 < 20:
        print("✅ QUICK TEST PASSED - Framework working correctly!")
        return True
    else:
        print("❌ QUICK TEST ISSUES - Check implementation")
        return False

def full_sparc_validation():
    """Run complete SPARC database validation"""
    
    print("🔬 FULL SPARC DATABASE VALIDATION")
    print("=" * 50)
    
    # Initialize components
    udt = UDTFramework()
    parser = SPARCDatabaseParser()
    validator = SPARCValidationSuite(udt)
    
    # Check for SPARC data
    sparc_file = Path("data/sparc_database/SPARC_Lelli2016c.mrt")
    if not sparc_file.exists():
        print(f"❌ SPARC database file not found: {sparc_file}")
        print("Please download SPARC data to data/sparc_database/ directory")
        print("Available at: http://astroweb.case.edu/SPARC/")
        return None, None
    
    print(f"✅ Found SPARC database: {sparc_file}")
    
    # Load and parse SPARC database
    try:
        start_time = time.time()
        sparc_df = parser.parse_mrt_file()
        parse_time = time.time() - start_time
        print(f"✅ Loaded {len(sparc_df)} galaxies in {parse_time:.1f}s")
    except Exception as e:
        print(f"❌ Error loading SPARC database: {e}")
        return None, None
    
    # Run validation
    print("Running validation on full sample...")
    start_time = time.time()
    results_df = validator.full_sample_validation(sparc_df)
    validation_time = time.time() - start_time
    
    print(f"✅ Validation completed in {validation_time:.1f}s")
    
    # Generate comprehensive report
    report = validator.generate_validation_report(results_df)
    
    # Display results
    print("\n📊 FULL VALIDATION RESULTS:")
    print("-" * 40)
    print(f"Total galaxies validated: {report['total_galaxies']}")
    print(f"Mean χ²: {report['mean_chi_squared']:.1f}")
    print(f"Median χ²: {report['median_chi_squared']:.1f}")
    print(f"Excellent fits (χ² < 4): {report['excellent_fit_rate']:.1f}%")
    print(f"Good+ fits (χ² < 10): {report['good_or_better_rate']:.1f}%")
    print(f"Mean velocity ratio: {report['mean_velocity_ratio']:.2f} ± {report['std_velocity_ratio']:.2f}")
    
    # Theory assessment
    if report['good_or_better_rate'] > 70:
        print("\n🎉 VALIDATION HIGHLY SUCCESSFUL!")
        print("UDT theory strongly supported by SPARC data")
        status = "excellent"
    elif report['good_or_better_rate'] > 50:
        print("\n✅ VALIDATION SUCCESSFUL")
        print("UDT shows good agreement with observations")
        status = "good"
    elif report['good_or_better_rate'] > 30:
        print("\n⚠️ VALIDATION PROMISING")
        print("UDT shows promise but needs refinement")
        status = "promising"
    else:
        print("\n❌ VALIDATION NEEDS WORK")
        print("Theory requires significant improvement")
        status = "poor"
    
    return results_df, report

def generate_all_plots(results_df, report):
    """Generate publication-quality plots"""
    
    print("\n🎨 GENERATING PUBLICATION PLOTS")
    print("=" * 40)
    
    output_dir = "results/figures/"
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    try:
        create_publication_plots(results_df, report, output_dir)
        print("✅ All plots generated successfully")
        print(f"📁 Plots saved to: {output_dir}")
        
        # List generated files
        plot_files = list(Path(output_dir).glob("*.png"))
        print("\nGenerated plots:")
        for plot_file in plot_files:
            print(f"  📊 {plot_file.name}")
            
    except Exception as e:
        print(f"❌ Error generating plots: {e}")
        return False
    
    return True

def create_validation_summary(results_df, report, output_file="results/validation_summary.md"):
    """Create markdown summary of validation results"""
    
    os.makedirs("results", exist_ok=True)
    
    with open(output_file, 'w') as f:
        f.write("# UDT Validation Summary\n\n")
        f.write(f"**Generated:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## Theory Parameters\n\n")
        f.write(f"- α = {report['theoretical_parameters']['alpha']} (geometric constant)\n")
        f.write(f"- β = {report['theoretical_parameters']['beta_galactic']} (galactic scale)\n")
        f.write(f"- R₀ = {report['theoretical_parameters']['R0_galactic_kpc']} kpc (characteristic radius)\n\n")
        
        f.write("## Validation Results\n\n")
        f.write(f"- **Total galaxies:** {report['total_galaxies']}\n")
        f.write(f"- **Mean χ²:** {report['mean_chi_squared']:.1f}\n")
        f.write(f"- **Excellent fits:** {report['excellent_fit_rate']:.1f}%\n")
        f.write(f"- **Good+ fits:** {report['good_or_better_rate']:.1f}%\n")
        f.write(f"- **Mean ratio:** {report['mean_velocity_ratio']:.2f} ± {report['std_velocity_ratio']:.2f}\n\n")
        
        f.write("## Best Performing Galaxies\n\n")
        best_galaxies = results_df.nsmallest(5, 'chi_squared')
        f.write("| Galaxy | V_obs (km/s) | V_pred (km/s) | Ratio | χ² |\n")
        f.write("|--------|--------------|---------------|-------|----|\n")
        for _, galaxy in best_galaxies.iterrows():
            f.write(f"| {galaxy['name']} | {galaxy['v_observed']:.1f} | "
                   f"{galaxy['v_predicted']:.1f} | {galaxy['ratio_pred_obs']:.2f} | {galaxy['chi_squared']:.1f} |\n")
        
        f.write("\n## Status\n\n")
        if report['good_or_better_rate'] > 70:
            f.write("✅ **VALIDATION SUCCESSFUL** - Theory strongly supported by data\n")
        elif report['good_or_better_rate'] > 50:
            f.write("✅ **VALIDATION GOOD** - Theory shows good agreement\n")
        else:
            f.write("⚠️ **VALIDATION PARTIAL** - Theory needs refinement\n")
    
    print(f"📄 Validation summary saved to: {output_file}")

def main():
    """Main validation pipeline"""
    
    parser = argparse.ArgumentParser(description='UDT Validation Pipeline')
    parser.add_argument('--quick-test', action='store_true', 
                       help='Run quick validation test only')
    parser.add_argument('--generate-plots', action='store_true',
                       help='Generate publication plots')
    parser.add_argument('--skip-validation', action='store_true',
                       help='Skip validation, load existing results')
    
    args = parser.parse_args()
    
    print("🔬 Universal Distance Dilation Theory - Validation Pipeline")
    print("=" * 60)
    print(f"Start time: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Quick test mode
    if args.quick_test:
        success = quick_validation_test()
        if not success:
            print("\n❌ Quick test failed - check implementation")
            sys.exit(1)
        print("\n🎯 Quick test completed successfully!")
        return
    
    # Load existing results if requested
    if args.skip_validation:
        results_file = "results/sparc_validation_results.csv"
        if Path(results_file).exists():
            print(f"📂 Loading existing results from {results_file}")
            results_df = pd.read_csv(results_file)
            udt = UDTFramework()
            validator = SPARCValidationSuite(udt)
            report = validator.generate_validation_report(results_df)
        else:
            print(f"❌ Results file not found: {results_file}")
            print("Run without --skip-validation to generate results")
            sys.exit(1)
    else:
        # Run full validation
        results_df, report = full_sparc_validation()
        
        if results_df is None:
            print("❌ Validation failed - check SPARC data availability")
            sys.exit(1)
    
    # Generate plots if requested or if validation was successful
    if args.generate_plots or (report and report['good_or_better_rate'] > 30):
        print("\n" + "="*50)
        plot_success = generate_all_plots(results_df, report)
        
        if plot_success:
            print("✅ Plot generation completed")
        else:
            print("❌ Plot generation failed")
    
    # Create summary
    create_validation_summary(results_df, report)
    
    # Final status
    print("\n" + "="*60)
    print("🏁 VALIDATION PIPELINE COMPLETED")
    
    if report['good_or_better_rate'] > 70:
        print("🎉 RESULT: UDT theory strongly validated against SPARC data!")
    elif report['good_or_better_rate'] > 50:
        print("✅ RESULT: UDT theory shows good agreement with data")
    else:
        print("⚠️ RESULT: UDT theory shows promise but needs refinement")
    
    print(f"📁 Results available in: results/")
    print(f"End time: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")

if __name__ == "__main__":
    main()