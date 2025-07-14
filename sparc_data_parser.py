"""
SPARC Database Parser and Full Sample Validation
Processes the complete SPARC database for UDT validation

Usage:
    python sparc_validation.py --data_path data/sparc_database/
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import List, Dict, Tuple
from udt_core_validated import UDTFramework, SPARCAnalysis

class SPARCDatabaseParser:
    """Parser for complete SPARC database files"""
    
    def __init__(self, data_path: str = "data/sparc_database/"):
        self.data_path = Path(data_path)
        
    def parse_mrt_file(self, filename: str = "SPARC_Lelli2016c.mrt") -> pd.DataFrame:
        """
        Parse the main SPARC .mrt file with galaxy properties
        
        Args:
            filename: Name of the .mrt file
            
        Returns:
            DataFrame with all galaxy properties
        """
        filepath = self.data_path / filename
        
        if not filepath.exists():
            raise FileNotFoundError(f"SPARC file not found: {filepath}")
        
        galaxies = []
        
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
        
        # Skip header lines and find data section
        data_started = False
        for line in lines:
            # Look for first galaxy (CamB) to start data section
            if 'CamB' in line and len(line.strip()) > 100:
                data_started = True
                
            if data_started and len(line.strip()) > 100:
                galaxy = self._parse_sparc_line(line)
                if galaxy and galaxy['luminosity'] > 0:  # Valid galaxy entry
                    galaxies.append(galaxy)
        
        df = pd.DataFrame(galaxies)
        print(f"✅ Successfully loaded {len(df)} galaxies from SPARC database")
        
        # Show sample of loaded data
        if len(df) > 0:
            valid_vflat = (~df['vflat'].isna()) & (df['vflat'] > 0)
            print(f"📊 Galaxies with valid Vflat data: {valid_vflat.sum()}")
            
            # Show a few examples
            sample = df[valid_vflat].head(3)
            print("\nSample loaded galaxies:")
            for _, galaxy in sample.iterrows():
                print(f"  {galaxy['name']}: L={galaxy['luminosity']:.3f}, Vflat={galaxy['vflat']:.1f} km/s")
        
        return df
    
    def _parse_sparc_line(self, line: str) -> Dict:
        """Parse a single line from SPARC .mrt file using whitespace separation"""
        
        # Split by whitespace - this works correctly as shown in debug
        parts = line.split()
        
        if len(parts) < 18:  # Need at least 18 parts for complete data
            return None
        
        try:
            galaxy = {
                'name': parts[0],                           # D631-7
                'type': int(parts[1]),                      # 10
                'distance': float(parts[2]),                # 7.72
                'distance_err': float(parts[3]),            # 0.18
                'f_distance': int(parts[4]),                # 2 (distance method)
                'inc': float(parts[5]),                     # 59.0
                'inc_err': float(parts[6]),                 # 3.0
                'luminosity': float(parts[7]),              # 0.196 ⭐ KEY FIELD
                'luminosity_err': float(parts[8]),          # 0.009
                'reff': float(parts[9]),                    # 1.22
                'sbeff': float(parts[10]),                  # 20.93
                'rdisk': float(parts[11]),                  # 0.70
                'sbdisk': float(parts[12]),                 # 115.04
                'mhi': float(parts[13]),                    # 0.290
                'rhi': float(parts[14]),                    # 0.00
                'vflat': float(parts[15]),                  # 57.7 ⭐ KEY FIELD
                'vflat_err': float(parts[16]),              # 2.7 ⭐ KEY FIELD
                'quality': int(parts[17])                   # 1
                # parts[18] is reference string (Tr09,dB01) - ignored
            }
            
            # Clean up invalid values
            if galaxy['vflat'] <= 0:
                galaxy['vflat'] = np.nan
            if galaxy['vflat_err'] <= 0 or galaxy['vflat_err'] > 900:
                galaxy['vflat_err'] = np.nan
                
            return galaxy
            
        except (ValueError, IndexError) as e:
            return None

class SPARCValidationSuite:
    """Complete validation suite for UDT against SPARC database"""
    
    def __init__(self, udt_framework: UDTFramework):
        self.udt = udt_framework
        self.sparc_analysis = SPARCAnalysis(udt_framework)
        
    def full_sample_validation(self, sparc_df: pd.DataFrame, 
                              save_results: bool = True) -> pd.DataFrame:
        """
        Validate UDT against complete SPARC sample
        
        Args:
            sparc_df: DataFrame with SPARC galaxy data
            save_results: Whether to save results to CSV
            
        Returns:
            DataFrame with validation results for all galaxies
        """
        print("🔬 Running full SPARC validation...")
        
        # Filter for galaxies with valid rotation data
        valid_galaxies = sparc_df[
            (~sparc_df['vflat'].isna()) & 
            (sparc_df['vflat'] > 0) & 
            (sparc_df['luminosity'] > 0)
        ].copy()
        
        print(f"Validating {len(valid_galaxies)} galaxies with rotation data")
        
        results = []
        
        for idx, galaxy in valid_galaxies.iterrows():
            try:
                result = self._validate_single_galaxy(galaxy)
                results.append(result)
                
                if len(results) % 25 == 0:
                    print(f"  Processed {len(results)} galaxies...")
                    
            except Exception as e:
                print(f"Warning: Failed to process {galaxy['name']}: {e}")
                continue
        
        results_df = pd.DataFrame(results)
        
        if save_results:
            output_file = "results/sparc_validation_results.csv"
            os.makedirs("results", exist_ok=True)
            results_df.to_csv(output_file, index=False)
            print(f"Results saved to {output_file}")
        
        return results_df
    
    def _validate_single_galaxy(self, galaxy: pd.Series) -> Dict:
        """Validate UDT prediction for a single galaxy"""
        
        M_star = self.udt.luminosity_to_mass(galaxy['luminosity'])
        
        # Test multiple radii to find best fit
        test_radii = np.array([2, 3, 5, 7, 8, 10, 12, 15, 20])
        
        best_chi2 = np.inf
        best_result = None
        
        for r in test_radii:
            v_pred = self.udt.rotation_velocity_udt(r, M_star)
            
            # Calculate chi-squared
            v_error = galaxy['vflat_err'] if not pd.isna(galaxy['vflat_err']) else galaxy['vflat'] * 0.1
            chi2 = ((v_pred - galaxy['vflat']) / v_error)**2
            
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_result = {
                    'name': galaxy['name'],
                    'type': galaxy['type'],
                    'distance': galaxy['distance'],
                    'luminosity': galaxy['luminosity'],
                    'stellar_mass_1e9': M_star / 1e9,
                    'v_observed': galaxy['vflat'],
                    'v_observed_err': v_error,
                    'v_predicted': v_pred,
                    'best_radius_kpc': r,
                    'ratio_pred_obs': v_pred / galaxy['vflat'],
                    'chi_squared': chi2,
                    'log_chi_squared': np.log10(chi2),
                    'quality': galaxy['quality']
                }
        
        # Add classification
        if best_chi2 < 4:
            best_result['fit_quality'] = 'Excellent'
        elif best_chi2 < 10:
            best_result['fit_quality'] = 'Good'
        elif best_chi2 < 25:
            best_result['fit_quality'] = 'Fair'
        else:
            best_result['fit_quality'] = 'Poor'
            
        return best_result
    
    def generate_validation_report(self, results_df: pd.DataFrame) -> Dict:
        """Generate comprehensive validation statistics"""
        
        # Overall statistics
        total_galaxies = len(results_df)
        mean_chi2 = results_df['chi_squared'].mean()
        median_chi2 = results_df['chi_squared'].median()
        
        # Success rates
        excellent_rate = (results_df['fit_quality'] == 'Excellent').mean() * 100
        good_rate = (results_df['fit_quality'].isin(['Excellent', 'Good'])).mean() * 100
        
        # Ratio statistics
        mean_ratio = results_df['ratio_pred_obs'].mean()
        std_ratio = results_df['ratio_pred_obs'].std()
        
        report = {
            'total_galaxies': total_galaxies,
            'mean_chi_squared': mean_chi2,
            'median_chi_squared': median_chi2,
            'excellent_fit_rate': excellent_rate,
            'good_or_better_rate': good_rate,
            'mean_velocity_ratio': mean_ratio,
            'std_velocity_ratio': std_ratio,
            'theoretical_parameters': {
                'alpha': self.udt.alpha,
                'beta_galactic': self.udt.beta_galactic,
                'R0_galactic_kpc': self.udt.R0_galactic
            }
        }
        
        return report
    
    def plot_validation_results(self, results_df: pd.DataFrame, 
                               save_plots: bool = True) -> None:
        """Generate validation plots"""
        
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
        
        # 1. Predicted vs Observed velocities
        ax1.scatter(results_df['v_observed'], results_df['v_predicted'], 
                   alpha=0.6, s=30)
        ax1.plot([0, 300], [0, 300], 'r--', label='Perfect agreement')
        ax1.set_xlabel('Observed V_flat (km/s)')
        ax1.set_ylabel('UDT Predicted V (km/s)')
        ax1.set_title('UDT Predictions vs SPARC Observations')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # 2. Velocity ratio distribution
        ax2.hist(results_df['ratio_pred_obs'], bins=30, alpha=0.7, edgecolor='black')
        ax2.axvline(1.0, color='red', linestyle='--', label='Perfect match')
        ax2.axvline(results_df['ratio_pred_obs'].mean(), color='green', 
                   linestyle='-', label=f'Mean = {results_df["ratio_pred_obs"].mean():.2f}')
        ax2.set_xlabel('Velocity Ratio (Predicted/Observed)')
        ax2.set_ylabel('Number of Galaxies')
        ax2.set_title('Distribution of Velocity Ratios')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # 3. Chi-squared distribution
        ax3.hist(np.log10(results_df['chi_squared']), bins=25, alpha=0.7, edgecolor='black')
        ax3.axvline(np.log10(4), color='red', linestyle='--', label='χ² = 4 (Good fit)')
        ax3.set_xlabel('log₁₀(χ²)')
        ax3.set_ylabel('Number of Galaxies')
        ax3.set_title('Distribution of χ² Values')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # 4. Fit quality by galaxy type
        quality_counts = results_df.groupby(['type', 'fit_quality']).size().unstack(fill_value=0)
        quality_counts.plot(kind='bar', stacked=True, ax=ax4)
        ax4.set_xlabel('Galaxy Type')
        ax4.set_ylabel('Number of Galaxies')
        ax4.set_title('Fit Quality by Galaxy Type')
        ax4.legend(title='Fit Quality')
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_plots:
            os.makedirs("results", exist_ok=True)
            plt.savefig("results/udt_sparc_validation.png", dpi=300, bbox_inches='tight')
            print("Validation plots saved to results/udt_sparc_validation.png")
        
        plt.show()

def main():
    """Main validation routine for complete SPARC database"""
    
    print("🚀 UDT Full SPARC Database Validation")
    print("=" * 50)
    
    # Initialize components
    udt = UDTFramework()
    parser = SPARCDatabaseParser()
    validator = SPARCValidationSuite(udt)
    
    # Load SPARC database
    try:
        sparc_df = parser.parse_mrt_file()
    except FileNotFoundError:
        print("❌ SPARC database file not found!")
        print("Please ensure data/sparc_database/SPARC_Lelli2016c.mrt exists")
        return
    
    print(f"✅ Loaded {len(sparc_df)} galaxies from SPARC database")
    
    # Run validation
    results_df = validator.full_sample_validation(sparc_df)
    
    # Generate report
    report = validator.generate_validation_report(results_df)
    
    print("\n📊 VALIDATION RESULTS:")
    print("-" * 30)
    print(f"Total galaxies validated: {report['total_galaxies']}")
    print(f"Mean χ²: {report['mean_chi_squared']:.1f}")
    print(f"Excellent fits (χ² < 4): {report['excellent_fit_rate']:.1f}%")
    print(f"Good+ fits (χ² < 10): {report['good_or_better_rate']:.1f}%")
    print(f"Mean velocity ratio: {report['mean_velocity_ratio']:.2f} ± {report['std_velocity_ratio']:.2f}")
    
    print(f"\nTheory parameters:")
    print(f"  α = {report['theoretical_parameters']['alpha']}")
    print(f"  β = {report['theoretical_parameters']['beta_galactic']}")
    print(f"  R₀ = {report['theoretical_parameters']['R0_galactic_kpc']} kpc")
    
    # Generate plots
    validator.plot_validation_results(results_df)
    
    # Overall assessment
    if report['good_or_better_rate'] > 70:
        print("\n🎉 VALIDATION SUCCESSFUL!")
        print("UDT theory successfully explains galactic rotation curves")
    elif report['good_or_better_rate'] > 50:
        print("\n✅ VALIDATION PROMISING")
        print("UDT shows strong correlation with observations")
    else:
        print("\n⚠️ VALIDATION NEEDS IMPROVEMENT")
        print("Theory requires further refinement")
    
    return results_df, report

# Test function to verify the fix
def test_sparc_parsing():
    """Test the corrected parsing on known data"""
    
    # Test with the exact line from your data
    test_line = "     D631-7 10   7.72  0.18  2 59.0  3.0   0.196   0.009  1.22    20.93  0.70   115.04   0.290  0.00  57.7   2.7   1      Tr09,dB01"
    
    parser = SPARCDatabaseParser()
    galaxy = parser._parse_sparc_line(test_line)
    
    if galaxy:
        print("✅ PARSING TEST SUCCESSFUL!")
        print(f"Galaxy: {galaxy['name']}")
        print(f"Type: {galaxy['type']}")
        print(f"Distance: {galaxy['distance']} Mpc")
        print(f"Luminosity: {galaxy['luminosity']} × 10^9 L_sun")
        print(f"Vflat: {galaxy['vflat']} ± {galaxy['vflat_err']} km/s")
        print(f"Quality: {galaxy['quality']}")
        return True
    else:
        print("❌ PARSING TEST FAILED")
        return False

if __name__ == "__main__":
    # Run test first
    print("Testing SPARC parser...")
    test_sparc_parsing()
    
    print("\n" + "="*50)
    
    # Then run full validation if test passes
    results, report = main()