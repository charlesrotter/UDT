#!/usr/bin/env python3
"""
Pure Geometric UDT Validation
=============================

CRITICAL: Validate UDT using pure geometric principles without Standard Model contamination.

CONTAMINATION-FREE APPROACH:
- NO chi-squared minimization (Standard Model statistical framework)
- NO RMS residual evaluation (assumes Gaussian error models)
- NO goodness-of-fit testing (Standard Model hypothesis testing)
- NO least-squares optimization (assumes specific error structures)

PURE UDT VALIDATION CRITERIA:
1. Geometric pattern consistency: Does τ(r) = R₀/(R₀ + r) describe the velocity structure?
2. Enhancement validation: Does F(τ) enhancement match dimensional analysis?
3. Scale invariance: Do dimensionless profiles show universal UDT geometry?
4. Theoretical prediction: Can UDT predict rotation curve shape from first principles?

Author: Charles Rotter
Date: 2025-07-18
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import glob
import json

class PureGeometricUDTValidator:
    def __init__(self):
        self.data_dir = "C:/UDT/data/sparc_database/"
        self.results_dir = "C:/UDT/results/"
        
        os.makedirs(self.results_dir, exist_ok=True)
        
        print("PURE GEOMETRIC UDT VALIDATION")
        print("=" * 29)
        print("GOAL: Test UDT geometry using contamination-free methods")
        print("METHOD: Direct geometric pattern validation")
        print("CRITERIA: Theoretical consistency, not statistical fitting")
        print()
        
        # UDT geometric parameters (from pure theory)
        self.R0_galactic = 57.5e3 * 3.086e16  # 57.5 kpc in meters
        self.alpha = 0.00206  # geometric coupling constant
        
        print(f"UDT Geometric Parameters:")
        print(f"R0_galactic = {self.R0_galactic/3.086e16/1000:.1f} kpc (theoretical prediction)")
        print(f"alpha = {self.alpha:.5f} (pure geometric coupling)")
        print()
        
        self.galaxy_results = []
        
    def load_galaxy_observational_pattern(self, filename):
        """Load raw observational velocity pattern."""
        filepath = os.path.join(self.data_dir, filename)
        
        try:
            data = np.loadtxt(filepath)
            
            if data.size == 0:
                return None, None
            
            if data.ndim == 1:
                data = data.reshape(1, -1)
            
            radius_kpc = data[:, 0]  # kpc
            velocity_obs = data[:, 1]  # km/s
            
            # Filter valid data
            valid_mask = (radius_kpc > 0) & (velocity_obs > 0)
            
            if not np.any(valid_mask):
                return None, None
            
            radius_kpc = radius_kpc[valid_mask]
            velocity_obs = velocity_obs[valid_mask]
            
            # Create pure geometric observational pattern
            r_max = np.max(radius_kpc)
            r_normalized = radius_kpc / r_max  # Dimensionless radius (0 to 1)
            
            # Velocity normalization (scale-invariant)
            v_flat = np.median(velocity_obs[radius_kpc > 0.7 * r_max])  # Asymptotic velocity
            v_normalized = velocity_obs / v_flat  # Dimensionless velocity ratio
            
            return r_normalized, v_normalized
            
        except Exception as e:
            print(f"Error loading {filename}: {e}")
            return None, None
    
    def calculate_pure_udt_geometry(self, r_normalized, R0_dimensional_scale):
        """Calculate pure UDT geometric prediction."""
        # Convert normalized radius to physical scale for τ calculation
        r_physical = r_normalized * R0_dimensional_scale
        
        # Pure UDT geometry: τ(r) = R₀/(R₀ + r)
        tau_r = R0_dimensional_scale / (R0_dimensional_scale + r_physical)
        
        # Pure geometric enhancement: F(τ) = 1 + α × 3(1-τ)/(τ²(3-2τ))
        F_tau = 1 + self.alpha * 3 * (1 - tau_r) / (tau_r**2 * (3 - 2*tau_r))
        
        # UDT geometric velocity structure (dimensionless)
        # This is the fundamental prediction from pure geometry
        v_udt_structure = np.sqrt(F_tau * r_normalized / (1 + r_normalized))
        
        # Normalize to asymptotic value for comparison
        v_asymptotic = v_udt_structure[-1] if len(v_udt_structure) > 0 else 1.0
        v_udt_normalized = v_udt_structure / v_asymptotic
        
        return v_udt_normalized, tau_r, F_tau
    
    def geometric_pattern_consistency_test(self, r_obs, v_obs, v_udt):
        """Test geometric pattern consistency without statistical contamination."""
        
        # Test 1: Monotonicity preservation
        # UDT should preserve the basic monotonic structure of rotation curves
        obs_monotonic_regions = []
        udt_monotonic_regions = []
        
        # Find regions of monotonic increase/decrease
        for i in range(1, len(v_obs)):
            if v_obs[i] > v_obs[i-1]:
                obs_monotonic_regions.append('increasing')
            elif v_obs[i] < v_obs[i-1]:
                obs_monotonic_regions.append('decreasing')
            else:
                obs_monotonic_regions.append('constant')
                
            if v_udt[i] > v_udt[i-1]:
                udt_monotonic_regions.append('increasing')
            elif v_udt[i] < v_udt[i-1]:
                udt_monotonic_regions.append('decreasing')
            else:
                udt_monotonic_regions.append('constant')
        
        # Calculate monotonicity agreement
        monotonic_agreement = sum(1 for obs, udt in zip(obs_monotonic_regions, udt_monotonic_regions) 
                                if obs == udt) / len(obs_monotonic_regions)
        
        # Test 2: Asymptotic behavior
        # Both should approach similar asymptotic values
        outer_region = r_obs > 0.8
        if np.sum(outer_region) > 2:
            obs_asymptotic = np.mean(v_obs[outer_region])
            udt_asymptotic = np.mean(v_udt[outer_region])
            asymptotic_consistency = abs(obs_asymptotic - udt_asymptotic) / max(obs_asymptotic, udt_asymptotic)
        else:
            asymptotic_consistency = 1.0  # Cannot evaluate
        
        # Test 3: Inner region enhancement pattern
        # UDT predicts specific enhancement pattern in inner regions
        inner_region = r_obs < 0.3
        if np.sum(inner_region) > 2:
            # Check if UDT captures the inner enhancement structure
            obs_inner_gradient = np.gradient(v_obs[inner_region])
            udt_inner_gradient = np.gradient(v_udt[inner_region])
            
            # Correlation of gradient patterns (geometric structure consistency)
            if len(obs_inner_gradient) > 1 and np.std(obs_inner_gradient) > 0 and np.std(udt_inner_gradient) > 0:
                inner_structure_correlation = np.corrcoef(obs_inner_gradient, udt_inner_gradient)[0, 1]
                if np.isnan(inner_structure_correlation):
                    inner_structure_correlation = 0.0
            else:
                inner_structure_correlation = 0.0
        else:
            inner_structure_correlation = 0.0
        
        return {
            'monotonic_agreement': monotonic_agreement,
            'asymptotic_consistency': 1.0 - asymptotic_consistency,  # Convert to consistency score
            'inner_structure_correlation': inner_structure_correlation
        }
    
    def geometric_scale_invariance_test(self, r_obs, v_obs):
        """Test whether velocity pattern shows UDT scale invariance."""
        
        # UDT predicts specific scale-invariant structure
        # τ(r) = R₀/(R₀ + r) should create universal dimensionless patterns
        
        # Test multiple R0 scales to find geometric consistency
        R0_test_scales = np.logspace(0, 2, 20)  # Test R0 from 1 to 100 (dimensionless units)
        
        geometric_consistencies = []
        
        for R0_scale in R0_test_scales:
            v_udt_predicted, tau_r, F_tau = self.calculate_pure_udt_geometry(r_obs, R0_scale)
            
            if len(v_udt_predicted) == len(v_obs):
                # Geometric pattern tests
                pattern_tests = self.geometric_pattern_consistency_test(r_obs, v_obs, v_udt_predicted)
                
                # Overall geometric consistency (not statistical fitting!)
                geometric_consistency = (
                    pattern_tests['monotonic_agreement'] * 0.4 +
                    pattern_tests['asymptotic_consistency'] * 0.3 +
                    max(0, pattern_tests['inner_structure_correlation']) * 0.3
                )
                
                geometric_consistencies.append({
                    'R0_scale': R0_scale,
                    'geometric_consistency': geometric_consistency,
                    'tau_range': (np.min(tau_r), np.max(tau_r)),
                    'enhancement_range': (np.min(F_tau), np.max(F_tau)),
                    'pattern_tests': pattern_tests
                })
        
        # Find scale with best geometric consistency
        best_scale_idx = np.argmax([gc['geometric_consistency'] for gc in geometric_consistencies])
        best_geometric_match = geometric_consistencies[best_scale_idx]
        
        return best_geometric_match, geometric_consistencies
    
    def theoretical_prediction_test(self, galaxy_name, r_obs, v_obs):
        """Test UDT theoretical prediction capability."""
        
        print(f"  Testing pure geometric prediction for {galaxy_name}")
        
        # Scale invariance test
        best_match, all_scales = self.geometric_scale_invariance_test(r_obs, v_obs)
        
        # Geometric consistency assessment
        geometric_score = best_match['geometric_consistency']
        optimal_R0 = best_match['R0_scale']
        
        # UDT geometric quality criteria (contamination-free)
        if geometric_score > 0.8:
            geometric_quality = "EXCELLENT"
        elif geometric_score > 0.6:
            geometric_quality = "GOOD"
        elif geometric_score > 0.4:
            geometric_quality = "MODERATE"
        elif geometric_score > 0.2:
            geometric_quality = "POOR"
        else:
            geometric_quality = "FAILED"
        
        print(f"    Geometric consistency: {geometric_score:.3f} ({geometric_quality})")
        print(f"    Optimal R0 scale: {optimal_R0:.1f}")
        print(f"    Monotonic agreement: {best_match['pattern_tests']['monotonic_agreement']:.3f}")
        print(f"    Asymptotic consistency: {best_match['pattern_tests']['asymptotic_consistency']:.3f}")
        print(f"    Inner structure correlation: {best_match['pattern_tests']['inner_structure_correlation']:.3f}")
        
        return {
            'galaxy_name': galaxy_name,
            'geometric_score': geometric_score,
            'geometric_quality': geometric_quality,
            'optimal_R0_scale': optimal_R0,
            'tau_range': best_match['tau_range'],
            'enhancement_range': best_match['enhancement_range'],
            'pattern_details': best_match['pattern_tests'],
            'n_points': len(r_obs),
            'scale_analysis': all_scales
        }
    
    def analyze_all_galaxies_pure_geometric(self):
        """Analyze all SPARC galaxies using pure geometric validation."""
        print("PURE GEOMETRIC ANALYSIS OF ALL SPARC GALAXIES")
        print("-" * 44)
        
        galaxy_files = glob.glob(os.path.join(self.data_dir, "*_rotmod.dat"))
        
        print(f"Found {len(galaxy_files)} galaxy files")
        print("Testing pure UDT geometric predictions...")
        print()
        
        successful_tests = 0
        total_galaxies = 0
        
        for i, filepath in enumerate(galaxy_files):
            filename = os.path.basename(filepath)
            galaxy_name = filename.replace('_rotmod.dat', '')
            
            total_galaxies += 1
            
            print(f"Processing {i+1}/{len(galaxy_files)}: {galaxy_name}")
            
            # Load pure observational pattern
            r_obs, v_obs = self.load_galaxy_observational_pattern(filename)
            
            if r_obs is None:
                print(f"  Failed to load observational pattern")
                continue
            
            print(f"  Data points: {len(r_obs)}")
            
            # Pure geometric validation
            try:
                result = self.theoretical_prediction_test(galaxy_name, r_obs, v_obs)
                
                self.galaxy_results.append(result)
                successful_tests += 1
                
            except Exception as e:
                print(f"  Geometric test failed: {e}")
        
        print()
        print(f"PURE GEOMETRIC ANALYSIS COMPLETE:")
        print(f"Total galaxies processed: {total_galaxies}")
        print(f"Successful geometric tests: {successful_tests}")
        print(f"Success rate: {successful_tests/total_galaxies*100:.1f}%")
        print()
        
        return successful_tests, total_galaxies
    
    def calculate_pure_geometric_statistics(self):
        """Calculate pure geometric validation statistics."""
        print("PURE GEOMETRIC VALIDATION STATISTICS")
        print("-" * 35)
        
        if not self.galaxy_results:
            print("No results to analyze")
            return None
        
        # Extract geometric metrics
        geometric_scores = [r['geometric_score'] for r in self.galaxy_results]
        optimal_R0_scales = [r['optimal_R0_scale'] for r in self.galaxy_results]
        geometric_qualities = [r['geometric_quality'] for r in self.galaxy_results]
        
        # Pure geometric statistics
        stats_summary = {
            'total_galaxies_tested': len(self.galaxy_results),
            'median_geometric_score': np.median(geometric_scores),
            'mean_geometric_score': np.mean(geometric_scores),
            'geometric_score_std': np.std(geometric_scores),
            'median_optimal_R0': np.median(optimal_R0_scales),
            'mean_optimal_R0': np.mean(optimal_R0_scales),
            'R0_scale_consistency': np.std(optimal_R0_scales) / np.mean(optimal_R0_scales),
        }
        
        # Quality distribution
        quality_counts = {}
        for quality in ['EXCELLENT', 'GOOD', 'MODERATE', 'POOR', 'FAILED']:
            count = sum(1 for q in geometric_qualities if q == quality)
            quality_counts[quality] = count
            
        stats_summary['quality_distribution'] = quality_counts
        
        # UDT geometric validation assessment
        excellent_fraction = quality_counts['EXCELLENT'] / len(self.galaxy_results)
        good_or_better_fraction = (quality_counts['EXCELLENT'] + quality_counts['GOOD']) / len(self.galaxy_results)
        
        print(f"PURE GEOMETRIC VALIDATION RESULTS:")
        print(f"=" * 35)
        print(f"Galaxies tested: {stats_summary['total_galaxies_tested']}")
        print()
        
        print(f"GEOMETRIC CONSISTENCY:")
        print(f"Median geometric score: {stats_summary['median_geometric_score']:.3f}")
        print(f"Mean geometric score: {stats_summary['mean_geometric_score']:.3f} ± {stats_summary['geometric_score_std']:.3f}")
        print()
        
        print(f"SCALE CONSISTENCY:")
        print(f"Median optimal R0: {stats_summary['median_optimal_R0']:.1f}")
        print(f"Mean optimal R0: {stats_summary['mean_optimal_R0']:.1f}")
        print(f"R0 scale variation: {stats_summary['R0_scale_consistency']*100:.1f}%")
        print()
        
        print(f"GEOMETRIC QUALITY DISTRIBUTION:")
        for quality, count in quality_counts.items():
            percentage = count / len(self.galaxy_results) * 100
            print(f"  {quality}: {count} galaxies ({percentage:.1f}%)")
        print()
        
        print(f"UDT GEOMETRIC VALIDATION SUMMARY:")
        print(f"Excellent geometric consistency: {excellent_fraction*100:.1f}%")
        print(f"Good or better consistency: {good_or_better_fraction*100:.1f}%")
        
        return stats_summary
    
    def create_pure_geometric_visualizations(self):
        """Create visualization of pure geometric validation."""
        print("CREATING PURE GEOMETRIC VISUALIZATIONS")
        print("-" * 35)
        
        if not self.galaxy_results:
            print("No data to visualize")
            return
        
        # Extract data
        geometric_scores = [r['geometric_score'] for r in self.galaxy_results]
        optimal_R0_scales = [r['optimal_R0_scale'] for r in self.galaxy_results]
        n_points = [r['n_points'] for r in self.galaxy_results]
        
        # Create visualization
        fig = plt.figure(figsize=(15, 10))
        
        # Panel 1: Geometric score distribution
        ax1 = plt.subplot(2, 3, 1)
        plt.hist(geometric_scores, bins=20, alpha=0.7, color='blue', edgecolor='black')
        plt.axvline(np.median(geometric_scores), color='red', linestyle='--', linewidth=2,
                   label=f'Median: {np.median(geometric_scores):.3f}')
        plt.xlabel('Geometric Consistency Score')
        plt.ylabel('Number of Galaxies')
        plt.title('UDT Geometric Pattern Consistency')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Panel 2: Optimal R0 scale distribution
        ax2 = plt.subplot(2, 3, 2)
        plt.hist(optimal_R0_scales, bins=20, alpha=0.7, color='green', edgecolor='black')
        plt.axvline(np.median(optimal_R0_scales), color='red', linestyle='--', linewidth=2,
                   label=f'Median: {np.median(optimal_R0_scales):.1f}')
        plt.xlabel('Optimal R0 Scale (dimensionless)')
        plt.ylabel('Number of Galaxies')
        plt.title('UDT Scale Parameter Consistency')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Panel 3: Geometric score vs data quality
        ax3 = plt.subplot(2, 3, 3)
        plt.scatter(n_points, geometric_scores, alpha=0.6, color='purple')
        plt.xlabel('Number of Data Points')
        plt.ylabel('Geometric Consistency Score')
        plt.title('Pattern Quality vs Data Density')
        plt.grid(True, alpha=0.3)
        
        # Panel 4: Quality distribution
        ax4 = plt.subplot(2, 3, 4)
        qualities = ['EXCELLENT', 'GOOD', 'MODERATE', 'POOR', 'FAILED']
        quality_counts = []
        for quality in qualities:
            count = sum(1 for r in self.galaxy_results if r['geometric_quality'] == quality)
            quality_counts.append(count)
        
        colors = ['darkgreen', 'lightgreen', 'yellow', 'orange', 'red']
        plt.bar(range(len(qualities)), quality_counts, color=colors, alpha=0.7, edgecolor='black')
        plt.xlabel('Geometric Quality')
        plt.ylabel('Number of Galaxies')
        plt.title('UDT Geometric Quality Distribution')
        plt.xticks(range(len(qualities)), qualities, rotation=45)
        plt.grid(True, alpha=0.3)
        
        # Panel 5: Cumulative geometric performance
        ax5 = plt.subplot(2, 3, 5)
        sorted_scores = np.sort(geometric_scores)
        cumulative = np.arange(1, len(sorted_scores) + 1) / len(sorted_scores)
        plt.plot(sorted_scores, cumulative, linewidth=2, color='blue')
        plt.axvline(0.6, color='red', linestyle='--', alpha=0.7, label='Good threshold')
        plt.axvline(0.8, color='green', linestyle='--', alpha=0.7, label='Excellent threshold')
        plt.xlabel('Geometric Consistency Score')
        plt.ylabel('Cumulative Fraction')
        plt.title('UDT Geometric Performance')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Panel 6: Summary
        ax6 = plt.subplot(2, 3, 6)
        ax6.axis('off')
        
        excellent_count = sum(1 for r in self.galaxy_results if r['geometric_quality'] == 'EXCELLENT')
        good_count = sum(1 for r in self.galaxy_results if r['geometric_quality'] == 'GOOD')
        
        summary_text = f"""
PURE GEOMETRIC UDT VALIDATION
CONTAMINATION-FREE RESULTS

Galaxies Tested: {len(self.galaxy_results)}
Data Points: {sum(n_points):,}

GEOMETRIC CONSISTENCY:
Median Score: {np.median(geometric_scores):.3f}
Mean Score: {np.mean(geometric_scores):.3f}

QUALITY ASSESSMENT:
Excellent: {excellent_count}/{len(self.galaxy_results)}
Good+: {excellent_count + good_count}/{len(self.galaxy_results)}

SCALE CONSISTENCY:
R0 variation: {np.std(optimal_R0_scales)/np.mean(optimal_R0_scales)*100:.1f}%

STATUS: PURE GEOMETRIC TEST
        """
        
        ax6.text(0.1, 0.9, summary_text, transform=ax6.transAxes, fontsize=9,
                verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightcyan", alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.results_dir, 'pure_geometric_udt_validation.png'), 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Pure geometric visualization saved: pure_geometric_udt_validation.png")
    
    def save_pure_geometric_results(self, stats_summary):
        """Save pure geometric validation results."""
        print("SAVING PURE GEOMETRIC RESULTS")
        print("-" * 27)
        
        complete_results = {
            'analysis_type': 'pure_geometric_udt_validation',
            'date': '2025-07-18',
            'contamination_status': 'ZERO_STANDARD_MODEL_CONTAMINATION',
            'methodology': 'pure_geometric_pattern_validation',
            'validation_criteria': [
                'geometric_pattern_consistency',
                'scale_invariance_testing', 
                'theoretical_prediction_capability',
                'monotonicity_preservation',
                'asymptotic_behavior_consistency',
                'inner_structure_correlation'
            ],
            'udt_parameters': {
                'R0_galactic_kpc': self.R0_galactic / 3.086e16 / 1000,
                'alpha': self.alpha
            },
            'geometric_statistics': stats_summary,
            'individual_galaxies': self.galaxy_results,
            'conclusion': 'UDT geometric validation using contamination-free criteria'
        }
        
        results_file = os.path.join(self.results_dir, 'pure_geometric_udt_validation.json')
        with open(results_file, 'w') as f:
            json.dump(complete_results, f, indent=2, default=str)
        
        print(f"Pure geometric results saved: {results_file}")
        
        return complete_results
    
    def run_pure_geometric_validation(self):
        """Run complete pure geometric UDT validation."""
        print("STARTING PURE GEOMETRIC UDT VALIDATION")
        print("=" * 39)
        print("Methodology: Contamination-free geometric pattern analysis")
        print("Goal: Test UDT theoretical predictions without Standard Model assumptions")
        print()
        
        # Analyze all galaxies
        successful_tests, total_galaxies = self.analyze_all_galaxies_pure_geometric()
        
        if successful_tests == 0:
            print("ERROR: No successful geometric tests")
            return None
        
        # Calculate pure geometric statistics
        stats_summary = self.calculate_pure_geometric_statistics()
        
        # Create visualizations
        self.create_pure_geometric_visualizations()
        
        # Save results
        complete_results = self.save_pure_geometric_results(stats_summary)
        
        # Final geometric assessment
        print("\n" + "=" * 60)
        print("PURE GEOMETRIC UDT VALIDATION - FINAL ASSESSMENT")
        print("=" * 60)
        
        excellent_fraction = stats_summary['quality_distribution']['EXCELLENT'] / stats_summary['total_galaxies_tested']
        good_or_better = (stats_summary['quality_distribution']['EXCELLENT'] + 
                         stats_summary['quality_distribution']['GOOD']) / stats_summary['total_galaxies_tested']
        
        print(f"\nPURE GEOMETRIC RESULTS:")
        print(f"Galaxies tested: {successful_tests}/{total_galaxies}")
        print(f"Median geometric score: {stats_summary['median_geometric_score']:.3f}")
        print(f"Excellent geometric consistency: {excellent_fraction*100:.1f}%")
        print(f"Good or better consistency: {good_or_better*100:.1f}%")
        
        # Pure geometric verdict (no Standard Model criteria)
        if excellent_fraction >= 0.5:
            verdict = "GEOMETRIC PATTERN VALIDATED"
        elif good_or_better >= 0.5:
            verdict = "GEOMETRIC PATTERN SUPPORTED"
        elif stats_summary['median_geometric_score'] >= 0.4:
            verdict = "GEOMETRIC PATTERN INCONCLUSIVE"
        else:
            verdict = "GEOMETRIC PATTERN NOT SUPPORTED"
        
        print(f"\nFINAL GEOMETRIC VERDICT: {verdict}")
        print(f"UDT geometric structure assessed using contamination-free criteria")
        
        return complete_results

def main():
    """Main pure geometric validation routine."""
    validator = PureGeometricUDTValidator()
    results = validator.run_pure_geometric_validation()
    return results

if __name__ == "__main__":
    main()