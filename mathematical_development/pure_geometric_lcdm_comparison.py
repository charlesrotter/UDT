#!/usr/bin/env python3
"""
Pure Geometric ΛCDM vs UDT Comparison
====================================

CRITICAL: Compare ΛCDM and UDT using identical contamination-free geometric criteria.

METHODOLOGY:
- Same geometric pattern tests for both theories
- NO statistical fitting contamination
- Direct comparison of theoretical predictions
- Geometric consistency scoring only

ΛCDM ROTATION CURVE PREDICTION:
- Stellar disk component: v_disk(r) ∝ sqrt(M_disk(r)/r)
- Dark matter halo: v_dm(r) = v_200 * sqrt(r/r_s / (1 + r/r_s)²) (NFW profile)
- Total: v_total² = v_disk² + v_dm²

UDT ROTATION CURVE PREDICTION:
- Pure geometric: v(r) ∝ sqrt(F(τ) * M(r)/r) where τ(r) = R₀/(R₀ + r)
- Enhancement: F(τ) = 1 + α × 3(1-τ)/(τ²(3-2τ))

Author: Charles Rotter
Date: 2025-07-18
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import glob
import json

class PureGeometricLCDMUDTComparison:
    def __init__(self):
        self.data_dir = "C:/UDT/data/sparc_database/"
        self.results_dir = "C:/UDT/results/"
        
        os.makedirs(self.results_dir, exist_ok=True)
        
        print("PURE GEOMETRIC ΛCDM vs UDT COMPARISON")
        print("=" * 37)
        print("GOAL: Compare theories using identical contamination-free criteria")
        print("METHOD: Geometric pattern consistency testing")
        print("CRITERIA: Theory-agnostic geometric validation")
        print()
        
        # UDT parameters
        self.R0_galactic = 57.5  # kpc (dimensionless scale)
        self.alpha = 0.00206
        
        # ΛCDM parameters (typical values)
        self.v_200_typical = 200.0  # km/s (halo circular velocity)
        self.r_s_typical = 20.0     # kpc (halo scale radius)
        
        print(f"Theory Parameters:")
        print(f"UDT: R₀ = {self.R0_galactic:.1f} kpc, α = {self.alpha:.5f}")
        print(f"ΛCDM: v₂₀₀ = {self.v_200_typical:.1f} km/s, r_s = {self.r_s_typical:.1f} kpc")
        print()
        
        self.comparison_results = []
        
    def load_galaxy_observational_pattern(self, filename):
        """Load raw observational pattern (same as UDT analysis)."""
        filepath = os.path.join(self.data_dir, filename)
        
        try:
            data = np.loadtxt(filepath)
            
            if data.size == 0:
                return None, None
            
            if data.ndim == 1:
                data = data.reshape(1, -1)
            
            radius_kpc = data[:, 0]
            velocity_obs = data[:, 1]
            
            # Filter valid data
            valid_mask = (radius_kpc > 0) & (velocity_obs > 0)
            
            if not np.any(valid_mask):
                return None, None
            
            radius_kpc = radius_kpc[valid_mask]
            velocity_obs = velocity_obs[valid_mask]
            
            # Create dimensionless pattern
            r_max = np.max(radius_kpc)
            r_normalized = radius_kpc / r_max
            
            v_flat = np.median(velocity_obs[radius_kpc > 0.7 * r_max])
            v_normalized = velocity_obs / v_flat
            
            return r_normalized, v_normalized
            
        except Exception as e:
            print(f"Error loading {filename}: {e}")
            return None, None
    
    def calculate_lcdm_geometric_prediction(self, r_normalized, v_200_scale, r_s_scale):
        """Calculate ΛCDM rotation curve prediction."""
        # Convert normalized radius to physical scale
        r_max_guess = 30.0  # kpc (typical galaxy size)
        r_physical = r_normalized * r_max_guess
        
        # ΛCDM NFW dark matter profile
        x = r_physical / r_s_scale
        v_dm_squared = v_200_scale**2 * (np.log(1 + x) - x/(1 + x)) / x
        
        # Stellar disk component (exponential disk)
        # Simplified: v_disk ∝ sqrt(r) for inner region, flat for outer
        v_disk_inner = v_200_scale * 0.6 * np.sqrt(r_normalized)  # Inner rising
        v_disk_outer = v_200_scale * 0.6 * np.ones_like(r_normalized)  # Outer flat
        
        # Smooth transition
        transition = 1 / (1 + np.exp(-10 * (r_normalized - 0.3)))
        v_disk_squared = (1 - transition) * v_disk_inner**2 + transition * v_disk_outer**2
        
        # Total ΛCDM velocity
        v_total_squared = v_disk_squared + np.maximum(0, v_dm_squared)
        v_lcdm = np.sqrt(v_total_squared)
        
        # Normalize to asymptotic value
        v_asymptotic = v_lcdm[-1] if len(v_lcdm) > 0 else 1.0
        v_lcdm_normalized = v_lcdm / v_asymptotic
        
        return v_lcdm_normalized
    
    def calculate_udt_geometric_prediction(self, r_normalized, R0_scale):
        """Calculate UDT rotation curve prediction (same as pure validation)."""
        # Convert normalized radius to physical scale
        r_physical = r_normalized * R0_scale
        
        # UDT geometry
        tau_r = R0_scale / (R0_scale + r_physical)
        F_tau = 1 + self.alpha * 3 * (1 - tau_r) / (tau_r**2 * (3 - 2*tau_r))
        
        # UDT velocity structure
        v_udt_structure = np.sqrt(F_tau * r_normalized / (1 + r_normalized))
        
        # Normalize
        v_asymptotic = v_udt_structure[-1] if len(v_udt_structure) > 0 else 1.0
        v_udt_normalized = v_udt_structure / v_asymptotic
        
        return v_udt_normalized
    
    def geometric_pattern_consistency_test(self, r_obs, v_obs, v_theory, theory_name):
        """Apply identical geometric tests to both theories."""
        
        # Test 1: Monotonicity preservation
        obs_monotonic = []
        theory_monotonic = []
        
        for i in range(1, len(v_obs)):
            if v_obs[i] > v_obs[i-1]:
                obs_monotonic.append('inc')
            elif v_obs[i] < v_obs[i-1]:
                obs_monotonic.append('dec')
            else:
                obs_monotonic.append('const')
                
            if v_theory[i] > v_theory[i-1]:
                theory_monotonic.append('inc')
            elif v_theory[i] < v_theory[i-1]:
                theory_monotonic.append('dec')
            else:
                theory_monotonic.append('const')
        
        monotonic_agreement = sum(1 for obs, th in zip(obs_monotonic, theory_monotonic) 
                                if obs == th) / len(obs_monotonic)
        
        # Test 2: Asymptotic behavior
        outer_region = r_obs > 0.8
        if np.sum(outer_region) > 2:
            obs_asymptotic = np.mean(v_obs[outer_region])
            theory_asymptotic = np.mean(v_theory[outer_region])
            asymptotic_consistency = 1.0 - abs(obs_asymptotic - theory_asymptotic) / max(obs_asymptotic, theory_asymptotic)
        else:
            asymptotic_consistency = 1.0
        
        # Test 3: Inner region structure correlation
        inner_region = r_obs < 0.3
        if np.sum(inner_region) > 2:
            obs_inner_grad = np.gradient(v_obs[inner_region])
            theory_inner_grad = np.gradient(v_theory[inner_region])
            
            if len(obs_inner_grad) > 1 and np.std(obs_inner_grad) > 0 and np.std(theory_inner_grad) > 0:
                inner_correlation = np.corrcoef(obs_inner_grad, theory_inner_grad)[0, 1]
                if np.isnan(inner_correlation):
                    inner_correlation = 0.0
            else:
                inner_correlation = 0.0
        else:
            inner_correlation = 0.0
        
        # Overall geometric consistency score
        geometric_score = (
            monotonic_agreement * 0.4 +
            asymptotic_consistency * 0.3 +
            max(0, inner_correlation) * 0.3
        )
        
        return {
            'theory': theory_name,
            'geometric_score': geometric_score,
            'monotonic_agreement': monotonic_agreement,
            'asymptotic_consistency': asymptotic_consistency,
            'inner_correlation': inner_correlation
        }
    
    def optimize_theory_parameters(self, r_obs, v_obs, theory_type):
        """Find optimal parameters for each theory using geometric criteria."""
        
        best_score = -1
        best_params = None
        best_prediction = None
        
        if theory_type == 'LCDM':
            # Test ΛCDM parameter space
            v_200_range = np.linspace(100, 300, 20)  # km/s
            r_s_range = np.linspace(10, 40, 20)      # kpc
            
            for v_200 in v_200_range:
                for r_s in r_s_range:
                    try:
                        v_theory = self.calculate_lcdm_geometric_prediction(r_obs, v_200, r_s)
                        
                        if len(v_theory) == len(v_obs):
                            test_result = self.geometric_pattern_consistency_test(r_obs, v_obs, v_theory, 'LCDM')
                            
                            if test_result['geometric_score'] > best_score:
                                best_score = test_result['geometric_score']
                                best_params = {'v_200': v_200, 'r_s': r_s}
                                best_prediction = v_theory.copy()
                    except:
                        continue
        
        elif theory_type == 'UDT':
            # Test UDT parameter space
            R0_range = np.logspace(0, 2, 20)  # 1 to 100 (dimensionless)
            
            for R0 in R0_range:
                try:
                    v_theory = self.calculate_udt_geometric_prediction(r_obs, R0)
                    
                    if len(v_theory) == len(v_obs):
                        test_result = self.geometric_pattern_consistency_test(r_obs, v_obs, v_theory, 'UDT')
                        
                        if test_result['geometric_score'] > best_score:
                            best_score = test_result['geometric_score']
                            best_params = {'R0': R0}
                            best_prediction = v_theory.copy()
                except:
                    continue
        
        return best_score, best_params, best_prediction
    
    def compare_galaxy_theories(self, galaxy_name, r_obs, v_obs):
        """Compare ΛCDM and UDT for single galaxy using geometric criteria."""
        
        print(f"  Comparing theories for {galaxy_name}")
        
        # Optimize both theories
        lcdm_score, lcdm_params, lcdm_prediction = self.optimize_theory_parameters(r_obs, v_obs, 'LCDM')
        udt_score, udt_params, udt_prediction = self.optimize_theory_parameters(r_obs, v_obs, 'UDT')
        
        # Determine geometric winner
        if lcdm_score > udt_score:
            winner = 'LCDM'
            advantage = lcdm_score - udt_score
        elif udt_score > lcdm_score:
            winner = 'UDT'
            advantage = udt_score - lcdm_score
        else:
            winner = 'TIE'
            advantage = 0.0
        
        print(f"    ΛCDM geometric score: {lcdm_score:.3f}")
        print(f"    UDT geometric score: {udt_score:.3f}")
        print(f"    Geometric winner: {winner} (advantage: {advantage:.3f})")
        
        return {
            'galaxy_name': galaxy_name,
            'lcdm_score': lcdm_score,
            'udt_score': udt_score,
            'winner': winner,
            'advantage': advantage,
            'lcdm_params': lcdm_params,
            'udt_params': udt_params,
            'n_points': len(r_obs)
        }
    
    def analyze_all_galaxies_comparison(self):
        """Compare ΛCDM vs UDT on all SPARC galaxies."""
        print("PURE GEOMETRIC ΛCDM vs UDT COMPARISON")
        print("-" * 39)
        
        galaxy_files = glob.glob(os.path.join(self.data_dir, "*_rotmod.dat"))
        
        print(f"Found {len(galaxy_files)} galaxy files")
        print("Comparing ΛCDM and UDT geometric predictions...")
        print()
        
        successful_comparisons = 0
        total_galaxies = 0
        
        lcdm_wins = 0
        udt_wins = 0
        ties = 0
        
        for i, filepath in enumerate(galaxy_files):
            filename = os.path.basename(filepath)
            galaxy_name = filename.replace('_rotmod.dat', '')
            
            total_galaxies += 1
            
            print(f"Processing {i+1}/{len(galaxy_files)}: {galaxy_name}")
            
            # Load observational pattern
            r_obs, v_obs = self.load_galaxy_observational_pattern(filename)
            
            if r_obs is None:
                print(f"  Failed to load observational pattern")
                continue
            
            print(f"  Data points: {len(r_obs)}")
            
            # Compare theories
            try:
                result = self.compare_galaxy_theories(galaxy_name, r_obs, v_obs)
                
                self.comparison_results.append(result)
                successful_comparisons += 1
                
                if result['winner'] == 'LCDM':
                    lcdm_wins += 1
                elif result['winner'] == 'UDT':
                    udt_wins += 1
                else:
                    ties += 1
                
            except Exception as e:
                print(f"  Comparison failed: {e}")
        
        print()
        print(f"GEOMETRIC COMPARISON COMPLETE:")
        print(f"Total galaxies processed: {total_galaxies}")
        print(f"Successful comparisons: {successful_comparisons}")
        print(f"ΛCDM geometric wins: {lcdm_wins}")
        print(f"UDT geometric wins: {udt_wins}")
        print(f"Ties: {ties}")
        print()
        
        return successful_comparisons, total_galaxies, lcdm_wins, udt_wins, ties
    
    def calculate_comparison_statistics(self):
        """Calculate comprehensive comparison statistics."""
        print("PURE GEOMETRIC COMPARISON STATISTICS")
        print("-" * 35)
        
        if not self.comparison_results:
            print("No comparison results to analyze")
            return None
        
        # Extract scores
        lcdm_scores = [r['lcdm_score'] for r in self.comparison_results]
        udt_scores = [r['udt_score'] for r in self.comparison_results]
        advantages = [r['advantage'] for r in self.comparison_results]
        
        # Count winners
        lcdm_wins = sum(1 for r in self.comparison_results if r['winner'] == 'LCDM')
        udt_wins = sum(1 for r in self.comparison_results if r['winner'] == 'UDT')
        ties = sum(1 for r in self.comparison_results if r['winner'] == 'TIE')
        
        total = len(self.comparison_results)
        
        stats_summary = {
            'total_comparisons': total,
            'lcdm_wins': lcdm_wins,
            'udt_wins': udt_wins,
            'ties': ties,
            'lcdm_win_rate': lcdm_wins / total,
            'udt_win_rate': udt_wins / total,
            'tie_rate': ties / total,
            'median_lcdm_score': np.median(lcdm_scores),
            'median_udt_score': np.median(udt_scores),
            'mean_lcdm_score': np.mean(lcdm_scores),
            'mean_udt_score': np.mean(udt_scores),
            'lcdm_score_std': np.std(lcdm_scores),
            'udt_score_std': np.std(udt_scores),
            'median_advantage': np.median(advantages),
            'mean_advantage': np.mean(advantages)
        }
        
        print(f"CONTAMINATION-FREE ΛCDM vs UDT COMPARISON:")
        print(f"=" * 43)
        print(f"Total galaxy comparisons: {total}")
        print()
        
        print(f"GEOMETRIC PERFORMANCE SCORES:")
        print(f"ΛCDM median: {stats_summary['median_lcdm_score']:.3f}")
        print(f"UDT median: {stats_summary['median_udt_score']:.3f}")
        print(f"ΛCDM mean: {stats_summary['mean_lcdm_score']:.3f} ± {stats_summary['lcdm_score_std']:.3f}")
        print(f"UDT mean: {stats_summary['mean_udt_score']:.3f} ± {stats_summary['udt_score_std']:.3f}")
        print()
        
        print(f"GEOMETRIC COMPETITION RESULTS:")
        print(f"ΛCDM wins: {lcdm_wins}/{total} ({stats_summary['lcdm_win_rate']*100:.1f}%)")
        print(f"UDT wins: {udt_wins}/{total} ({stats_summary['udt_win_rate']*100:.1f}%)")
        print(f"Ties: {ties}/{total} ({stats_summary['tie_rate']*100:.1f}%)")
        print()
        
        # Determine overall geometric winner
        if stats_summary['udt_win_rate'] > stats_summary['lcdm_win_rate']:
            overall_winner = "UDT"
            win_margin = stats_summary['udt_win_rate'] - stats_summary['lcdm_win_rate']
        elif stats_summary['lcdm_win_rate'] > stats_summary['udt_win_rate']:
            overall_winner = "ΛCDM"
            win_margin = stats_summary['lcdm_win_rate'] - stats_summary['udt_win_rate']
        else:
            overall_winner = "TIE"
            win_margin = 0.0
        
        print(f"OVERALL GEOMETRIC WINNER: {overall_winner}")
        if win_margin > 0:
            print(f"Win margin: {win_margin*100:.1f} percentage points")
        
        return stats_summary
    
    def create_comparison_visualizations(self):
        """Create visualization comparing ΛCDM and UDT geometric performance."""
        print("CREATING COMPARISON VISUALIZATIONS")
        print("-" * 31)
        
        if not self.comparison_results:
            print("No data to visualize")
            return
        
        # Extract data
        lcdm_scores = [r['lcdm_score'] for r in self.comparison_results]
        udt_scores = [r['udt_score'] for r in self.comparison_results]
        advantages = [r['advantage'] for r in self.comparison_results]
        
        fig = plt.figure(figsize=(15, 12))
        
        # Panel 1: Score comparison scatter
        ax1 = plt.subplot(2, 3, 1)
        plt.scatter(lcdm_scores, udt_scores, alpha=0.6, color='purple')
        plt.plot([0, 1], [0, 1], 'k--', alpha=0.5, label='Equal performance')
        plt.xlabel('ΛCDM Geometric Score')
        plt.ylabel('UDT Geometric Score')
        plt.title('Direct Geometric Performance Comparison')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Panel 2: Score distributions
        ax2 = plt.subplot(2, 3, 2)
        plt.hist(lcdm_scores, bins=20, alpha=0.6, color='red', label='ΛCDM', density=True)
        plt.hist(udt_scores, bins=20, alpha=0.6, color='blue', label='UDT', density=True)
        plt.axvline(np.median(lcdm_scores), color='red', linestyle='--', 
                   label=f'ΛCDM median: {np.median(lcdm_scores):.3f}')
        plt.axvline(np.median(udt_scores), color='blue', linestyle='--',
                   label=f'UDT median: {np.median(udt_scores):.3f}')
        plt.xlabel('Geometric Score')
        plt.ylabel('Density')
        plt.title('Geometric Score Distributions')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Panel 3: Win distribution
        ax3 = plt.subplot(2, 3, 3)
        winners = [r['winner'] for r in self.comparison_results]
        winner_counts = {'ΛCDM': winners.count('ΛCDM'), 
                        'UDT': winners.count('UDT'), 
                        'TIE': winners.count('TIE')}
        
        colors = ['red', 'blue', 'gray']
        plt.bar(winner_counts.keys(), winner_counts.values(), color=colors, alpha=0.7)
        plt.ylabel('Number of Galaxies')
        plt.title('Geometric Competition Winners')
        plt.grid(True, alpha=0.3)
        
        # Panel 4: Advantage distribution
        ax4 = plt.subplot(2, 3, 4)
        plt.hist(advantages, bins=30, alpha=0.7, color='green', edgecolor='black')
        plt.axvline(0, color='black', linestyle='-', linewidth=2, label='Equal performance')
        plt.axvline(np.median(advantages), color='red', linestyle='--', linewidth=2,
                   label=f'Median advantage: {np.median(advantages):.3f}')
        plt.xlabel('Geometric Advantage (UDT - ΛCDM)')
        plt.ylabel('Number of Galaxies')
        plt.title('Geometric Performance Advantage')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Panel 5: Cumulative comparison
        ax5 = plt.subplot(2, 3, 5)
        sorted_lcdm = np.sort(lcdm_scores)
        sorted_udt = np.sort(udt_scores)
        cum_lcdm = np.arange(1, len(sorted_lcdm) + 1) / len(sorted_lcdm)
        cum_udt = np.arange(1, len(sorted_udt) + 1) / len(sorted_udt)
        
        plt.plot(sorted_lcdm, cum_lcdm, color='red', linewidth=2, label='ΛCDM')
        plt.plot(sorted_udt, cum_udt, color='blue', linewidth=2, label='UDT')
        plt.xlabel('Geometric Score')
        plt.ylabel('Cumulative Fraction')
        plt.title('Cumulative Geometric Performance')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Panel 6: Summary statistics
        ax6 = plt.subplot(2, 3, 6)
        ax6.axis('off')
        
        lcdm_wins = sum(1 for r in self.comparison_results if r['winner'] == 'ΛCDM')
        udt_wins = sum(1 for r in self.comparison_results if r['winner'] == 'UDT')
        ties = sum(1 for r in self.comparison_results if r['winner'] == 'TIE')
        
        if udt_wins > lcdm_wins:
            overall_winner = "UDT"
        elif lcdm_wins > udt_wins:
            overall_winner = "ΛCDM"
        else:
            overall_winner = "TIE"
        
        summary_text = f"""
CONTAMINATION-FREE COMPARISON
ΛCDM vs UDT

Galaxies Compared: {len(self.comparison_results)}

GEOMETRIC SCORES:
ΛCDM: {np.median(lcdm_scores):.3f} median
UDT: {np.median(udt_scores):.3f} median

COMPETITION RESULTS:
ΛCDM wins: {lcdm_wins}
UDT wins: {udt_wins}
Ties: {ties}

OVERALL WINNER: {overall_winner}

STATUS: PURE GEOMETRIC TEST
        """
        
        ax6.text(0.1, 0.9, summary_text, transform=ax6.transAxes, fontsize=10,
                verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow", alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.results_dir, 'pure_geometric_lcdm_udt_comparison.png'), 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Comparison visualization saved: pure_geometric_lcdm_udt_comparison.png")
    
    def save_comparison_results(self, stats_summary):
        """Save comparison results."""
        print("SAVING COMPARISON RESULTS")
        print("-" * 23)
        
        complete_results = {
            'analysis_type': 'pure_geometric_lcdm_udt_comparison',
            'date': '2025-07-18',
            'contamination_status': 'ZERO_STANDARD_MODEL_STATISTICAL_CONTAMINATION',
            'methodology': 'identical_geometric_criteria_for_both_theories',
            'comparison_summary': stats_summary,
            'individual_comparisons': self.comparison_results,
            'conclusion': 'Head-to-head geometric comparison of ΛCDM and UDT'
        }
        
        results_file = os.path.join(self.results_dir, 'pure_geometric_lcdm_udt_comparison.json')
        with open(results_file, 'w') as f:
            json.dump(complete_results, f, indent=2, default=str)
        
        print(f"Comparison results saved: {results_file}")
        
        return complete_results
    
    def run_pure_geometric_comparison(self):
        """Run complete ΛCDM vs UDT geometric comparison."""
        print("STARTING PURE GEOMETRIC ΛCDM vs UDT COMPARISON")
        print("=" * 48)
        print("Methodology: Identical contamination-free geometric criteria")
        print("Goal: Determine which theory better describes galactic geometry")
        print()
        
        # Compare all galaxies
        successful, total, lcdm_wins, udt_wins, ties = self.analyze_all_galaxies_comparison()
        
        if successful == 0:
            print("ERROR: No successful comparisons")
            return None
        
        # Calculate statistics
        stats_summary = self.calculate_comparison_statistics()
        
        # Create visualizations
        self.create_comparison_visualizations()
        
        # Save results
        complete_results = self.save_comparison_results(stats_summary)
        
        # Final verdict
        print("\n" + "=" * 60)
        print("PURE GEOMETRIC ΛCDM vs UDT COMPARISON - FINAL VERDICT")
        print("=" * 60)
        
        print(f"\nCOMPARISON RESULTS:")
        print(f"Galaxies compared: {successful}/{total}")
        print(f"ΛCDM geometric wins: {lcdm_wins} ({lcdm_wins/successful*100:.1f}%)")
        print(f"UDT geometric wins: {udt_wins} ({udt_wins/successful*100:.1f}%)")
        print(f"Ties: {ties} ({ties/successful*100:.1f}%)")
        
        # Determine winner
        if udt_wins > lcdm_wins:
            winner = "UDT"
            margin = udt_wins - lcdm_wins
        elif lcdm_wins > udt_wins:
            winner = "ΛCDM"
            margin = lcdm_wins - udt_wins
        else:
            winner = "TIE"
            margin = 0
        
        print(f"\nGEOMETRIC WINNER: {winner}")
        if margin > 0:
            print(f"Victory margin: {margin} galaxies ({margin/successful*100:.1f} percentage points)")
        
        print(f"Both theories evaluated using identical contamination-free criteria")
        
        return complete_results

def main():
    """Main comparison routine."""
    comparator = PureGeometricLCDMUDTComparison()
    results = comparator.run_pure_geometric_comparison()
    return results

if __name__ == "__main__":
    main()