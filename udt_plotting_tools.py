"""
UDT Plotting Tools - Publication Quality Figures
Generate professional plots for UDT theory validation and results

Usage:
    from udt_plotting_tools import UDTPlotter
    plotter = UDTPlotter()
    plotter.plot_rotation_curve_comparison(galaxy_data, udt_prediction)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from udt_core_validated import UDTFramework, SPARCAnalysis

class UDTPlotter:
    """Publication-quality plotting tools for UDT theory"""
    
    def __init__(self, style: str = 'publication'):
        """
        Initialize plotter with publication settings
        
        Args:
            style: 'publication', 'presentation', or 'default'
        """
        self.setup_style(style)
        
    def setup_style(self, style: str):
        """Configure matplotlib for publication quality"""
        
        if style == 'publication':
            plt.rcParams.update({
                'font.size': 12,
                'font.family': 'serif',
                'font.serif': ['Computer Modern Roman'],
                'text.usetex': False,  # Set to True if LaTeX available
                'figure.figsize': (8, 6),
                'figure.dpi': 300,
                'axes.linewidth': 1.2,
                'axes.labelsize': 14,
                'axes.titlesize': 16,
                'xtick.labelsize': 12,
                'ytick.labelsize': 12,
                'legend.fontsize': 12,
                'lines.linewidth': 2,
                'grid.alpha': 0.3
            })
        elif style == 'presentation':
            plt.rcParams.update({
                'font.size': 14,
                'axes.labelsize': 16,
                'axes.titlesize': 18,
                'xtick.labelsize': 14,
                'ytick.labelsize': 14,
                'legend.fontsize': 14,
                'lines.linewidth': 3,
                'figure.figsize': (10, 7)
            })
    
    def plot_rotation_curve_comparison(self, galaxy: Dict, 
                                     sparc_analysis: SPARCAnalysis,
                                     save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot rotation curve: UDT prediction vs observations
        
        Args:
            galaxy: Galaxy data dictionary
            sparc_analysis: SPARC analysis object
            save_path: Optional path to save figure
            
        Returns:
            matplotlib Figure object
        """
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), 
                                      gridspec_kw={'height_ratios': [3, 1]})
        
        # Generate UDT prediction
        radii = np.linspace(0.5, 25, 100)
        _, velocities = sparc_analysis.predict_rotation_curve(galaxy, radii)
        
        # Main plot
        ax1.plot(radii, velocities, 'b-', linewidth=3, label='UDT Prediction', zorder=3)
        
        # Add observed point if available
        if 'vflat' in galaxy and not np.isnan(galaxy['vflat']):
            # Plot at effective radius as approximation
            r_eff = galaxy.get('reff', 8.0)  # Default to 8 kpc if not available
            v_err = galaxy.get('vflat_err', galaxy['vflat'] * 0.1)
            
            ax1.errorbar(r_eff, galaxy['vflat'], yerr=v_err, 
                        fmt='ro', markersize=8, capsize=5, capthick=2,
                        label='SPARC Observation', zorder=4)
        
        # Formatting
        ax1.set_ylabel('Rotation Velocity (km/s)', fontsize=14)
        ax1.set_title(f'{galaxy["name"]} - UDT Rotation Curve Prediction', fontsize=16)
        ax1.legend(loc='upper right')
        ax1.grid(True, alpha=0.3)
        ax1.set_xlim(0, 25)
        ax1.set_ylim(0, max(velocities) * 1.1)
        
        # Residual plot
        if 'vflat' in galaxy and not np.isnan(galaxy['vflat']):
            udt_framework = sparc_analysis.udt
            M_star = udt_framework.luminosity_to_mass(galaxy['luminosity'])
            v_predicted = udt_framework.rotation_velocity_udt(r_eff, M_star)
            residual = v_predicted - galaxy['vflat']
            
            ax2.errorbar(r_eff, residual, yerr=v_err, 
                        fmt='ro', markersize=6, capsize=3)
            ax2.axhline(0, color='black', linestyle='--', alpha=0.7)
            ax2.set_ylabel('Residual\n(km/s)', fontsize=12)
            ax2.set_xlabel('Radius (kpc)', fontsize=14)
            ax2.grid(True, alpha=0.3)
            ax2.set_xlim(0, 25)
        else:
            ax2.set_xlabel('Radius (kpc)', fontsize=14)
            ax2.text(0.5, 0.5, 'No observational data for residuals', 
                    transform=ax2.transAxes, ha='center', va='center',
                    fontsize=12, style='italic')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        return fig
    
    def plot_velocity_scaling_analysis(self, results_df: pd.DataFrame,
                                     save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot velocity ratio vs galaxy properties for scaling analysis
        """
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
        
        # 1. Velocity ratio vs stellar mass
        ax1.scatter(results_df['stellar_mass_1e9'], results_df['ratio_pred_obs'], 
                   alpha=0.6, s=40, c=results_df['chi_squared'], cmap='viridis_r')
        ax1.axhline(1.0, color='red', linestyle='--', alpha=0.7)
        ax1.set_xlabel('Stellar Mass (10⁹ M☉)')
        ax1.set_ylabel('Velocity Ratio (Predicted/Observed)')
        ax1.set_title('UDT Performance vs Stellar Mass')
        ax1.set_xscale('log')
        ax1.grid(True, alpha=0.3)
        
        # 2. Velocity ratio vs luminosity
        ax2.scatter(results_df['luminosity'], results_df['ratio_pred_obs'], 
                   alpha=0.6, s=40, c=results_df['chi_squared'], cmap='viridis_r')
        ax2.axhline(1.0, color='red', linestyle='--', alpha=0.7)
        ax2.set_xlabel('3.6μm Luminosity (10⁹ L☉)')
        ax2.set_ylabel('Velocity Ratio (Predicted/Observed)')
        ax2.set_title('UDT Performance vs Luminosity')
        ax2.set_xscale('log')
        ax2.grid(True, alpha=0.3)
        
        # 3. Best-fit radius distribution
        ax3.hist(results_df['best_radius_kpc'], bins=15, alpha=0.7, edgecolor='black')
        ax3.axvline(results_df['best_radius_kpc'].mean(), color='red', 
                   linestyle='--', label=f'Mean = {results_df["best_radius_kpc"].mean():.1f} kpc')
        ax3.set_xlabel('Best-fit Radius (kpc)')
        ax3.set_ylabel('Number of Galaxies')
        ax3.set_title('Distribution of Best-fit Radii')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # 4. Performance by galaxy type
        type_performance = results_df.groupby('type')['ratio_pred_obs'].agg(['mean', 'std']).reset_index()
        ax4.errorbar(type_performance['type'], type_performance['mean'], 
                    yerr=type_performance['std'], fmt='o-', capsize=5)
        ax4.axhline(1.0, color='red', linestyle='--', alpha=0.7)
        ax4.set_xlabel('Galaxy Type')
        ax4.set_ylabel('Mean Velocity Ratio')
        ax4.set_title('UDT Performance by Galaxy Type')
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        return fig
    
    def plot_chi_squared_analysis(self, results_df: pd.DataFrame,
                                save_path: Optional[str] = None) -> plt.Figure:
        """
        Comprehensive chi-squared analysis plots
        """
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
        
        # 1. Chi-squared histogram with fit quality regions
        chi2_values = results_df['chi_squared']
        ax1.hist(np.log10(chi2_values), bins=25, alpha=0.7, edgecolor='black', density=True)
        
        # Add quality regions
        ax1.axvspan(-2, np.log10(4), alpha=0.2, color='green', label='Excellent (χ² < 4)')
        ax1.axvspan(np.log10(4), np.log10(10), alpha=0.2, color='yellow', label='Good (4 ≤ χ² < 10)')
        ax1.axvspan(np.log10(10), np.log10(25), alpha=0.2, color='orange', label='Fair (10 ≤ χ² < 25)')
        ax1.axvspan(np.log10(25), 3, alpha=0.2, color='red', label='Poor (χ² ≥ 25)')
        
        ax1.set_xlabel('log₁₀(χ²)')
        ax1.set_ylabel('Probability Density')
        ax1.set_title('Distribution of χ² Values')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # 2. Cumulative chi-squared distribution
        sorted_chi2 = np.sort(chi2_values)
        cumulative = np.arange(1, len(sorted_chi2) + 1) / len(sorted_chi2)
        ax2.plot(sorted_chi2, cumulative * 100, 'b-', linewidth=2)
        ax2.axvline(4, color='red', linestyle='--', alpha=0.7, label='χ² = 4')
        ax2.axvline(10, color='orange', linestyle='--', alpha=0.7, label='χ² = 10')
        ax2.set_xlabel('χ²')
        ax2.set_ylabel('Cumulative Percentage')
        ax2.set_title('Cumulative χ² Distribution')
        ax2.set_xscale('log')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # 3. Chi-squared vs observed velocity
        ax3.scatter(results_df['v_observed'], results_df['chi_squared'], alpha=0.6, s=30)
        ax3.set_xlabel('Observed Velocity (km/s)')
        ax3.set_ylabel('χ²')
        ax3.set_title('Fit Quality vs Observed Velocity')
        ax3.set_yscale('log')
        ax3.grid(True, alpha=0.3)
        
        # 4. Fit quality pie chart
        quality_counts = results_df['fit_quality'].value_counts()
        colors = ['green', 'yellow', 'orange', 'red']
        ax4.pie(quality_counts.values, labels=quality_counts.index, autopct='%1.1f%%',
               colors=colors[:len(quality_counts)], startangle=90)
        ax4.set_title('Overall Fit Quality Distribution')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        return fig
    
    def plot_parameter_exploration(self, udt_framework: UDTFramework,
                                 test_galaxy: Dict = None,
                                 save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot UDT parameter space exploration
        """
        if test_galaxy is None:
            test_galaxy = {'name': 'Test', 'luminosity': 0.1, 'vflat': 50.0}
            
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
        
        M_star = udt_framework.luminosity_to_mass(test_galaxy['luminosity'])
        
        # 1. Velocity vs radius for different β values
        radii = np.linspace(1, 20, 50)
        beta_values = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
        
        for beta in beta_values:
            velocities = [udt_framework.rotation_velocity_udt(r, M_star, 
                                                           udt_framework.alpha, 
                                                           udt_framework.R0_galactic, 
                                                           beta) for r in radii]
            ax1.plot(radii, velocities, label=f'β = {beta}', linewidth=2)
        
        if 'vflat' in test_galaxy:
            ax1.axhline(test_galaxy['vflat'], color='red', linestyle='--', 
                       label=f'Observed ({test_galaxy["vflat"]} km/s)')
        
        ax1.set_xlabel('Radius (kpc)')
        ax1.set_ylabel('Velocity (km/s)')
        ax1.set_title('Velocity Profiles for Different β Values')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # 2. Velocity vs radius for different R₀ values
        R0_values = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0]
        
        for R0 in R0_values:
            velocities = [udt_framework.rotation_velocity_udt(r, M_star, 
                                                           udt_framework.alpha, 
                                                           R0, 
                                                           udt_framework.beta_galactic) for r in radii]
            ax2.plot(radii, velocities, label=f'R₀ = {R0} kpc', linewidth=2)
        
        if 'vflat' in test_galaxy:
            ax2.axhline(test_galaxy['vflat'], color='red', linestyle='--', 
                       label=f'Observed ({test_galaxy["vflat"]} km/s)')
        
        ax2.set_xlabel('Radius (kpc)')
        ax2.set_ylabel('Velocity (km/s)')
        ax2.set_title('Velocity Profiles for Different R₀ Values')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # 3. Distance dilation factor vs radius
        beta_values = [1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
        radii_extended = np.linspace(0.1, 50, 100)
        
        for beta in beta_values:
            D_factors = [udt_framework.distance_dilation_factor(r, 
                                                              udt_framework.R0_galactic, 
                                                              beta) for r in radii_extended]
            ax3.plot(radii_extended, D_factors, label=f'β = {beta}', linewidth=2)
        
        ax3.set_xlabel('Radius (kpc)')
        ax3.set_ylabel('Distance Dilation Factor D(r)')
        ax3.set_title('Distance Dilation vs Radius')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        ax3.set_xscale('log')
        
        # 4. Parameter space heatmap (β vs R₀)
        beta_range = np.linspace(2.0, 5.0, 20)
        R0_range = np.linspace(0.5, 5.0, 20)
        
        test_radius = 8.0  # kpc
        if 'vflat' in test_galaxy:
            target_velocity = test_galaxy['vflat']
            
            chi2_grid = np.zeros((len(beta_range), len(R0_range)))
            
            for i, beta in enumerate(beta_range):
                for j, R0 in enumerate(R0_range):
                    v_pred = udt_framework.rotation_velocity_udt(test_radius, M_star, 
                                                               udt_framework.alpha, R0, beta)
                    chi2_grid[i, j] = ((v_pred - target_velocity) / (target_velocity * 0.1))**2
            
            im = ax4.imshow(chi2_grid, extent=[R0_range.min(), R0_range.max(),
                                             beta_range.min(), beta_range.max()],
                          aspect='auto', origin='lower', cmap='viridis_r')
            
            # Mark current parameters
            ax4.scatter(udt_framework.R0_galactic, udt_framework.beta_galactic, 
                       color='red', s=100, marker='*', label='Current UDT')
            
            ax4.set_xlabel('R₀ (kpc)')
            ax4.set_ylabel('β')
            ax4.set_title('χ² Landscape in Parameter Space')
            ax4.legend()
            
            # Add colorbar
            cbar = plt.colorbar(im, ax=ax4)
            cbar.set_label('χ²')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        return fig
    
    def create_summary_figure(self, results_df: pd.DataFrame, 
                            validation_report: Dict,
                            save_path: Optional[str] = None) -> plt.Figure:
        """
        Create a comprehensive summary figure for publication
        """
        fig = plt.figure(figsize=(16, 12))
        
        # Create grid layout
        gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.3)
        
        # Main correlation plot (large)
        ax_main = fig.add_subplot(gs[0:2, 0:2])
        scatter = ax_main.scatter(results_df['v_observed'], results_df['v_predicted'], 
                                 c=np.log10(results_df['chi_squared']), 
                                 cmap='viridis_r', s=50, alpha=0.7)
        
        # Perfect correlation line
        min_v, max_v = 0, max(results_df['v_observed'].max(), results_df['v_predicted'].max())
        ax_main.plot([min_v, max_v], [min_v, max_v], 'r--', linewidth=2, alpha=0.8,
                    label='Perfect Agreement')
        
        ax_main.set_xlabel('SPARC Observed Velocity (km/s)', fontsize=14)
        ax_main.set_ylabel('UDT Predicted Velocity (km/s)', fontsize=14)
        ax_main.set_title('UDT Theory Validation Against SPARC Database', fontsize=16, fontweight='bold')
        ax_main.legend()
        ax_main.grid(True, alpha=0.3)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax_main)
        cbar.set_label('log₁₀(χ²)', fontsize=12)
        
        # Statistics text box
        stats_text = f"""Validation Statistics:
Total Galaxies: {validation_report['total_galaxies']}
Mean χ²: {validation_report['mean_chi_squared']:.1f}
Excellent Fits: {validation_report['excellent_fit_rate']:.1f}%
Good+ Fits: {validation_report['good_or_better_rate']:.1f}%
Mean Ratio: {validation_report['mean_velocity_ratio']:.2f}

Theory Parameters:
α = {validation_report['theoretical_parameters']['alpha']}
β = {validation_report['theoretical_parameters']['beta_galactic']}
R₀ = {validation_report['theoretical_parameters']['R0_galactic_kpc']} kpc"""
        
        ax_main.text(0.02, 0.98, stats_text, transform=ax_main.transAxes, 
                    fontsize=10, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Ratio histogram
        ax_hist = fig.add_subplot(gs[0, 2])
        ax_hist.hist(results_df['ratio_pred_obs'], bins=20, alpha=0.7, edgecolor='black')
        ax_hist.axvline(1.0, color='red', linestyle='--', label='Perfect match')
        ax_hist.axvline(results_df['ratio_pred_obs'].mean(), color='green', 
                       linestyle='-', label='Mean')
        ax_hist.set_xlabel('Velocity Ratio')
        ax_hist.set_ylabel('Count')
        ax_hist.set_title('Ratio Distribution')
        ax_hist.legend(fontsize=8)
        ax_hist.grid(True, alpha=0.3)
        
        # Chi-squared distribution
        ax_chi2 = fig.add_subplot(gs[1, 2])
        ax_chi2.hist(np.log10(results_df['chi_squared']), bins=20, alpha=0.7, edgecolor='black')
        ax_chi2.axvline(np.log10(4), color='red', linestyle='--', label='Good fit')
        ax_chi2.set_xlabel('log₁₀(χ²)')
        ax_chi2.set_ylabel('Count')
        ax_chi2.set_title('Fit Quality Distribution')
        ax_chi2.legend(fontsize=8)
        ax_chi2.grid(True, alpha=0.3)
        
        # Fit quality pie chart
        ax_pie = fig.add_subplot(gs[0, 3])
        quality_counts = results_df['fit_quality'].value_counts()
        colors = ['green', 'yellow', 'orange', 'red']
        ax_pie.pie(quality_counts.values, labels=quality_counts.index, autopct='%1.1f%%',
                  colors=colors[:len(quality_counts)], startangle=90)
        ax_pie.set_title('Fit Quality')
        
        # Performance by galaxy type
        ax_type = fig.add_subplot(gs[1, 3])
        type_performance = results_df.groupby('type')['ratio_pred_obs'].mean()
        ax_type.bar(type_performance.index, type_performance.values, alpha=0.7)
        ax_type.axhline(1.0, color='red', linestyle='--', alpha=0.7)
        ax_type.set_xlabel('Galaxy Type')
        ax_type.set_ylabel('Mean Ratio')
        ax_type.set_title('Performance by Type')
        ax_type.grid(True, alpha=0.3)
        
        # Sample rotation curves
        ax_curves = fig.add_subplot(gs[2, :])
        
        # Select a few representative galaxies
        sample_galaxies = results_df.nlargest(3, 'ratio_pred_obs').iloc[:3]
        
        udt = UDTFramework()
        sparc = SPARCAnalysis(udt)
        
        colors = ['blue', 'green', 'purple']
        radii = np.linspace(1, 20, 50)
        
        for i, (_, galaxy) in enumerate(sample_galaxies.iterrows()):
            galaxy_dict = {'name': galaxy['name'], 'luminosity': galaxy['luminosity']}
            _, velocities = sparc.predict_rotation_curve(galaxy_dict, radii)
            
            ax_curves.plot(radii, velocities, color=colors[i], linewidth=2,
                          label=f"{galaxy['name']} (UDT)")
            
            # Add observed point
            ax_curves.scatter(galaxy['best_radius_kpc'], galaxy['v_observed'], 
                            color=colors[i], s=80, marker='o', 
                            label=f"{galaxy['name']} (Obs)")
        
        ax_curves.set_xlabel('Radius (kpc)', fontsize=14)
        ax_curves.set_ylabel('Velocity (km/s)', fontsize=14)
        ax_curves.set_title('Sample UDT Rotation Curve Predictions', fontsize=14)
        ax_curves.legend(ncol=3, fontsize=10)
        ax_curves.grid(True, alpha=0.3)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        return fig

def create_publication_plots(results_df: pd.DataFrame, 
                           validation_report: Dict,
                           output_dir: str = "results/figures/") -> None:
    """
    Generate all publication-quality plots
    
    Args:
        results_df: Validation results DataFrame
        validation_report: Validation statistics
        output_dir: Directory to save plots
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    plotter = UDTPlotter(style='publication')
    
    print("🎨 Generating publication plots...")
    
    # 1. Summary figure
    fig1 = plotter.create_summary_figure(results_df, validation_report,
                                        f"{output_dir}/udt_validation_summary.png")
    print("  ✅ Summary figure saved")
    
    # 2. Scaling analysis
    fig2 = plotter.plot_velocity_scaling_analysis(results_df,
                                                 f"{output_dir}/udt_scaling_analysis.png")
    print("  ✅ Scaling analysis saved")
    
    # 3. Chi-squared analysis
    fig3 = plotter.plot_chi_squared_analysis(results_df,
                                           f"{output_dir}/udt_chi_squared_analysis.png")
    print("  ✅ Chi-squared analysis saved")
    
    # 4. Parameter exploration
    udt = UDTFramework()
    test_galaxy = {'name': 'Representative', 'luminosity': 0.1, 'vflat': 50.0}
    fig4 = plotter.plot_parameter_exploration(udt, test_galaxy,
                                            f"{output_dir}/udt_parameter_exploration.png")
    print("  ✅ Parameter exploration saved")
    
    plt.close('all')  # Clean up memory
    print(f"📊 All plots saved to {output_dir}")

if __name__ == "__main__":
    # Example usage
    print("UDT Plotting Tools - Example Usage")
    
    # This would typically be called after running validation
    # create_publication_plots(results_df, validation_report)