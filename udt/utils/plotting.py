"""
Plotting utilities for UDT visualizations.

Provides consistent plotting functions for rotation curves,
Hubble diagrams, and other UDT analyses.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def setup_plot_style():
    """Set up consistent plot styling for UDT figures."""
    plt.style.use('seaborn-v0_8-darkgrid')
    plt.rcParams['figure.figsize'] = (10, 8)
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['axes.titlesize'] = 16
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['lines.linewidth'] = 2


def plot_rotation_curve(radius, velocity, velocity_error, 
                       v_predicted=None, galaxy_name="Galaxy",
                       save_path=None):
    """
    Plot galaxy rotation curve with UDT fit.
    
    Parameters
    ----------
    radius : array-like
        Galactocentric radius in kpc
    velocity : array-like
        Observed velocity in km/s
    velocity_error : array-like
        Velocity errors in km/s
    v_predicted : array-like, optional
        UDT model prediction
    galaxy_name : str
        Name of the galaxy
    save_path : str, optional
        Path to save figure
    """
    setup_plot_style()
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), 
                                   gridspec_kw={'height_ratios': [3, 1]})
    
    # Main plot
    ax1.errorbar(radius, velocity, yerr=velocity_error, 
                fmt='o', color='black', markersize=6, 
                label='Observed', capsize=3, alpha=0.7)
    
    if v_predicted is not None:
        ax1.plot(radius, v_predicted, 'r-', linewidth=2.5, 
                label='UDT Model', alpha=0.9)
    
    ax1.set_xlabel('Radius (kpc)')
    ax1.set_ylabel('Rotation Velocity (km/s)')
    ax1.set_title(f'{galaxy_name} - UDT Rotation Curve Fit')
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)
    
    # Residuals plot
    if v_predicted is not None:
        residuals = velocity - v_predicted
        ax2.errorbar(radius, residuals, yerr=velocity_error,
                    fmt='o', color='black', markersize=5, alpha=0.7)
        ax2.axhline(y=0, color='red', linestyle='--', alpha=0.5)
        ax2.set_xlabel('Radius (kpc)')
        ax2.set_ylabel('Residuals (km/s)')
        ax2.grid(True, alpha=0.3)
        
        # Add RMS to plot
        rms = np.sqrt(np.mean(residuals**2))
        ax2.text(0.95, 0.95, f'RMS = {rms:.1f} km/s',
                transform=ax2.transAxes, ha='right', va='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    
    return fig


def plot_hubble_diagram(z, m_obs, m_err=None, m_predicted=None,
                       title="UDT Hubble Diagram", save_path=None):
    """
    Plot supernova Hubble diagram with UDT fit.
    
    Parameters
    ----------
    z : array-like
        Redshifts
    m_obs : array-like
        Observed magnitudes
    m_err : array-like, optional
        Magnitude errors
    m_predicted : array-like, optional
        UDT model predictions
    title : str
        Plot title
    save_path : str, optional
        Path to save figure
    """
    setup_plot_style()
    
    fig = plt.figure(figsize=(12, 9))
    gs = gridspec.GridSpec(3, 1, height_ratios=[3, 1, 1])
    
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])
    
    # Main Hubble diagram
    if m_err is not None:
        ax1.errorbar(z, m_obs, yerr=m_err, fmt='o', 
                    color='black', markersize=5, alpha=0.6,
                    label='Observed', capsize=2)
    else:
        ax1.scatter(z, m_obs, color='black', s=30, alpha=0.6, label='Observed')
    
    if m_predicted is not None:
        # Sort for smooth line
        idx = np.argsort(z)
        ax1.plot(z[idx], m_predicted[idx], 'r-', linewidth=2.5,
                label='UDT Model', alpha=0.9)
    
    ax1.set_xscale('log')
    ax1.set_xlabel('Redshift z')
    ax1.set_ylabel('Apparent Magnitude')
    ax1.set_title(title)
    ax1.legend(loc='upper left')
    ax1.grid(True, alpha=0.3)
    
    # Residuals
    if m_predicted is not None:
        residuals = m_obs - m_predicted
        
        # Linear residuals
        ax2.scatter(z, residuals, color='black', s=20, alpha=0.6)
        ax2.axhline(y=0, color='red', linestyle='--', alpha=0.5)
        ax2.set_xscale('log')
        ax2.set_ylabel('Δm (mag)')
        ax2.grid(True, alpha=0.3)
        
        # Residuals vs magnitude
        ax3.scatter(m_obs, residuals, color='black', s=20, alpha=0.6)
        ax3.axhline(y=0, color='red', linestyle='--', alpha=0.5)
        ax3.set_xlabel('Apparent Magnitude')
        ax3.set_ylabel('Δm (mag)')
        ax3.grid(True, alpha=0.3)
        
        # Add statistics
        rms = np.sqrt(np.mean(residuals**2))
        ax2.text(0.95, 0.95, f'RMS = {rms:.3f} mag',
                transform=ax2.transAxes, ha='right', va='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    
    return fig


def plot_scale_comparison(galactic_scale, cosmic_scale, save_path=None):
    """
    Visualize the scale hierarchy between galactic and cosmic R₀.
    
    Parameters
    ----------
    galactic_scale : float
        R₀ for galaxies in kpc
    cosmic_scale : float
        R₀ for cosmology in Mpc
    save_path : str, optional
        Path to save figure
    """
    setup_plot_style()
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Convert to common units (kpc)
    cosmic_kpc = cosmic_scale * 1000
    
    # Create log scale visualization
    scales = [galactic_scale, cosmic_kpc]
    labels = [f'Galactic R₀\n{galactic_scale:.1f} kpc', 
              f'Cosmic R₀\n{cosmic_scale:.0f} Mpc']
    
    y_pos = [0.5, 0.5]
    colors = ['blue', 'red']
    
    for i, (scale, label, color) in enumerate(zip(scales, labels, colors)):
        ax.scatter(scale, y_pos[i], s=500, color=color, alpha=0.7, zorder=5)
        ax.text(scale, y_pos[i] + 0.1, label, ha='center', va='bottom',
               fontsize=14, weight='bold')
    
    # Add scale ratio
    ratio = cosmic_kpc / galactic_scale
    ax.text(np.sqrt(galactic_scale * cosmic_kpc), 0.5, 
           f'Scale Ratio\n{ratio:.0f}:1', 
           ha='center', va='center',
           bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5),
           fontsize=12)
    
    ax.set_xscale('log')
    ax.set_xlim(10, 1e7)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Characteristic Scale R₀ (kpc)', fontsize=14)
    ax.set_title('UDT Scale Hierarchy: Galactic vs Cosmic', fontsize=16)
    ax.grid(True, alpha=0.3, axis='x')
    ax.set_yticks([])
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    
    return fig