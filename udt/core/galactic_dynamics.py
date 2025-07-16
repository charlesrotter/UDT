"""
Galactic dynamics functions for UDT framework.

Implements rotation curve predictions and fitting functions
based on temporal geometry applied to galactic scales.
"""

import numpy as np
from scipy.optimize import minimize
from .temporal_geometry import temporal_dilation_function, enhancement_factor


def pure_temporal_velocity(r, R0_gal, V_scale):
    """
    Calculate rotation velocity using pure temporal geometry.
    
    Applies the temporal enhancement factor to base Newtonian dynamics
    to predict observed rotation curves without dark matter.
    
    Parameters
    ----------
    r : array-like
        Galactocentric radius in kpc
    R0_gal : float
        Galactic characteristic scale parameter in kpc
    V_scale : float
        Velocity scale parameter in km/s
        
    Returns
    -------
    v : array-like
        Predicted rotation velocity in km/s
    """
    # Ensure positive radius
    r = np.maximum(r, 1e-6)
    
    # Calculate temporal enhancement
    enhancement = enhancement_factor(r, R0_gal)
    
    # Base velocity profile (asymptotic to V_scale)
    # This represents the baryonic contribution
    v_base = V_scale * np.sqrt(r / (r + R0_gal/3))
    
    # Apply temporal enhancement
    return v_base * np.sqrt(enhancement)


def fit_galaxy_rotation_curve(radius, velocity, velocity_error, 
                            R0_bounds=(10, 100), V_scale_bounds=(50, 300)):
    """
    Fit UDT model to observed galaxy rotation curve.
    
    Parameters
    ----------
    radius : array-like
        Galactocentric radius in kpc
    velocity : array-like
        Observed rotation velocity in km/s
    velocity_error : array-like
        Velocity measurement errors in km/s
    R0_bounds : tuple, optional
        Bounds for R0_gal parameter search
    V_scale_bounds : tuple, optional
        Bounds for V_scale parameter search
        
    Returns
    -------
    dict
        Fit results including:
        - R0_gal: Best-fit characteristic scale
        - V_scale: Best-fit velocity scale
        - chi2: Chi-squared value
        - rms: RMS residual
        - success: Whether fit converged
    """
    # Objective function for optimization
    def objective(params):
        R0_gal, V_scale = params
        v_pred = pure_temporal_velocity(radius, R0_gal, V_scale)
        residuals = (velocity - v_pred) / velocity_error
        return np.sum(residuals**2)
    
    # Initial guess - use characteristic values
    R0_initial = np.median(radius)
    V_scale_initial = np.median(velocity)
    initial_params = [R0_initial, V_scale_initial]
    
    # Bounds
    bounds = [R0_bounds, V_scale_bounds]
    
    # Optimize
    result = minimize(objective, initial_params, bounds=bounds, method='L-BFGS-B')
    
    # Extract best-fit parameters
    R0_gal_fit, V_scale_fit = result.x
    
    # Calculate goodness of fit
    v_pred = pure_temporal_velocity(radius, R0_gal_fit, V_scale_fit)
    residuals = velocity - v_pred
    chi2 = result.fun
    rms = np.sqrt(np.mean(residuals**2))
    
    return {
        'R0_gal': R0_gal_fit,
        'V_scale': V_scale_fit,
        'chi2': chi2,
        'rms': rms,
        'success': result.success,
        'v_predicted': v_pred
    }


def calculate_mass_discrepancy(radius, v_observed, v_baryonic, R0_gal):
    """
    Calculate the mass discrepancy explained by temporal enhancement.
    
    Parameters
    ----------
    radius : array-like
        Galactocentric radius in kpc
    v_observed : array-like
        Observed rotation velocity in km/s
    v_baryonic : array-like
        Expected baryonic velocity in km/s
    R0_gal : float
        Galactic characteristic scale parameter
        
    Returns
    -------
    dict
        Contains:
        - enhancement: Temporal enhancement factor
        - mass_ratio: Apparent mass ratio (observed/baryonic)
        - v_temporal: Predicted velocity from temporal model
    """
    # Calculate enhancement
    enhancement = enhancement_factor(radius, R0_gal)
    
    # Mass ratio from velocity ratio
    mass_ratio = (v_observed / v_baryonic)**2
    
    # Predicted velocity with enhancement
    v_temporal = v_baryonic * np.sqrt(enhancement)
    
    return {
        'enhancement': enhancement,
        'mass_ratio': mass_ratio,
        'v_temporal': v_temporal,
        'temporal_mass_ratio': enhancement
    }