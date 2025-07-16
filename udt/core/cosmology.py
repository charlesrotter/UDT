"""
Cosmological functions for UDT framework.

Implements supernova magnitude predictions and cosmological
distance calculations based on temporal geometry.
"""

import numpy as np
from scipy.optimize import minimize
from .temporal_geometry import temporal_redshift, distance_from_redshift


def pure_temporal_magnitude(z, R0, M_B):
    """
    Calculate apparent magnitude using pure temporal geometry.
    
    In UDT, the luminosity distance is directly related to redshift:
    d_L = z × R₀
    
    Parameters
    ----------
    z : float or array-like
        Observed redshift
    R0 : float
        Cosmological characteristic scale parameter in Mpc
    M_B : float
        Absolute B-band magnitude
        
    Returns
    -------
    m_B : float or array-like
        Apparent B-band magnitude
    """
    # Luminosity distance in kpc
    d_L_kpc = distance_from_redshift(z, R0)
    
    # Convert to parsecs for distance modulus
    d_L_pc = d_L_kpc * 1000
    
    # Distance modulus: m - M = 5 log10(d_L/10 pc)
    return M_B + 5 * np.log10(d_L_pc / 10)


def temporal_distance_modulus(z, R0):
    """
    Calculate distance modulus from redshift in UDT framework.
    
    μ = 5 log10(d_L/10 pc) where d_L = z × R₀
    
    Parameters
    ----------
    z : float or array-like
        Observed redshift
    R0 : float
        Cosmological characteristic scale parameter in Mpc
        
    Returns
    -------
    mu : float or array-like
        Distance modulus
    """
    d_L_kpc = distance_from_redshift(z, R0)
    d_L_pc = d_L_kpc * 1000
    return 5 * np.log10(d_L_pc / 10)


def fit_supernova_hubble_diagram(z_obs, m_obs, m_err, 
                               R0_bounds=(1000, 10000), 
                               M_B_bounds=(-20, -17)):
    """
    Fit UDT model to supernova Hubble diagram.
    
    Parameters
    ----------
    z_obs : array-like
        Observed redshifts
    m_obs : array-like
        Observed apparent magnitudes
    m_err : array-like
        Magnitude errors
    R0_bounds : tuple, optional
        Bounds for R0 parameter search (in Mpc)
    M_B_bounds : tuple, optional
        Bounds for absolute magnitude
        
    Returns
    -------
    dict
        Fit results including:
        - R0: Best-fit characteristic scale in Mpc
        - M_B: Best-fit absolute magnitude
        - chi2: Chi-squared value
        - rms: RMS magnitude residual
        - success: Whether fit converged
    """
    # Objective function
    def objective(params):
        R0, M_B = params
        m_pred = pure_temporal_magnitude(z_obs, R0, M_B)
        residuals = (m_obs - m_pred) / m_err
        return np.sum(residuals**2)
    
    # Initial guess
    initial_R0 = 3000  # Mpc, based on previous analyses
    initial_M_B = -19  # Typical SNIa absolute magnitude
    initial_params = [initial_R0, initial_M_B]
    
    # Bounds
    bounds = [R0_bounds, M_B_bounds]
    
    # Optimize
    result = minimize(objective, initial_params, bounds=bounds, method='L-BFGS-B')
    
    # Extract best-fit parameters
    R0_fit, M_B_fit = result.x
    
    # Calculate goodness of fit
    m_pred = pure_temporal_magnitude(z_obs, R0_fit, M_B_fit)
    residuals = m_obs - m_pred
    chi2 = result.fun
    rms = np.sqrt(np.mean(residuals**2))
    
    return {
        'R0': R0_fit,
        'M_B': M_B_fit,
        'chi2': chi2,
        'rms': rms,
        'success': result.success,
        'm_predicted': m_pred,
        'residuals': residuals
    }


def calculate_hubble_parameter(R0):
    """
    Calculate effective Hubble parameter from UDT framework.
    
    In UDT, the apparent Hubble parameter emerges from temporal geometry:
    H_eff = c / R₀
    
    Parameters
    ----------
    R0 : float
        Characteristic scale parameter in Mpc
        
    Returns
    -------
    H0 : float
        Effective Hubble parameter in km/s/Mpc
    """
    c_km_s = 299792.458  # Speed of light in km/s
    return c_km_s / R0


def compare_with_standard_cosmology(z, R0_udt, M_B_udt, 
                                  H0_standard=70, Om_standard=0.3):
    """
    Compare UDT predictions with standard ΛCDM cosmology.
    
    Parameters
    ----------
    z : array-like
        Redshift values for comparison
    R0_udt : float
        UDT characteristic scale in Mpc
    M_B_udt : float
        UDT absolute magnitude
    H0_standard : float, optional
        Standard cosmology Hubble constant
    Om_standard : float, optional
        Standard cosmology matter density
        
    Returns
    -------
    dict
        Comparison results including magnitudes and residuals
    """
    # UDT predictions
    m_udt = pure_temporal_magnitude(z, R0_udt, M_B_udt)
    
    # Standard cosmology (simplified ΛCDM)
    # This is approximate - full calculation would use integration
    c = 299792.458  # km/s
    d_L_standard = c * z / H0_standard * (1 + 0.5 * (1 - Om_standard) * z)
    d_L_pc = d_L_standard * 1e6  # Convert Mpc to pc
    m_standard = M_B_udt + 5 * np.log10(d_L_pc / 10)
    
    return {
        'z': z,
        'm_udt': m_udt,
        'm_standard': m_standard,
        'delta_m': m_udt - m_standard,
        'd_L_udt': distance_from_redshift(z, R0_udt),
        'd_L_standard': d_L_standard
    }