"""
Core temporal geometry functions for UDT framework.

Implements the fundamental temporal dilation function τ(r) = R₀/(R₀ + r)
and derived quantities that form the basis of UDT theory.
"""

import numpy as np


def temporal_dilation_function(r, R0):
    """
    Calculate the temporal dilation factor τ(r).
    
    This is the fundamental function of UDT theory:
    τ(r) = R₀/(R₀ + r)
    
    Parameters
    ----------
    r : float or array-like
        Distance from origin
    R0 : float
        Characteristic scale parameter
        
    Returns
    -------
    tau : float or array-like
        Temporal dilation factor at distance r
    """
    return R0 / (R0 + r)


def effective_light_speed(r, R0, c0=299792.458):
    """
    Calculate position-dependent effective speed of light.
    
    c_eff(r) = c₀ × τ(r) = c₀ × R₀/(R₀ + r)
    
    Parameters
    ----------
    r : float or array-like
        Distance from origin
    R0 : float
        Characteristic scale parameter
    c0 : float, optional
        Base speed of light (km/s), default is standard value
        
    Returns
    -------
    c_eff : float or array-like
        Effective light speed at distance r
    """
    tau = temporal_dilation_function(r, R0)
    return c0 * tau


def enhancement_factor(r, R0):
    """
    Calculate the enhancement factor for galactic dynamics.
    
    Enhancement = 1/τ² = (1 + r/R₀)²
    
    This factor enhances gravitational effects in galactic rotation curves.
    
    Parameters
    ----------
    r : float or array-like
        Distance from galactic center
    R0 : float
        Characteristic scale parameter (galactic scale)
        
    Returns
    -------
    enhancement : float or array-like
        Enhancement factor at distance r
    """
    tau = temporal_dilation_function(r, R0)
    return 1 / (tau**2)


def temporal_redshift(r, R0):
    """
    Calculate redshift from temporal geometry.
    
    In UDT, redshift is purely a temporal dilation effect:
    z = r/R₀
    
    Parameters
    ----------
    r : float or array-like
        Distance from origin
    R0 : float
        Characteristic scale parameter (cosmological scale)
        
    Returns
    -------
    z : float or array-like
        Redshift at distance r
    """
    return r / R0


def distance_from_redshift(z, R0):
    """
    Calculate distance from observed redshift.
    
    Inverse of temporal_redshift:
    r = z × R₀
    
    Parameters
    ----------
    z : float or array-like
        Observed redshift
    R0 : float
        Characteristic scale parameter (cosmological scale)
        
    Returns
    -------
    r : float or array-like
        Distance corresponding to redshift z
    """
    return z * R0