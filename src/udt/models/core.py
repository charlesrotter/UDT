"""
Core UDT Models

Unified UDT classes for cosmological and galactic dynamics.
These are the validated implementations used throughout the analysis.
"""

import numpy as np


class UDTCosmology:
    """
    Universal Distance Dilation Theory cosmological model.
    
    Implements τ(r) = R₀/(R₀ + r) temporal geometry for 
    cosmological distance-redshift relationships.
    """
    
    def __init__(self, R0):
        """
        Initialize UDT cosmology model.
        
        Parameters
        ----------
        R0 : float
            Characteristic scale parameter (Mpc)
        """
        self.R0 = R0
        
    def temporal_dilation_function(self, r):
        """Calculate temporal dilation factor τ(r) = R₀/(R₀ + r)"""
        return self.R0 / (self.R0 + r)
    
    def redshift_from_distance(self, r):
        """Calculate redshift z = r/R₀"""
        return r / self.R0
        
    def distance_from_redshift(self, z):
        """Calculate distance r = z × R₀"""
        return z * self.R0
        
    def temporal_distance_modulus(self, z):
        """
        Calculate distance modulus in UDT temporal geometry.
        
        μ = 5 log₁₀(d_L/10pc) where d_L = (1+z) × r = (1+z)² × R₀
        """
        r = self.distance_from_redshift(z)
        d_L = (1 + z) * r  # Luminosity distance in UDT
        d_L_pc = d_L * 1e6  # Convert Mpc to pc
        return 5 * np.log10(d_L_pc / 10)


class UDTGalacticDynamics:
    """
    Universal Distance Dilation Theory galactic dynamics model.
    
    Implements enhancement factor 1/τ² = (1 + r/R₀)² for 
    galactic rotation curve fitting.
    """
    
    def __init__(self, R0):
        """
        Initialize UDT galactic dynamics model.
        
        Parameters
        ----------
        R0 : float
            Characteristic scale parameter (kpc)
        """
        self.R0 = R0
        
    def temporal_dilation_function(self, r):
        """Calculate temporal dilation factor τ(r) = R₀/(R₀ + r)"""
        return self.R0 / (self.R0 + r)
        
    def enhancement_factor(self, r):
        """Calculate enhancement factor 1/τ² = (1 + r/R₀)²"""
        tau = self.temporal_dilation_function(r)
        return 1 / (tau**2)
        
    def velocity_curve(self, r, M_stars, M_gas=0):
        """
        Calculate rotation velocity with UDT enhancement.
        
        v² = G(M_stars + M_gas) × enhancement_factor(r) / r
        
        Parameters
        ----------
        r : array-like
            Galactocentric radius (kpc)
        M_stars : float or array-like
            Stellar mass enclosed within radius r (solar masses)
        M_gas : float or array-like, optional
            Gas mass enclosed within radius r (solar masses)
            
        Returns
        -------
        v : array-like
            Rotation velocity (km/s)
        """
        G = 4.302e-6  # km²/s² kpc/Msun (gravitational constant)
        
        M_total = M_stars + M_gas
        enhancement = self.enhancement_factor(r)
        
        v_squared = G * M_total * enhancement / r
        return np.sqrt(v_squared)
        
    def fit_galaxy_data(self, r_obs, v_obs, v_err, M_stars, M_gas=0):
        """
        Fit galaxy rotation curve and return chi-squared.
        
        Parameters
        ----------
        r_obs : array-like
            Observed radii (kpc)
        v_obs : array-like  
            Observed velocities (km/s)
        v_err : array-like
            Velocity uncertainties (km/s)
        M_stars : array-like
            Stellar mass profile (solar masses)
        M_gas : array-like, optional
            Gas mass profile (solar masses)
            
        Returns
        -------
        chi2 : float
            Chi-squared value
        dof : int
            Degrees of freedom
        """
        v_pred = self.velocity_curve(r_obs, M_stars, M_gas)
        
        chi2 = np.sum(((v_obs - v_pred) / v_err)**2)
        dof = len(r_obs) - 1  # 1 parameter (R0) fitted
        
        return chi2, dof