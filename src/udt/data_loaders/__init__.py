"""
UDT Data Loaders Package

Real observational data loading utilities with built-in integrity checking
and artifact correction. No synthetic data generation.
"""

from .sparc import load_sparc_data, load_sparc_galaxy
from .supernova import load_supernova_data, load_csp_data, load_pantheon_data
from .cmb import load_planck_data
from .bao import load_bao_data
from .ligo import load_ligo_data, load_gw150914_parameters

__all__ = [
    "load_sparc_data", "load_sparc_galaxy",
    "load_supernova_data", "load_csp_data", "load_pantheon_data", 
    "load_planck_data", "load_bao_data",
    "load_ligo_data", "load_gw150914_parameters"
]