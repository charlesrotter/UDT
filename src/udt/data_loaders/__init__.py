"""
UDT Data Loaders Package

Real observational data loading utilities with built-in integrity checking
and artifact correction. No synthetic data generation.

MANDATORY VALIDATION REQUIRED:
All data loading and analysis must pass through validation gate system.
"""

from .sparc import load_sparc_data, load_sparc_galaxy
from .supernova import load_supernova_data, load_csp_data, load_pantheon_data
from .cmb import load_planck_data
from .bao import load_bao_data
from .ligo import load_ligo_data, load_gw150914_parameters
from .muon_g2 import load_muon_g2_data, predict_udt_muon_anomaly

# Import validation gate system
from ..diagnostics.mandatory_validation_gate import (
    validate_before_analysis, 
    check_validation_status,
    requires_validation,
    ValidationError
)

__all__ = [
    "load_sparc_data", "load_sparc_galaxy",
    "load_supernova_data", "load_csp_data", "load_pantheon_data", 
    "load_planck_data", "load_bao_data",
    "load_ligo_data", "load_gw150914_parameters",
    "load_muon_g2_data", "predict_udt_muon_anomaly",
    "validate_before_analysis", "check_validation_status", 
    "requires_validation", "ValidationError"
]