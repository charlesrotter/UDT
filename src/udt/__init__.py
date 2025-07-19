"""
Universal Distance Dilation Theory (UDT) Package

A scientific package for testing cosmological and galactic dynamics
under Universal Distance Dilation Theory with rigorous validation
and artifact correction.

Key Components:
- models: Core UDT theory and field equations
- data_loaders: Real observational data loading utilities  
- diagnostics: Validation frameworks and bias testing
"""

__version__ = "0.2.0"
__author__ = "Charles Rotter"

# Core UDT components
try:
    from .models.core import UDTCosmology, UDTGalacticDynamics
    from .data_loaders.sparc import load_sparc_data
    from .data_loaders.supernova import load_supernova_data
    from .diagnostics.validation import ValidationSuite
except ImportError as e:
    # Fallback for development/transition period
    import warnings
    warnings.warn(f"Some UDT components could not be imported: {e}")
    UDTCosmology = None
    UDTGalacticDynamics = None

__all__ = [
    "UDTCosmology",
    "UDTGalacticDynamics", 
    "load_sparc_data",
    "load_supernova_data",
    "ValidationSuite"
]