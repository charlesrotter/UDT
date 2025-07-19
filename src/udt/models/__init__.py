"""
UDT Models Package

Core theoretical models and field equations for Universal Distance Dilation Theory.
All models are validated against real observational data with mandatory artifact correction.
"""

from .core import UDTCosmology, UDTGalacticDynamics
from .field_equations import solve_udt_field_equations

__all__ = ["UDTCosmology", "UDTGalacticDynamics", "solve_udt_field_equations"]