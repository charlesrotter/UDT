"""Universal Distance Dilation Theory (UDT) Package

A theoretical physics framework proposing a unified temporal geometry
for galactic and cosmological phenomena.
"""

__version__ = "0.1.0"
__author__ = "UDT Research Team"

from .core.temporal_geometry import (
    temporal_dilation_function,
    effective_light_speed,
    enhancement_factor
)

from .core.galactic_dynamics import (
    pure_temporal_velocity,
    fit_galaxy_rotation_curve
)

from .core.cosmology import (
    pure_temporal_magnitude,
    temporal_distance_modulus
)

__all__ = [
    'temporal_dilation_function',
    'effective_light_speed',
    'enhancement_factor',
    'pure_temporal_velocity',
    'fit_galaxy_rotation_curve',
    'pure_temporal_magnitude',
    'temporal_distance_modulus'
]