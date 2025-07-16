# Reference Implementations

This directory contains original implementations that served as the foundation for the current modular UDT framework.

## `original_implementations/`

### `csp_udt_temporal.py`
**Original breakthrough implementation of pure temporal cosmology**

- **Purpose**: Proof-of-concept implementation that first achieved excellent cosmological fits
- **Key Results**: 
  - R₀ = 3,000 Mpc (cosmic temporal scale)
  - RMS = 0.122 mag (excellent fit to sample data)
  - Δχ² = +339 vs ΛCDM (strongly preferred)
- **Status**: Reference only - functionality has been integrated into modular framework
- **Usage**: For validation, comparison, and historical documentation

**Important**: This script uses sample data and contains Unicode characters that may cause encoding issues. It has been preserved exactly as it was during the breakthrough analysis.

## Current Implementation

The functionality from these reference implementations has been integrated into the modular framework:

- **Core theory**: `udt/core/temporal_geometry.py`
- **Cosmological functions**: `udt/core/cosmology.py`
- **Analysis scripts**: `scripts/analyze_supernovae.py`
- **Galactic dynamics**: `scripts/analyze_sparc_galaxies.py`

## Usage Guidelines

1. **For Research**: Use reference implementations to understand the theoretical development
2. **For Production**: Use the modular framework in `udt/` and `scripts/`
3. **For Validation**: Compare new implementations against reference results
4. **For Documentation**: Historical record of the breakthrough implementation

## Key Differences

| Aspect | Reference Implementation | Current Framework |
|--------|-------------------------|-------------------|
| Data | Sample data generation | Real datasets (Pantheon+, SPARC) |
| Structure | Monolithic script | Modular packages |
| Units | Mixed kpc/Mpc handling | Consistent Mpc throughout |
| Encoding | Unicode characters | ASCII-compatible |
| Maintenance | Frozen reference | Actively maintained |

Last updated: 2025-01-16