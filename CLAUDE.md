# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is the **Universal Distance Dilation Theory (UDT)** research project, formerly known as Information Curvature Theory. It's a theoretical physics framework proposing a parameter-free geometric explanation for galactic rotation curves and cosmological phenomena without requiring dark matter or dark energy.

## Key Commands

### Setup
```bash
# Install all dependencies
pip install -r requirements.txt
```

### Running Tests
```bash
# Run the single unit test file
python tests/test_galactic_dynamics.py
```

### Main Analysis Scripts
```bash
# Run temporal unification analysis on SPARC galaxy data
python temporal_unification_breakthrough.py

# Run CSP supernova temporal analysis
python csp_udt_temporal.py
```

## Architecture and Theory

### Core Theory Parameters
- **β = 2.5**: Universal geometric exponent (from dimensional analysis)
- **α = 1/8 (0.125)**: Tesseract projection factor from 4D to 3D
- **R₀ = 1.5 kpc**: Characteristic galactic scale
- **τ(r) = R₀/(R₀ + r)**: Universal temporal geometry

### Key Theoretical Framework
The project implements a "pure temporal" universe model where:
- c = ∞ (infinite speed, instant causality)
- Distance ≡ Velocity ≡ Acceleration (unified through dilation)
- No expansion - redshift is purely temporal geometric effect
- Enhancement factor 1/τ² for galactic dynamics

### Data Sources
- **SPARC Database**: Galaxy rotation curves in `data/sparc_database/`
- **CSP DR3**: Supernova photometry in `data/CSP_Photometry_DR3/`
- **Pantheon+**: Supernova data in `data/Pantheon_SH0ES.dat`

### Code Organization
The project uses standalone analysis scripts rather than a centralized module structure:
- `temporal_unification_breakthrough.py`: Main SPARC galaxy analysis using temporal framework
- `csp_udt_temporal.py`: CSP supernova analysis with pure temporal cosmology
- Individual analysis scripts at root level for specific investigations
- Results stored in `results/` with subdirectories for different analyses

### Validation Status
- Phase 1 validation confirmed β = 2.512 ± 0.020 (target: 2.5)
- SPARC galaxies: χ² ~ 12.5 (some fits show χ² ~ 200.7 in newer runs)
- Successfully fits rotation curves without dark matter
- Active development with ongoing parameter optimization

## Important Development Notes

1. **No formal package structure**: The project doesn't use a standard Python package layout. Core implementations are embedded within analysis scripts rather than separate modules.

2. **Missing core files**: Several files referenced in README (like `udt_galactic_rotation.py`, `udt_relativistic_metric.py`) don't exist in the current working directory. The functionality appears to be implemented directly in the analysis scripts.

3. **Test coverage**: Only one basic unit test exists. When adding new functionality, consider adding tests to `tests/test_galactic_dynamics.py`.

4. **Data dependencies**: Scripts expect specific data files in `data/` subdirectories. Ensure these exist before running analyses.

5. **Theory evolution**: The project has evolved from "Information Curvature Theory" to "Universal Distance Dilation Theory" - some documentation may still reference the old name.