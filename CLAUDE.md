# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is the **Universal Distance Dilation Theory (UDT)** research project, formerly known as Information Curvature Theory. It's a theoretical physics framework proposing a unified temporal geometry for galactic and cosmological phenomena. UDT introduces a position-dependent temporal dilation function that provides an alternative explanation for observations currently attributed to dark matter and dark energy.

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
- **τ(r) = R₀/(R₀ + r)**: Universal temporal geometry function
- **R₀**: Characteristic scale parameter
  - Galactic scale: ~38 kpc (from SPARC analysis)
  - Cosmological scale: ~3,000 Mpc (from supernova analysis)
  - Scale ratio: ~79,000:1
- **Enhancement Factor**: 1/τ² = (1 + r/R₀)² for galactic dynamics

### Key Theoretical Framework
The project implements a temporal geometry model where:
- Position-dependent effective speed of light: c_eff(r) = c₀ × R₀/(R₀ + r)
- Galactic dynamics enhancement: V²_observed(r) = V²_baryonic(r) × (1 + r/R₀)²
- Cosmological distance relation: d_L = z × R₀ (z interpreted as temporal dilation)
- No expansion - redshift is purely temporal geometric effect
- Single temporal geometry function applied across multiple scales

### Data Sources
- **SPARC Database**: Galaxy rotation curves in `data/sparc_database/`
- **CSP DR3**: Supernova photometry in `data/CSP_Photometry_DR3/`
- **Pantheon+**: Supernova data in `data/Pantheon_SH0ES.dat`

### Data Contamination Prevention
When performing model comparisons, it's critical to use raw observational data to avoid circular reasoning. Key guidelines:

**SPARC Raw Data (Use These):**
- V_obs(R): Direct rotation velocities from spectroscopy
- R: Galactocentric radius (verify distance source)
- err_V: Measurement uncertainties
- Spitzer 3.6μm surface brightness for stellar mass

**SPARC Contaminated Data (Avoid):**
- V_bar/V_obs ratios (pre-calculated with assumptions)
- Dark matter fractions (assumes baryonic models)
- Absolute magnitudes (distance-dependent)

**Pantheon+ Raw Data (Use These):**
- z_cmb: CMB frame redshift
- m_b_corrected: Peak apparent magnitude
- x1, c: SALT2 parameters

**Pantheon+ Contaminated Data (Avoid):**
- d_L: Luminosity distance (assumes cosmology)
- Hubble residuals (deviation from ΛCDM)
- mu_model: Model-predicted distance modulus

**Testing Protocol:**
1. Extract identical raw datasets for all models
2. Apply identical quality cuts and corrections
3. Document all processing assumptions
4. Compare models on exact same data points

See `docs/Data Contamination Prevention Guide.md` for complete details.

### Code Organization
The project uses standalone analysis scripts rather than a centralized module structure:
- `temporal_unification_breakthrough.py`: Main SPARC galaxy analysis using temporal framework
- `csp_udt_temporal.py`: CSP supernova analysis with pure temporal cosmology
- Individual analysis scripts at root level for specific investigations
- Results stored in `results/` with subdirectories for different analyses

### Validation Status
- SPARC galaxies: 171/175 successful fits (97.7% success rate)
- Mean RMS deviation: 31.3 ± 34.3 km/s for galactic rotation curves
- CSP supernovae: χ² difference of +12,660 relative to ΛCDM
- RMS magnitude residual: 1.168 (comparable to ΛCDM: 1.171)
- Successfully fits rotation curves without dark matter parameters

## Important Development Notes

1. **No formal package structure**: The project doesn't use a standard Python package layout. Core implementations are embedded within analysis scripts rather than separate modules.

2. **Missing core files**: Several files referenced in README (like `udt_galactic_rotation.py`, `udt_relativistic_metric.py`) don't exist in the current working directory. The functionality appears to be implemented directly in the analysis scripts.

3. **Test coverage**: Only one basic unit test exists. When adding new functionality, consider adding tests to `tests/test_galactic_dynamics.py`.

4. **Data dependencies**: Scripts expect specific data files in `data/` subdirectories. Ensure these exist before running analyses.

5. **Theory evolution**: The project has evolved from "Information Curvature Theory" to "Universal Distance Dilation Theory" - some documentation may still reference the old name.

6. **Current focus**: The theory now emphasizes a unified temporal geometry function τ(r) = R₀/(R₀ + r) applied across both galactic and cosmological scales, with position-dependent effective light speed maintaining causal consistency.