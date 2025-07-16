# Universal Distance Dilation Theory (UDT)

A theoretical physics framework proposing a unified temporal geometry for galactic and cosmological phenomena. UDT introduces a position-dependent temporal dilation function that provides an alternative explanation for observations currently attributed to dark matter and dark energy.

## Overview

The Universal Distance Dilation Theory presents a geometric framework based on temporal dilation:
- Single temporal geometry function τ(r) = R₀/(R₀ + r) applied across multiple scales
- Analysis of 175 SPARC galaxies with 171 successful fits (97.7% success rate)
- Comparison with standard cosmology using 133 Type Ia supernovae data
- Position-dependent effective speed of light maintaining causal consistency
- Extension of equivalence principles to include temporal-spatial relationships

## Theoretical Framework

### Core Equations

**Temporal Geometry Function:**
```
τ(r) = R₀/(R₀ + r)
```

**Galactic Dynamics Enhancement:**
```
Enhancement Factor = 1/τ² = (1 + r/R₀)²
V²_observed(r) = V²_baryonic(r) × (1 + r/R₀)²
```

**Cosmological Distance Relation:**
```
d_L = z × R₀
z = r/R₀ (interpreted as temporal dilation)
```

**Effective Speed of Light:**
```
c_eff(r) = c₀ × R₀/(R₀ + r)
```

### Characteristic Scales
- Galactic scale: R₀ ~ 38 kpc (from SPARC analysis)
- Cosmological scale: R₀ ~ 3,000 Mpc (from supernova analysis)
- Scale ratio: ~79,000:1

## Observational Analysis

### Galactic Rotation Curves (SPARC Database)
- Dataset: 175 disk galaxies with Spitzer photometry
- Successful fits: 171/175 galaxies (97.7%)
- Mean RMS deviation: 31.3 ± 34.3 km/s
- Method: Pure temporal geometry without dark matter parameters

### Type Ia Supernovae (Carnegie Supernova Project)
- Dataset: 133 supernovae (redshift range: 0.0046 - 0.0835)
- Comparative performance: χ² difference of +12,660 relative to ΛCDM
- RMS magnitude residual: 1.168 (comparable to ΛCDM: 1.171)
- Parameters: R₀ = 3,000 Mpc, M_B = -17.6

## Implementation

### Installation
```bash
# Install dependencies
pip install -r requirements.txt
```

### Running Analyses
```bash
# Test package imports and functionality
python test_imports.py

# Galactic rotation curve analysis (new organized script)
python scripts/analyze_sparc_galaxies.py --plot

# Supernova cosmology analysis (new organized script)
python scripts/analyze_supernovae.py --dataset csp --plot

# Legacy analysis scripts (for reference)
python temporal_unification_breakthrough.py
python csp_udt_temporal.py
```

### Core Functions
```python
# Import from organized package structure
from udt.core.temporal_geometry import temporal_dilation_function, enhancement_factor
from udt.core.galactic_dynamics import pure_temporal_velocity, fit_galaxy_rotation_curve
from udt.core.cosmology import pure_temporal_magnitude, fit_supernova_hubble_diagram

# Basic usage examples
tau = temporal_dilation_function(r=50, R0=38)  # kpc
enhancement = enhancement_factor(r=50, R0=38)
v_predicted = pure_temporal_velocity(radius, R0_gal=38, V_scale=200)
m_predicted = pure_temporal_magnitude(z=0.05, R0=3000, M_B=-19)
```

## Data Sources
- **SPARC Database**: `data/sparc_database/` - Galactic rotation curves
- **CSP DR3**: `data/CSP_Photometry_DR3/` - Supernova photometry
- **Pantheon+**: `data/Pantheon_SH0ES.dat` - Extended supernova catalog

## Data Contamination Prevention

When testing UDT against standard models, it's critical to avoid circular reasoning by using data that has been pre-processed through competing models. See `docs/Data Contamination Prevention Guide.md` for detailed guidance.

### Quick Reference

#### SPARC Database - Use Raw Data:
- ✅ **V_obs(R)**: Direct rotation velocities from HI/Hα spectroscopy
- ✅ **R**: Galactocentric radius (check distance method)
- ✅ **err_V**: Observational uncertainties
- ✅ **Spitzer 3.6μm**: Raw infrared photometry
- ❌ **Avoid**: V_bar/V_obs ratios, dark matter fractions, absolute magnitudes

#### Pantheon+ Database - Use Raw Data:
- ✅ **z_cmb**: CMB frame redshift (kinematic corrections only)
- ✅ **m_b_corrected**: Peak apparent magnitude
- ✅ **x1, c**: SALT2 light curve parameters
- ❌ **Avoid**: d_L (luminosity distance), Hubble residuals, mu_model

### Key Principles
1. Always trace back to direct observational measurements
2. Check assumptions in any processed data
3. Use identical raw datasets for model comparisons
4. Verify geometric/kinematic distances over Hubble flow distances
5. Document all data processing steps

## Project Structure
```
├── udt/                                 # Main UDT package
│   ├── core/                           # Core theory implementations
│   │   ├── temporal_geometry.py        # Fundamental τ(r) functions
│   │   ├── galactic_dynamics.py        # Galaxy rotation curves
│   │   └── cosmology.py                # Supernova/cosmological
│   ├── analysis/                       # Analysis tools
│   └── utils/                          # Utilities (data, plotting)
├── scripts/                            # New organized analysis scripts
│   ├── analyze_sparc_galaxies.py       # SPARC galaxy analysis
│   └── analyze_supernovae.py           # Supernova analysis
├── temporal_unification_breakthrough.py # Legacy SPARC analysis
├── csp_udt_temporal.py                 # Legacy supernova analysis
├── test_imports.py                     # Package functionality test
├── data/                               # Observational datasets
├── results/                            # Analysis outputs
└── docs/                               # Documentation
```

## Key Findings

### Theoretical Implications
- Proposes geometric origin for galactic rotation curve anomalies
- Suggests alternative interpretation of cosmological redshift
- Introduces position-dependent effective light speed
- Extends equivalence principle to temporal-spatial domain

### Observational Correspondence
- Successfully fits majority of SPARC galaxy rotation curves
- Provides comparable fit quality to standard cosmology for nearby supernovae
- Maintains consistency with special and general relativity in local limit
- Preserves causality through effective light speed structure

## Comparison with Standard Model

| Aspect | ΛCDM | UDT |
|--------|------|-----|
| Free Parameters | Multiple | Single scale parameter |
| Dark Matter | Required | Not required |
| Dark Energy | Λ constant | Geometric effect |
| Redshift Interpretation | Expansion | Temporal dilation |
| Light Speed | Constant | Position-dependent |

## Current Limitations and Future Work

### Limitations
- Analysis limited to z < 0.1 for supernovae
- CMB and BAO analyses not yet performed
- Gravitational lensing predictions require development
- Quantum mechanical implications unexplored

### Proposed Extensions
1. High-redshift supernova analysis
2. CMB power spectrum predictions
3. Gravitational lensing formulation
4. Laboratory tests of c_eff variations

## Documentation
- `docs/UDT Framework.md`: Detailed theoretical development
- `CLAUDE.md`: Project development guidelines
- Analysis results in `results/` directory

## References

See individual analysis scripts for detailed citations to observational datasets and related theoretical work.

---

*Note: This is an active research project exploring alternative cosmological models. Results should be interpreted within the context of ongoing theoretical physics research.*
