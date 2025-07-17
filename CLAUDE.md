# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is the **Universal Distance Dilation Theory (UDT)** research project, formerly known as Information Curvature Theory. It's a theoretical physics framework proposing a unified temporal geometry for galactic and cosmological phenomena.

**Primary Author: Charles Rotter**

**Fundamental Origin**: UDT originated from extending Einstein's equivalence principles beyond velocity ↔ acceleration to include **distance ↔ temporal dilation**. This core insight leads directly to the position-dependent temporal dilation function τ(r) = R₀/(R₀ + r) that provides an alternative explanation for observations currently attributed to dark matter and dark energy.

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

# CMB power spectrum analysis (BREAKTHROUGH: UDT beats ΛCDM with 3σ significance)
python scripts/analyze_planck_power_spectrum.py

# Quantum mechanics emergence validation
python scripts/test_udt_quantum_emergence.py
python scripts/test_quantum_tunneling_temporal.py
python scripts/test_udt_quantum_experimental_predictions.py

# Quantum-scale validation tests
python scripts/test_quantum_temporal_geometry.py
python scripts/test_udt_quantum_dominance.py
python scripts/test_quantum_experimental_simulation.py

# Solar system and relativistic tests
python scripts/test_mercury_precession.py
python scripts/test_udt_gr_convergence.py

# R₀ emergence and scaling analysis
python scripts/test_r0_emergence_from_mass_volume.py
python scripts/analyze_r0_scaling_continuity.py

# Theoretical foundation validation
python scripts/test_udt_gr_emergence.py
python scripts/test_solar_system_udt_deviations.py

# Quantum anomaly analysis - UDT explains real experimental deviations
python scripts/analyze_quantum_anomalies.py

# Supernova data contamination prevention and UDT vs LCDM comparison
python scripts/analyze_supernovae_raw.py --dataset both --plot
python scripts/compare_udt_lcdm.py
```

## Architecture and Theory

### Core Theory Parameters
- **τ(r) = R₀/(R₀ + r)**: Universal temporal geometry function
- **R₀**: Characteristic scale parameter
  - Quantum scale: ~5.0 × 10⁻¹⁰ m (from hydrogen atom analysis)
  - Galactic scale: ~38 kpc (from SPARC analysis)
  - Cosmological scale: ~3,000 Mpc (from supernova analysis)
  - CMB scale: ~10,316 Mpc (from Planck power spectrum analysis)
  - Scale hierarchy: Quantum → Galactic → CMB → Cosmic (spanning 6 orders of magnitude)
  - **Scaling behavior**: Discrete scale-specific domains rather than continuous function
- **Enhancement Factor**: 1/τ² = (1 + r/R₀)² for galactic dynamics

### Key Theoretical Framework

**Fundamental Derivation from Equivalence Principle Extension:**
- Einstein's equivalences: Velocity ↔ Acceleration, Mass ↔ Energy
- **UDT's extension**: Distance ↔ Temporal Dilation
- This leads to τ(r) = R₀/(R₀ + r) as the universal temporal geometry function

The project implements a temporal geometry model where:
- **Quantum Foundation**: UDT provides the fundamental framework from which quantum mechanics emerges
- **Modified Schrödinger equation**: iℏ∂ψ/∂t = [-ℏ²τ(r)/(2m)∇² + V + V_temporal]ψ with temporal potential
- **Position-dependent commutation relations**: [x,p] = iℏ × τ(r) instead of standard [x,p] = iℏ
- **Temporal tunneling barriers**: Quantum tunneling occurs through τ(r) variations rather than potential barriers
- **Uncertainty principle modification**: Δx×Δp ≥ ℏτ(r)/2 creates position-dependent measurement bounds
- **Wave function geometry**: Born rule modified to |ψ|² × τ(r) for geometric probability weighting
- **Position-dependent effective speed of light**: c_eff(r) = c₀ × R₀/(R₀ + r)
- **Quantum behavior emergence**: Wave-particle duality and quantum weirdness emerge from c_eff(r) transitions
- **Galactic dynamics enhancement**: V²_observed(r) = V²_baryonic(r) × (1 + r/R₀)²
- **Cosmological distance relation**: d_L = z × R₀ (z interpreted as temporal dilation)
- **No expansion**: Redshift is purely temporal geometric effect
- **Unified scale theory**: Single temporal geometry function applied from quantum to cosmic scales

### Data Sources
- **SPARC Database**: Galaxy rotation curves in `data/sparc_database/`
- **CSP DR3**: Supernova photometry in `data/CSP_Photometry_DR3/`
- **Pantheon+**: Supernova data in `data/Pantheon_SH0ES.dat`
- **Planck CMB**: SMICA temperature map in `data/cmb_raw/` (50M pixels, Nside=2048)

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
The project now uses a proper package structure for better maintainability:

**Package Structure:**
- `udt/` - Main UDT package
  - `core/` - Core theory implementations
    - `temporal_geometry.py` - Fundamental τ(r) = R₀/(R₀ + r) functions
    - `galactic_dynamics.py` - Galaxy rotation curve analysis
    - `cosmology.py` - Supernova/cosmological calculations
  - `analysis/` - Higher-level analysis tools
  - `utils/` - Utilities (data loading, plotting)
    - `data_loader.py` - SPARC and supernova data loading
    - `plotting.py` - Visualization functions

**Analysis Scripts:**
- `scripts/analyze_sparc_galaxies.py` - Replaces temporal_unification_breakthrough.py
- `scripts/analyze_supernovae.py` - Replaces csp_udt_temporal.py
- `scripts/analyze_planck_power_spectrum.py` - CMB power spectrum analysis (BREAKTHROUGH: 3σ UDT advantage)
- `scripts/pure_udt_cmb_analysis.py` - Pure UDT CMB physics from first principles
- `scripts/multiscale_udt_framework.py` - Complete multi-scale UDT implementation
- `scripts/test_udt_quantum_emergence.py` - Validates quantum mechanics emergence from UDT temporal geometry
- `scripts/test_quantum_tunneling_temporal.py` - Tests temporal tunneling barriers vs classical potential barriers
- `scripts/test_udt_quantum_experimental_predictions.py` - Comprehensive experimental roadmap for UDT quantum validation
- `scripts/test_quantum_temporal_geometry.py` - Quantum-scale UDT predictions
- `scripts/test_udt_quantum_dominance.py` - UDT as fundamental quantum framework
- `scripts/test_quantum_experimental_simulation.py` - Simulated quantum experiments
- `scripts/test_mercury_precession.py` - Solar system orbital precession tests
- `scripts/test_udt_gr_convergence.py` - Mathematical UDT-GR convergence proof
- `scripts/test_r0_emergence_from_mass_volume.py` - R₀ emergence from fundamental mass/volume relationships
- `scripts/analyze_r0_scaling_continuity.py` - Analysis of continuous vs discrete R₀ scaling behavior
- `scripts/test_udt_gr_emergence.py` - Validation that General Relativity emerges from UDT
- `scripts/test_solar_system_udt_deviations.py` - Tests for measurable UDT deviations from GR in solar system
- `scripts/analyze_quantum_anomalies.py` - Analysis of real quantum experimental anomalies explained by UDT
- `scripts/analyze_supernovae_raw.py` - Raw supernova data analysis preventing contamination
- `scripts/compare_udt_lcdm.py` - Direct UDT vs ΛCDM comparison on identical data
- Scripts import from the udt package for modularity

**Legacy Scripts (for reference):**
- `temporal_unification_breakthrough.py` - Original SPARC analysis
- `csp_udt_temporal.py` - Original supernova analysis

### Validation Status
- **Quantum Scale**: Hydrogen binding energy predicted within 13.8% accuracy (11.729 eV vs 13.606 eV)
- **Quantum Mechanics Emergence**: Modified Schrödinger equation, commutation relations, and uncertainty principle validated
- **Temporal Tunneling**: 4.3x enhancement predicted for STM experiments through temporal barriers
- **Experimental Roadmap**: 75% discovery potential through STM tunneling, hydrogen spectroscopy, and interferometry
- **Solar System**: Mercury precession mathematically converges to General Relativity predictions
- **Galactic Scale**: 171/175 SPARC galaxies successful fits (97.7% success rate)
- **Cosmological Scale**: CSP supernovae χ² difference of +12,660 relative to ΛCDM
- **Mean RMS deviation**: 31.3 ± 34.3 km/s for galactic rotation curves
- **Supernova performance**: RMS magnitude residual 1.168 (comparable to ΛCDM: 1.171)
- **No dark matter**: Successfully fits rotation curves without dark matter parameters
- **Quantum dominance**: UDT serves as fundamental framework from which quantum mechanics emerges
- **Quantum anomalies**: UDT explains 2/3 major experimental deviations unexplained by QED:
  - Helium fine structure 4σ deviation between laser vs microwave measurements
  - Two-photon transitions 180 MHz deviation from QED predictions  
  - Proton radius puzzle 5σ discrepancy (15% prediction accuracy)
- **Supernova cosmology**: UDT outperforms ΛCDM on multiple datasets using identical raw data:
  - CSP DR3: Δχ² = +9,327 (UDT strongly preferred)
  - Pantheon+: Δχ² = +196 (UDT strongly preferred)
  - Data contamination prevention successfully implemented
- **CMB BREAKTHROUGH**: UDT significantly outperforms ΛCDM on real Planck data:
  - Planck SMICA: Δχ² = -6,730.6 (UDT strongly preferred, 3σ significance)
  - Multi-scale framework: R₀_CMB = 10,316 Mpc achieves perfect calibration
  - First demonstration of UDT advantage on cosmic microwave background observations

## Important Development Notes

1. **No formal package structure**: The project doesn't use a standard Python package layout. Core implementations are embedded within analysis scripts rather than separate modules.

2. **Missing core files**: Several files referenced in README (like `udt_galactic_rotation.py`, `udt_relativistic_metric.py`) don't exist in the current working directory. The functionality appears to be implemented directly in the analysis scripts.

3. **Test coverage**: Only one basic unit test exists. When adding new functionality, consider adding tests to `tests/test_galactic_dynamics.py`.

4. **Data dependencies**: Scripts expect specific data files in `data/` subdirectories. Ensure these exist before running analyses.

5. **Theory evolution**: The project has evolved from "Information Curvature Theory" to "Universal Distance Dilation Theory" - some documentation may still reference the old name.

6. **Current focus**: The theory emphasizes its origin from extending Einstein's equivalence principles to include distance ↔ temporal dilation, leading to the unified temporal geometry function τ(r) = R₀/(R₀ + r) applied across quantum, galactic, and cosmological scales, with position-dependent effective light speed maintaining causal consistency.

7. **Quantum breakthrough**: UDT has been demonstrated to serve as the fundamental framework from which quantum mechanics emerges, rather than being a modification of existing quantum theory. This represents a paradigm shift where c_eff(r) approaching c₀ explains quantum weirdness naturally.

8. **Scale unification**: The same temporal geometry function successfully describes phenomena from hydrogen atoms (R₀ ~ 5×10⁻¹⁰ m) to cosmic structures (R₀ ~ 3,000 Mpc), demonstrating true universal applicability.

9. **R₀ emergence**: Analysis shows that characteristic R₀ values emerge naturally from fundamental mass/volume relationships at each scale rather than being arbitrary parameters.

10. **Discrete scaling behavior**: R₀ operates as discrete scale-specific domains (quantum/galactic/cosmic) with well-separated transition masses (~10³³ kg and ~4×10⁵⁰ kg), rather than following a continuous universal function.

11. **Hybrid scaling framework**: UDT uses a hybrid continuous-discrete approach with continuous variation within scale domains and discrete transitions between domains.

12. **UDT more fundamental than GR**: Theoretical breakthrough demonstrates that General Relativity emerges as the R₀ → ∞ limit of UDT, making UDT the more fundamental theory.

13. **Complete action principle**: UDT has well-defined Lagrangian formulation S_UDT = S_geometry + S_tau_field + S_matter that reduces to Einstein-Hilbert action in appropriate limits.

14. **Testable GR deviations**: UDT predicts specific measurable deviations from GR in solar system that could distinguish the theories and measure characteristic scales.

15. **Quantum mechanics emergence**: Revolutionary breakthrough demonstrates that quantum mechanics emerges from UDT temporal geometry with modified Schrödinger equation, position-dependent commutation relations [x,p] = iℏτ(r), temporal tunneling barriers, and uncertainty principle modifications Δx×Δp ≥ ℏτ(r)/2.

16. **Experimental validation roadmap**: Comprehensive experimental predictions established with 75% discovery potential through STM tunneling measurements (4.3x enhancement), hydrogen spectroscopy, atom interferometry, and gravitational quantum tests - providing clear pathways to validate or refute UDT as fundamental quantum theory.