# Migration Guide: UDT v0.1 → v0.2

This guide explains how to update code for the new UDT package structure.

## Overview of Changes

UDT v0.2 introduces a professional package structure following scientific Python standards:

### Old Structure (v0.1)
```
UDT/
├── udt/                    # Old package location
├── scripts/                # Mixed with package code
├── mathematical_development/
├── CLAUDE.md              # All documentation in one file
```

### New Structure (v0.2)
```
UDT/
├── src/udt/               # New package location
│   ├── models/            # Core theory
│   ├── data_loaders/      # Data utilities
│   └── diagnostics/       # Validation & testing
├── docs/                  # Structured documentation
├── examples/              # Jupyter notebooks
├── scripts/               # Clean analysis scripts
├── quarantine/            # Isolated problematic methods
```

## Import Changes

### Core Models
```python
# OLD (v0.1)
from udt.core.temporal_geometry import temporal_dilation_function
from udt.core.galactic_dynamics import fit_galaxy_rotation_curve

# NEW (v0.2)
from src.udt.models.core import UDTCosmology, UDTGalacticDynamics
# OR if package is installed
from udt.models.core import UDTCosmology, UDTGalacticDynamics
```

### Data Loading
```python
# OLD (v0.1)
from udt.utils.data_loader import load_sparc_database

# NEW (v0.2)
from src.udt.data_loaders.sparc import load_sparc_data
# OR if package is installed
from udt.data_loaders.sparc import load_sparc_data
```

### Validation and Testing
```python
# OLD (v0.1)
# No unified validation framework

# NEW (v0.2)
from src.udt.diagnostics.validation import ValidationSuite
from src.udt.diagnostics.artifact_correction import UnbiasedArtifactCorrection
```

## API Changes

### Galaxy Analysis
```python
# OLD (v0.1) - function-based
from udt.core.galactic_dynamics import fit_galaxy_rotation_curve
result = fit_galaxy_rotation_curve(radius, velocity, velocity_error)

# NEW (v0.2) - class-based
from udt.models.core import UDTGalacticDynamics
model = UDTGalacticDynamics(R0=25.0)
chi2, dof = model.fit_galaxy_data(radius, velocity, velocity_error, mass_profile)
```

### Cosmological Analysis
```python
# OLD (v0.1) - function-based
from udt.core.cosmology import temporal_distance_modulus
mu = temporal_distance_modulus(redshift, R0)

# NEW (v0.2) - class-based
from udt.models.core import UDTCosmology
cosmo = UDTCosmology(R0=4000.0)
mu = cosmo.temporal_distance_modulus(redshift)
```

## Installation Changes

### Development Installation
```bash
# OLD (v0.1)
pip install -e .

# NEW (v0.2) - improved dependency handling
pip install --prefer-binary numpy scipy matplotlib
pip install --no-build-isolation -e .[dev]
```

### Package Structure
```bash
# OLD (v0.1) - packages found at root
[build-system]
[tool.setuptools.packages.find]
where = ["."]

# NEW (v0.2) - packages in src/
[tool.setuptools.packages.find]
where = ["src"]
[tool.setuptools.package-dir]
"" = "src"
```

## CLI Changes

### New Command Line Interface
```bash
# NEW (v0.2) - unified CLI
udt-validate validate --sparc-data data/sparc_database --max-galaxies 50
udt-validate integrity --data-dir data/ --verify-manifest
```

## Documentation Changes

### File Locations
- **OLD**: All in `CLAUDE.md`
- **NEW**: Structured in `docs/`
  - `docs/index.md` - Main documentation
  - `docs/theory.md` - Theoretical framework
  - `docs/data.md` - Data guide
  - `docs/quarantine.md` - Quarantine system

### Examples
- **OLD**: No structured examples
- **NEW**: Interactive Jupyter notebooks in `examples/`

## Script Migration

### Existing Scripts
Most existing scripts will continue to work with minor modifications:

```python
# Add at top of script for transition period
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))
```

### Recommended Updates
1. **Update imports** to use new package structure
2. **Use class-based APIs** instead of standalone functions
3. **Add data integrity checks** using new diagnostics tools
4. **Apply artifact correction** for supernova/cosmological analyses

## Testing Changes

### Test Structure
```bash
# OLD (v0.1) - minimal testing
tests/test_galactic_dynamics.py

# NEW (v0.2) - comprehensive testing
tests/
├── test_models/
├── test_data_loaders/
├── test_diagnostics/
└── test_integration/
```

### CI/CD Improvements
- Data integrity verification with SHA-256 manifests
- Bias testing for artifact correction methods
- Synthetic data contamination detection
- Package import validation

## Data Changes

### Data Organization
```bash
# OLD (v0.1) - mixed data locations
data/
├── various files...

# NEW (v0.2) - structured with integrity checking
data/
├── sample/              # Small files for CI
├── sparc_database/      # Galaxy data
├── supernova/           # SN data  
├── cmb_planck/         # CMB data
└── manifest_sha256.txt  # Integrity manifest
```

### Data Integrity
- **NEW**: SHA-256 manifests for all data files
- **NEW**: Automated integrity verification in CI
- **NEW**: Synthetic data detection and prevention

## Quarantine System

### New Safety Features
- **Problematic methods isolated** in `quarantine/` directory
- **Clear warnings** prevent accidental use
- **Validated alternatives** for all quarantined methods

```bash
quarantine/
├── biased_methods/      # Methodologically flawed
├── synthetic_data/      # Artificial data use
└── deprecated_approaches/ # Outdated methods
```

## Migration Checklist

For existing code:

- [ ] Update import statements to new package structure
- [ ] Replace function calls with class-based APIs
- [ ] Add data integrity verification
- [ ] Apply artifact correction for cosmological data
- [ ] Update test scripts to use new validation framework
- [ ] Check for quarantined method usage
- [ ] Update documentation references

For new development:

- [ ] Use `src/udt/` package structure
- [ ] Follow class-based API patterns
- [ ] Include mandatory artifact correction
- [ ] Add comprehensive bias testing
- [ ] Use data integrity manifests
- [ ] Write examples as Jupyter notebooks

## Support

If you encounter issues during migration:

1. Check this migration guide
2. Review `docs/` for updated documentation
3. Run `udt-validate integrity --verify-manifest` to check data
4. Test imports: `python -c "from udt.models.core import UDTCosmology"`
5. Open GitHub issue with specific error details

## Timeline

- **v0.1**: Legacy structure (deprecated)
- **v0.2**: New professional structure (current)
- **Future**: Legacy imports may be removed in v0.3+

We recommend migrating to the new structure as soon as possible to benefit from improved validation, testing, and scientific integrity features.