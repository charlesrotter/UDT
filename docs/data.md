# Data Guide

This guide explains how to obtain and verify the observational datasets used in UDT validation.

## Data Sources

All UDT analyses use **real observational data only**. No synthetic data generation is used in any validation.

### SPARC Galaxy Database

**Source**: [SPARC Database](http://astroweb.cwru.edu/SPARC/)  
**Citation**: Lelli, McGaugh & Schombert (2016)  
**Contents**: Rotation curves for 175 galaxies

#### Download Instructions
```bash
# Download SPARC database
wget http://astroweb.cwru.edu/SPARC/MassModels_Lelli2016c.mrt
wget http://astroweb.cwru.edu/SPARC/RotationCurves_Lelli2016c.tar.gz

# Extract rotation curves
tar -xzf RotationCurves_Lelli2016c.tar.gz

# Organize data
mkdir -p data/sparc_database
mv MassModels_Lelli2016c.mrt data/sparc_database/
mv RotationCurves_Lelli2016c/* data/sparc_database/
```

#### Data Format
- **File format**: Machine-readable tables (.mrt)
- **Columns**: Radius (kpc), Velocity (km/s), Velocity Error (km/s)
- **Galaxy count**: 175 galaxies
- **Quality**: High-quality HI and Hα observations

### Supernova Data

#### Pantheon+ Compilation
**Source**: [Pantheon+ SH0ES](https://github.com/PantheonPlusSH0ES/DataRelease)  
**Citation**: Brout et al. (2022)  
**Contents**: 1701 Type Ia supernovae

```bash
# Download Pantheon+ data
wget https://github.com/PantheonPlusSH0ES/DataRelease/raw/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat

# Organize
mkdir -p data/supernova
mv PantheonSH0ES.dat data/supernova/
```

#### Carnegie Supernova Project DR3
**Source**: [CSP DR3](https://csp.obs.carnegiescience.edu/data)  
**Citation**: Krisciunas et al. (2017)  
**Contents**: Near-infrared supernova photometry

```bash
# Download CSP DR3 data
mkdir -p data/CSP_Photometry_DR3
# Manual download required from CSP website
```

### CMB Data

#### Planck Power Spectrum
**Source**: [Planck Legacy Archive](http://pla.esac.esa.int/pla/)  
**Citation**: Planck Collaboration (2020)  
**Contents**: CMB temperature and polarization power spectra

```bash
# Download Planck power spectrum data
mkdir -p data/planck_cmb
# Download COM_PowerSpect_CMB files from PLA
```

### BAO Data

#### Uncorrelated BAO Compilation
**Source**: [BAO Compilations](https://github.com/sfschen/BAOfit)  
**Citation**: Chen et al. (2019)  
**Contents**: Baryon acoustic oscillation measurements

```bash
# Download BAO compilation
wget https://raw.githubusercontent.com/sfschen/BAOfit/master/data/uncorBAO.txt
mv uncorBAO.txt data/bao/
```

## Data Integrity Verification

### SHA-256 Manifests

All data files are verified using SHA-256 checksums:

```bash
# Create data integrity manifest
udt-validate integrity --data-dir data/ --create-manifest

# Verify data integrity
udt-validate integrity --data-dir data/ --verify-manifest
```

### Expected Data Structure
```
data/
├── sparc_database/
│   ├── MassModels_Lelli2016c.mrt
│   ├── NGC0024.mrt
│   ├── NGC0055.mrt
│   └── ... (175 galaxy files)
├── supernova/
│   ├── PantheonSH0ES.dat
│   └── CSP_Photometry_DR3/
├── planck_cmb/
│   └── COM_PowerSpect_CMB_*.dat
├── bao/
│   └── uncorBAO.txt
└── manifest_sha256.txt
```

## Data Quality Assurance

### Artifact Correction

**IMPORTANT**: All supernova and BAO data requires artifact correction to remove ΛCDM processing contamination.

#### Automatic Correction
The UDT validation suite automatically applies validated artifact correction:

```python
from src.udt.diagnostics.artifact_correction import UnbiasedArtifactCorrection

corrector = UnbiasedArtifactCorrection()
corrected_data = corrector.apply_conservative_bias_correction(raw_data)
```

#### Manual Verification
You can verify the correction is working:

```bash
# Run bias testing on correction methods
python mathematical_development/artifact_correction_bias_testing.py
```

### Synthetic Data Detection

The framework automatically scans for synthetic data contamination:

```python
from src.udt.diagnostics.integrity import DataIntegrityChecker

checker = DataIntegrityChecker()
authenticity = checker.verify_data_authenticity(data_array, "SPARC rotation curves")
print(f"Authenticity score: {authenticity['authenticity_score']}%")
```

## Data Limitations and Caveats

### SPARC Database
- **Sample selection**: Primarily nearby, gas-rich galaxies
- **Systematic uncertainties**: Distance measurements, inclination corrections
- **Resolution limits**: Inner galaxy regions may be beam-smeared

### Supernova Data
- **ΛCDM contamination**: K-corrections and calibrations assume ΛCDM
- **Selection effects**: Malmquist bias affects high-redshift samples
- **Dust corrections**: May introduce systematic uncertainties

### CMB Data
- **Foreground subtraction**: Galactic contamination removal
- **Beam effects**: Instrument response corrections
- **Systematic uncertainties**: Calibration and pointing errors

### BAO Data
- **Model dependence**: Sound horizon calibration assumes ΛCDM
- **Non-linear corrections**: May introduce model bias
- **Scale dependence**: Assumptions about cosmological parameters

## Data Access and Redistribution

### Licensing
- **SPARC**: Open access for scientific use
- **Pantheon+**: Open access with citation requirement
- **Planck**: ESA data policy applies
- **BAO compilations**: Open access

### Attribution Requirements
When using this data, cite the original sources:

```bibtex
@article{sparc_database,
  title={SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves},
  author={Lelli, Federico and McGaugh, Stacy S and Schombert, James M},
  journal={The Astronomical Journal},
  year={2016}
}

@article{pantheon_plus,
  title={The Pantheon+ Analysis: Cosmological Constraints},
  author={Brout, Dillon and others},
  journal={The Astrophysical Journal},
  year={2022}
}
```

## Troubleshooting

### Common Issues

#### Data Download Failures
- Check network connectivity
- Verify URLs haven't changed
- Use alternative mirrors if available

#### Integrity Check Failures
- Re-download affected files
- Check for partial downloads
- Verify file permissions

#### Artifact Correction Issues
- Ensure all required dependencies installed
- Check for consistent data formats
- Validate input data ranges

### Getting Help

If you encounter data-related issues:

1. Check the [troubleshooting guide](troubleshooting.md)
2. Search [existing issues](https://github.com/charlesrotter/UDT/issues)
3. Open a new issue with detailed error information

## Data Update Protocol

### When to Update Data
- New releases of observational datasets
- Discovery of data quality issues
- Updates to artifact correction methods

### Update Procedure
1. Download new data files
2. Update SHA-256 manifest
3. Verify artifact correction still works
4. Update documentation
5. Re-run validation suite

### Version Control
All data updates are tracked in the repository:
- Data files: SHA-256 manifest changes
- Methodology: Code version control
- Results: Archived in separate branches