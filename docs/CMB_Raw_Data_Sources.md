# CMB Raw Data Sources for UDT Analysis

**Author: Charles Rotter**
**Date: 2025-01-17**

## Overview

This document catalogs available raw CMB data sources for UDT analysis, following the discovery of potential data contamination in the current Planck power spectrum analysis. The goal is to identify truly raw, model-independent CMB observables.

## Current Problem

The existing CMB analysis has critical issues:
- **Catastrophic fit quality**: χ²/dof = 35,563 (should be ~1)
- **Major discrepancies**: First peak prediction 29% off (284 vs 220)
- **Potential contamination**: Using processed power spectra instead of raw maps
- **Over-simplified physics**: Basic geometric scaling inadequate for CMB

## Available Raw Data Sources

### 1. Planck Mission Data

**Best Source**: ESA Planck Legacy Archive
- **URL**: https://pla.esac.esa.int/
- **Alternative**: NASA IRSA https://irsa.ipac.caltech.edu/data/Planck/

**Recommended Dataset**: 
- **File**: `COM_CMB_IQU-smica-nosz_2048_R3.00_full.fits`
- **Description**: SMICA CMB temperature map, DR3, with SZ sources projected out
- **Resolution**: Nside=2048 HEALPix format (~5 arcmin resolution)
- **Format**: FITS with HEALPix binary table extension
- **Size**: ~0.6 GB for intensity map, ~1.9 GB for full IQU
- **Coordinates**: Galactic, NESTED ordering
- **Units**: K_CMB

**Advantages**:
- Highest sensitivity and resolution
- Multiple cleaning methods (SMICA, COMMANDER, NILC, SEVEM)
- Well-documented foreground removal
- Standard reference for CMB community

**Potential Issues**:
- Foreground cleaning may introduce model dependencies
- Processing pipeline assumptions need verification

### 2. WMAP Mission Data

**Source**: NASA LAMBDA Archive
- **URL**: https://lambda.gsfc.nasa.gov/product/wmap/current/
- **Alternative**: https://wmap.gsfc.nasa.gov/resources/cmbimages.html

**Recommended Dataset**:
- **Description**: 9-year temperature and polarization maps
- **Resolution**: Nside=512 HEALPix format (~13 arcmin resolution)
- **Frequencies**: 5 bands (23-94 GHz)
- **Format**: FITS with HEALPix format

**Advantages**:
- Independent from Planck (different instrument, analysis pipeline)
- Multiple frequency bands for foreground assessment
- Long-established, well-validated dataset
- Simpler processing pipeline (less potential for contamination)

**Disadvantages**:
- Lower resolution than Planck
- Higher noise levels
- Limited sensitivity to high-l features

### 3. Atacama Cosmology Telescope (ACT) DR6

**Source**: NASA LAMBDA Archive
- **URL**: https://lambda.gsfc.nasa.gov/product/act/act_dr6.02/
- **Alternative**: NERSC Globus https://app.globus.org/file-manager?destination_id=53b2a147-ae9d-4bbf-9d18-3b46d133d4bb&destination_path=/act_dr6/

**Recommended Dataset**:
- **Description**: DR6 temperature and polarization maps
- **Resolution**: Arcminute resolution
- **Frequencies**: 98, 150, 220 GHz
- **Coverage**: 19,000 square degrees
- **Depth**: 10 μK arcmin median

**Advantages**:
- Very recent data (2017-2022)
- High resolution ground-based observations
- Independent systematic checks vs space missions
- Multiple frequency bands

**Disadvantages**:
- Partial sky coverage (not full-sky)
- Ground-based atmospheric contamination
- More complex noise properties

### 4. South Pole Telescope (SPT) Latest Data

**Source**: University of Chicago SPT Archive
- **URL**: https://pole.uchicago.edu/public/data/camphuis25/
- **Alternative**: NASA LAMBDA https://lambda.gsfc.nasa.gov/product/spt/

**Recommended Dataset**:
- **Description**: SPT-3G temperature and E-mode polarization maps (2024/2025)
- **Frequencies**: 90, 150, 220 GHz
- **Detectors**: 16,000-detector array
- **Quality**: "Unprecedentedly deep maps"

**Advantages**:
- Latest technology (SPT-3G camera)
- Highest sensitivity ground-based CMB measurements
- Independent systematic checks
- Fresh data with minimal legacy contamination

**Disadvantages**:
- Limited sky coverage
- Very recent (limited community validation)
- Ground-based complications

## Data Quality Assessment Criteria

### Raw Data Indicators
✅ **Clean**: Temperature maps in HEALPix format  
✅ **Clean**: Multi-frequency observations for foreground checks  
✅ **Clean**: Documented systematic error analysis  
⚠️ **Suspect**: Pre-computed power spectra  
❌ **Contaminated**: ΛCDM best-fit columns  
❌ **Contaminated**: Processed angular power spectra  

### Model Independence Tests
1. **Cross-dataset consistency**: Compare results across Planck/WMAP/ACT/SPT
2. **Frequency dependence**: Verify results are frequency-independent (CMB vs foregrounds)
3. **Systematic error budgets**: Ensure errors are properly characterized
4. **Null tests**: Check for unexpected correlations or biases

## Recommended Analysis Strategy

### Phase 1: Data Acquisition
1. **Download Planck SMICA map** (primary dataset)
2. **Download WMAP 9-year map** (independent check)
3. **Download ACT DR6 maps** (high-resolution check)
4. **Verify data integrity** and format consistency

### Phase 2: Cross-Validation
1. **Compare temperature maps** across different missions
2. **Check frequency consistency** within multi-band datasets
3. **Verify foreground cleaning** doesn't introduce UDT-relevant biases
4. **Assess systematic uncertainties** for UDT analysis

### Phase 3: UDT Analysis Development
1. **Develop proper CMB physics** for UDT framework
2. **Implement radiative transfer** with temporal geometry
3. **Test on simulated data** before applying to real observations
4. **Validate methodology** with known CMB features

## Success Metrics

### Minimum Requirements
- **Realistic fit quality**: χ²/dof ≈ 1-2
- **Acoustic peak accuracy**: <5% error on first peak position
- **Cross-dataset consistency**: Results agree across missions
- **Physical reasonableness**: UDT modifications are small corrections

### Validation Benchmarks
- **Compare with ΛCDM**: Should match within statistical uncertainties for standard features
- **Check new UDT predictions**: Look for specific signatures of temporal geometry
- **Systematic error budget**: Ensure UDT deviations exceed systematic uncertainties

## Next Steps

1. **Investigate Planck SMICA download**: Verify data access and format
2. **Download test dataset**: Start with smaller file to validate analysis pipeline
3. **Develop HEALPix analysis tools**: Python tools for map analysis
4. **Implement basic UDT CMB physics**: Start with simple temporal geometry effects

## Conclusion

Multiple independent CMB datasets are available for UDT analysis. The key is to use raw temperature maps rather than processed power spectra, and to cross-validate results across different missions. The catastrophic failure of the current analysis suggests we need to start completely fresh with proper raw data and realistic CMB physics.

Priority order for investigation:
1. **Planck SMICA** (highest quality, full sky)
2. **WMAP** (independent validation)
3. **ACT DR6** (high resolution check)
4. **SPT-3G** (latest technology validation)