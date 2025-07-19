# Quick Start Guide for Independent AI Verification

## Purpose
This guide enables an independent AI to quickly assess the UDT project for scientific integrity and verify major claims within 30-60 minutes.

## Step 1: Essential File Checks (5 minutes)

### Check Data Authenticity
```bash
# Count SPARC database files (should be ~175)
ls C:/UDT/data/sparc_database/*.dat | wc -l

# Examine sample galaxy data
head -n 20 C:/UDT/data/sparc_database/NGC3198_rotmod.dat
```

**Expected**: Real observational data with columns for radius, velocity, errors
**Red Flag**: Missing files, synthetic markers, or perfect mathematical patterns

### Check for Synthetic Data Functions
```bash
grep -r "_create_sample_galaxies\|generate.*data\|synthetic" C:/UDT/scripts/ C:/UDT/mathematical_development/
```

**Expected**: No results (archived code should be in archive/ directory only)
**Red Flag**: Active synthetic data generation in analysis scripts

## Step 2: Core Implementation Verification (10 minutes)

### Verify Enhancement Factor Implementation
**File**: `C:/UDT/udt/core/galactic_dynamics.py`

**Look for function**: `enhancement_factor(r, R0_gal)` or similar

**Should implement**: F(τ) = 1 + α × 3(1-τ)/(τ²(3-2τ)) where τ = R₀/(R₀ + r)

**Red Flag**: Hard-coded values, missing function, or implementation not matching formula

### Check Data Loading
**File**: `C:/UDT/udt/utils/data_loader.py`

**Look for function**: `load_sparc_database(data_dir)`

**Should do**: Load actual .dat files from filesystem

**Red Flag**: Functions that create synthetic data when files missing

## Step 3: Quick Analysis Test (10 minutes)

### Test SPARC Analysis on Small Sample
```bash
cd C:/UDT
python scripts/analyze_sparc_galaxies.py --max-galaxies 3 --output-dir quick_test
```

**Expected Output**:
- Analysis completes without errors
- R₀ values in range 10-100 kpc
- RMS residuals 3-15 km/s (not perfect)
- Some variation between galaxies

**Red Flags**:
- Script fails or uses fallback synthetic data
- All galaxies give identical results
- Perfect fits (RMS ≈ 0)
- Unrealistic parameter values

## Step 4: Verify Major Claims (10 minutes)

### Check CLAUDE.md Claims Against Evidence

**CLAIM**: "Real SPARC database (175 galaxies)"
**VERIFY**: File count from Step 1 should be ~175
**STATUS**: ✅ Pass / ❌ Fail

**CLAIM**: "Perfect 5.0/5 score on SPARC data"  
**VERIFY**: Run comprehensive validation if time permits
**STATUS**: ✅ Pass / ❌ Fail / ⏸ Skip

**CLAIM**: "LIGO timing prediction matches observations"
**VERIFY**: Check if quantum_validation/udt_ligo_final_analysis.py exists and uses real parameters
**STATUS**: ✅ Pass / ❌ Fail / ⏸ Skip

**CLAIM**: "Archived synthetic data fallbacks"
**VERIFY**: Check archive/deprecated_fallbacks/ exists with temporal_unification_breakthrough.py
**STATUS**: ✅ Pass / ❌ Fail

## Step 5: Scientific Integrity Assessment (5 minutes)

### High-Priority Red Flags
Check for these critical issues:

1. **Synthetic Data in Active Code**
   ```bash
   grep -n "create.*galaxies\|synthetic" C:/UDT/scripts/analyze_sparc_galaxies.py
   ```
   **Expected**: No matches
   
2. **Hard-Coded Results**
   ```bash
   grep -n "return.*success.*True" C:/UDT/scripts/analyze_sparc_galaxies.py
   ```
   **Expected**: No hard-coded success values

3. **Circular Validation**
   Look for: Code that uses UDT formulas to generate data, then fits UDT to that data
   **Expected**: Analysis fits UDT to independently obtained observational data

## Quick Assessment Decision Tree

### ✅ GREEN LIGHT (Proceed with confidence)
- Real SPARC data files present and loaded correctly
- Enhancement factor implementation matches mathematical formula  
- Test analysis produces reasonable, varied results
- No evidence of synthetic data in active analysis
- Major claims supported by file/code evidence

### ⚠️ YELLOW LIGHT (Proceed with caution)
- Minor implementation issues or missing features
- Some documentation inconsistencies 
- Test analysis works but with concerns about methodology
- Claims partially supported but need verification

### 🛑 RED LIGHT (Major concerns identified)
- Missing or synthetic SPARC data
- Enhancement factor implementation incorrect or missing
- Test analysis fails or produces suspicious results  
- Evidence of hard-coded results or circular validation
- Major claims unsupported by actual evidence

## 30-Second Assessment Protocol

If extremely limited time, check only these essentials:

1. **Data files exist**: `ls C:/UDT/data/sparc_database/*.dat | wc -l` returns ~175
2. **No active synthetic data**: `grep -r "_create_sample_galaxies" C:/UDT/scripts/` returns nothing
3. **Core function exists**: `grep -n "def enhancement_factor" C:/UDT/udt/core/galactic_dynamics.py` finds function
4. **Analysis script exists**: `ls C:/UDT/scripts/analyze_sparc_galaxies.py` file exists

**If all 4 pass**: Likely legitimate research (GREEN LIGHT)
**If any fail**: Requires investigation (YELLOW/RED LIGHT)

## Sample Quick Report Template

```
# 30-Minute UDT Verification Report

## Overall Assessment: [GREEN/YELLOW/RED LIGHT]

## Key Checks:
- SPARC data files: [✅ Pass / ❌ Fail] - [File count]
- Synthetic data search: [✅ Pass / ❌ Fail] - [Findings]  
- Enhancement factor: [✅ Pass / ❌ Fail] - [Implementation status]
- Test analysis: [✅ Pass / ❌ Fail] - [Output summary]
- Major claims: [✅ Pass / ❌ Fail] - [Verification results]

## Critical Issues Found: [None / List issues]

## Recommendations: [Continue verification / Investigate concerns / Major red flags identified]

## Verification Time: [X minutes]
```

## Next Steps Based on Assessment

### GREEN LIGHT → Extended Verification
- Proceed to full verification checklist  
- Run comprehensive statistical validation
- Examine theoretical framework in detail
- Assess manuscript readiness

### YELLOW LIGHT → Targeted Investigation  
- Focus on specific concerns identified
- Clarify ambiguous claims or implementations
- Request additional documentation
- Limited additional verification

### RED LIGHT → Critical Investigation
- Document all major issues found
- Examine for scientific misconduct patterns
- Recommend external review if serious concerns
- Do not endorse results without resolution

## Important Notes

- This is a **screening tool**, not a comprehensive review
- GREEN LIGHT means "worth deeper investigation," not "validated science"  
- RED LIGHT indicates serious concerns requiring expert review
- Focus on **objective verification** of claims vs evidence
- Report findings clearly and constructively