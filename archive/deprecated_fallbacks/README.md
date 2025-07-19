# Deprecated Fallback Code Archive

This directory contains code that was found to contain deceptive synthetic data fallbacks that could mislead future analysis.

## Archived Files

### `temporal_unification_breakthrough.py`
- **Issue**: Contains `_create_sample_galaxies()` function that generates synthetic data matching UDT predictions
- **Problem**: Lines 317-325 create artificial rotation curves using the exact UDT formula being tested
- **Risk**: Could cause false validation if fallback code is accidentally used
- **Archived**: 2025-07-19
- **Reason**: Prevents contamination of future validation efforts

## Code Fragment Identified as Problematic

```python
# Lines 317-325: Artificial data created to match UDT predictions
tau_r = R0_gal / (R0_gal + radius)  # τ(r) = R₀/(R₀ + r)  
enhancement = 1 / (tau_r**2)        # 1/τ² enhancement
v_temporal = v_base * np.sqrt(enhancement)  # UDT formula used to CREATE the data
```

## Replacement Strategy

The validated analysis now uses:
- **Real data**: Actual SPARC rotation curves from observational database
- **Proper UDT package**: `udt/core/galactic_dynamics.py` with clean implementations
- **No fallbacks**: Production code fails gracefully rather than using synthetic data

## Prevention Measures

1. **Archived location**: Code moved out of active development paths
2. **Documentation**: Clear marking of synthetic data generation
3. **Code review**: Future implementations checked for similar patterns
4. **Validation protocols**: Always verify data sources in analysis

## Historical Note

This archival is part of a contamination audit that discovered synthetic data fallbacks were being confused with validated results. The main UDT validation has been confirmed to use real observational data.