"""
BAO Data Loader

Functions for loading Baryon Acoustic Oscillation data.
Real observational data with model-independent analysis capabilities.
"""

import numpy as np
import pandas as pd
from pathlib import Path


def load_bao_data(data_file):
    """
    Load BAO (Baryon Acoustic Oscillation) compilation data.
    
    Parameters
    ----------
    data_file : str
        Path to BAO data file
        
    Returns
    -------
    list
        List of BAO measurements with redshift, parameter, value, and error
    """
    if not Path(data_file).exists():
        raise ValueError(f"BAO data file not found: {data_file}")
    
    try:
        with open(data_file, 'r') as f:
            lines = f.readlines()
        
        bao_data = []
        
        for line in lines:
            # Skip comments and empty lines
            if line.startswith('#') or not line.strip() or line.startswith('File written') or line.startswith('If you use'):
                continue
                
            parts = line.strip().split()
            if len(parts) >= 4:
                try:
                    z = float(parts[0])  # Redshift
                    value = float(parts[1])  # Measurement value
                    error = float(parts[2])  # Measurement error
                    parameter = parts[3]  # D_V/r_d, D_A/r_d, etc.
                    
                    # Additional info if available
                    survey = parts[6] if len(parts) > 6 else 'Unknown'
                    
                    if z > 0 and value > 0 and error > 0:
                        bao_data.append({
                            'z': z,
                            'parameter': parameter,
                            'value': value,
                            'error': error,
                            'survey': survey
                        })
                        
                except ValueError:
                    continue
        
        if len(bao_data) == 0:
            raise ValueError("No valid BAO data found")
        
        print(f"Loaded {len(bao_data)} BAO measurements from {data_file}")
        return bao_data
        
    except Exception as e:
        raise ValueError(f"Error loading BAO data: {e}")


def load_uncorrelated_bao_data(data_directory):
    """
    Load uncorrelated BAO dataset (uncorBAO.txt or similar).
    
    Parameters
    ----------
    data_directory : str
        Path to directory containing BAO data
        
    Returns
    -------
    list
        List of uncorrelated BAO measurements
    """
    data_path = Path(data_directory)
    
    # Look for BAO data files
    bao_patterns = ["uncorBAO.txt", "bao_compilation.dat", "*bao*.txt", "*BAO*.dat"]
    
    bao_files = []
    for pattern in bao_patterns:
        bao_files.extend(data_path.glob(pattern))
    
    if len(bao_files) == 0:
        raise ValueError(f"No BAO data files found in {data_directory}")
    
    # Use the first matching file
    return load_bao_data(bao_files[0])


def validate_bao_data_integrity(data_directory):
    """
    Validate BAO data integrity and check for model contamination.
    
    Parameters
    ----------
    data_directory : str
        Path to BAO data directory
        
    Returns
    -------
    dict
        Data integrity and contamination assessment
    """
    try:
        bao_data = load_uncorrelated_bao_data(data_directory)
        
        redshifts = [point['z'] for point in bao_data]
        parameters = [point['parameter'] for point in bao_data]
        values = [point['value'] for point in bao_data]
        errors = [point['error'] for point in bao_data]
        
        integrity_summary = {
            'status': 'valid',
            'n_measurements': len(bao_data),
            'redshift_range': [min(redshifts), max(redshifts)],
            'parameters': list(set(parameters)),
            'median_error_percent': np.median(np.array(errors) / np.array(values)) * 100,
            'data_source': 'Real BAO observations'
        }
        
        # Check for potential ΛCDM contamination
        contamination_indicators = []
        
        # Check if sound horizon is fixed to ΛCDM value
        dv_rd_values = [p['value'] for p in bao_data if 'D_V' in p['parameter']]
        if len(dv_rd_values) > 0:
            # If all D_V/r_d values assume r_d ≈ 147 Mpc, this indicates ΛCDM contamination
            implied_rd = np.mean(dv_rd_values) / 15  # Rough estimate
            if 140 < implied_rd < 155:
                contamination_indicators.append('Sound horizon may be ΛCDM-calibrated')
        
        if contamination_indicators:
            integrity_summary['contamination_warning'] = contamination_indicators
            integrity_summary['model_independent_analysis'] = 'recommended'
        
        return integrity_summary
        
    except Exception as e:
        return {
            'status': 'error',
            'error_message': str(e),
            'data_source': 'Unknown'
        }