"""
CMB Data Loader

Functions for loading Planck CMB power spectrum data.
Real observational data only.
"""

import numpy as np
import pandas as pd
from pathlib import Path


def load_planck_data(data_directory, dataset='SMICA'):
    """
    Load Planck CMB power spectrum data.
    
    Parameters
    ----------
    data_directory : str
        Path to Planck data directory
    dataset : str
        Planck dataset to load ('SMICA', 'temperature', etc.)
        
    Returns
    -------
    dict
        CMB power spectrum data with ell, C_ell, and errors
    """
    data_path = Path(data_directory)
    
    # Look for Planck power spectrum files
    planck_patterns = [
        f"*{dataset}*.dat",
        f"*{dataset}*.txt", 
        "COM_PowerSpect*.dat",
        "planck*.dat"
    ]
    
    planck_files = []
    for pattern in planck_patterns:
        planck_files.extend(data_path.glob(pattern))
    
    if len(planck_files) == 0:
        raise ValueError(f"No Planck {dataset} files found in {data_directory}")
    
    # Use the first matching file
    planck_file = planck_files[0]
    
    try:
        # Load power spectrum data
        with open(planck_file, 'r') as f:
            lines = f.readlines()
        
        # Skip header lines
        data_lines = [line for line in lines if not line.startswith('#') and line.strip()]
        
        ell = []
        c_ell = []
        c_ell_error = []
        
        for line in data_lines:
            parts = line.strip().split()
            if len(parts) >= 3:
                try:
                    l = int(parts[0])
                    cl = float(parts[1])
                    cl_err = float(parts[2]) if len(parts) > 2 else cl * 0.01  # 1% default error
                    
                    if l > 0 and cl > 0:  # Valid power spectrum point
                        ell.append(l)
                        c_ell.append(cl)
                        c_ell_error.append(cl_err)
                        
                except ValueError:
                    continue
        
        if len(ell) == 0:
            raise ValueError(f"No valid power spectrum data in {planck_file}")
        
        cmb_data = {
            'ell': np.array(ell),
            'c_ell': np.array(c_ell),
            'c_ell_error': np.array(c_ell_error),
            'dataset': dataset,
            'source_file': str(planck_file),
            'n_points': len(ell)
        }
        
        print(f"Loaded {len(ell)} CMB power spectrum points from {planck_file.name}")
        return cmb_data
        
    except Exception as e:
        raise ValueError(f"Error loading Planck data from {planck_file}: {e}")


def validate_cmb_data_integrity(data_directory):
    """
    Validate CMB data integrity.
    
    Parameters
    ----------
    data_directory : str
        Path to CMB data directory
        
    Returns
    -------
    dict
        Data integrity summary
    """
    try:
        cmb_data = load_planck_data(data_directory)
        
        ell = cmb_data['ell']
        c_ell = cmb_data['c_ell']
        
        integrity_summary = {
            'status': 'valid',
            'n_points': len(ell),
            'ell_range': [int(ell.min()), int(ell.max())],
            'c_ell_range': [float(c_ell.min()), float(c_ell.max())],
            'dataset': cmb_data['dataset'],
            'data_source': 'Real Planck observations'
        }
        
        # Check for reasonable power spectrum structure
        if ell.max() < 100:
            integrity_summary['warning'] = 'Limited ell range - may be partial data'
        
        return integrity_summary
        
    except Exception as e:
        return {
            'status': 'error',
            'error_message': str(e),
            'data_source': 'Unknown'
        }