"""
SPARC Galaxy Data Loader

Validated functions for loading real SPARC galaxy rotation curve data.
No synthetic data generation - uses only real observational data.
"""

import os
import glob
import numpy as np
import pandas as pd
from pathlib import Path


def load_sparc_galaxy(galaxy_file):
    """
    Load individual SPARC galaxy rotation curve data.
    
    Parameters
    ----------
    galaxy_file : str
        Path to galaxy data file
        
    Returns
    -------
    dict or None
        Galaxy data including name, radius, velocity, errors
        Returns None if loading fails
    """
    try:
        # Read the file
        with open(galaxy_file, 'r') as f:
            lines = f.readlines()
        
        # Extract galaxy name from filename
        galaxy_name = os.path.basename(galaxy_file).replace('.mrt', '').replace('.txt', '')
        
        # Parse data lines (skip comments)
        data_lines = [line for line in lines if not line.startswith('#') and line.strip()]
        
        radius = []
        velocity = []
        velocity_error = []
        
        for line in data_lines:
            parts = line.strip().split()
            if len(parts) >= 3:
                try:
                    r = float(parts[0])  # Radius (kpc)
                    v = float(parts[1])  # Velocity (km/s)
                    v_err = float(parts[2])  # Velocity error (km/s)
                    
                    # Data quality checks - only use positive values
                    if r > 0 and v > 0 and v_err > 0:
                        radius.append(r)
                        velocity.append(v)
                        velocity_error.append(v_err)
                except ValueError:
                    continue
        
        if len(radius) > 0:
            return {
                'name': galaxy_name,
                'radius': np.array(radius),
                'velocity': np.array(velocity),
                'velocity_error': np.array(velocity_error),
                'n_points': len(radius),
                'r_max': np.max(radius)
            }
        else:
            return None
            
    except Exception as e:
        print(f"Error loading {galaxy_file}: {e}")
        return None


def load_sparc_data(data_directory, max_galaxies=None):
    """
    Load SPARC galaxy database from directory.
    
    Parameters
    ----------
    data_directory : str
        Path to SPARC data directory
    max_galaxies : int, optional
        Maximum number of galaxies to load
        
    Returns
    -------
    list
        List of galaxy data dictionaries
        
    Raises
    ------
    ValueError
        If data directory doesn't exist or no valid data files found
    """
    data_path = Path(data_directory)
    if not data_path.exists():
        raise ValueError(f"SPARC data directory not found: {data_directory}")
    
    # Look for galaxy data files
    patterns = ['*.mrt', '*.txt', '*.dat']
    galaxy_files = []
    
    for pattern in patterns:
        galaxy_files.extend(data_path.glob(pattern))
    
    if len(galaxy_files) == 0:
        raise ValueError(f"No SPARC data files found in {data_directory}")
    
    # Load galaxies
    galaxies = []
    loaded_count = 0
    
    for galaxy_file in sorted(galaxy_files):
        if max_galaxies and loaded_count >= max_galaxies:
            break
            
        galaxy_data = load_sparc_galaxy(galaxy_file)
        if galaxy_data is not None:
            galaxies.append(galaxy_data)
            loaded_count += 1
    
    if len(galaxies) == 0:
        raise ValueError(f"No valid galaxy data loaded from {data_directory}")
    
    print(f"Loaded {len(galaxies)} SPARC galaxies from {data_directory}")
    return galaxies


def get_sparc_galaxy_names(data_directory):
    """
    Get list of available galaxy names in SPARC database.
    
    Parameters
    ----------
    data_directory : str
        Path to SPARC data directory
        
    Returns
    -------
    list
        List of galaxy names
    """
    try:
        galaxies = load_sparc_data(data_directory)
        return [g['name'] for g in galaxies]
    except Exception as e:
        print(f"Error loading galaxy names: {e}")
        return []


def validate_sparc_data_integrity(data_directory):
    """
    Validate SPARC data integrity and return summary statistics.
    
    Parameters
    ----------
    data_directory : str
        Path to SPARC data directory
        
    Returns
    -------
    dict
        Data integrity summary
    """
    try:
        galaxies = load_sparc_data(data_directory)
        
        n_galaxies = len(galaxies)
        total_points = sum(g['n_points'] for g in galaxies)
        point_counts = [g['n_points'] for g in galaxies]
        max_radii = [g['r_max'] for g in galaxies]
        
        integrity_summary = {
            'status': 'valid',
            'n_galaxies': n_galaxies,
            'total_data_points': total_points,
            'avg_points_per_galaxy': np.mean(point_counts),
            'min_points': np.min(point_counts),
            'max_points': np.max(point_counts),
            'avg_max_radius_kpc': np.mean(max_radii),
            'data_source': 'Real SPARC observations'
        }
        
        return integrity_summary
        
    except Exception as e:
        return {
            'status': 'error',
            'error_message': str(e),
            'data_source': 'Unknown'
        }