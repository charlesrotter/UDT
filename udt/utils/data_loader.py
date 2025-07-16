"""
Data loading utilities for UDT analyses.

Provides functions to load and parse SPARC galaxy data and
supernova observations from various formats.
"""

import os
import glob
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')


def load_sparc_galaxy(galaxy_file):
    """
    Load individual SPARC galaxy rotation curve data.
    
    Parameters
    ----------
    galaxy_file : str
        Path to galaxy data file
        
    Returns
    -------
    dict
        Galaxy data including name, radius, velocity, errors
    """
    try:
        # Read the file
        with open(galaxy_file, 'r') as f:
            lines = f.readlines()
        
        # Extract galaxy name from filename
        galaxy_name = os.path.basename(galaxy_file).replace('.mrt', '').replace('.txt', '')
        
        # Parse data lines
        data_lines = [line for line in lines if not line.startswith('#') and line.strip()]
        
        radius = []
        velocity = []
        velocity_error = []
        
        for line in data_lines:
            parts = line.strip().split()
            if len(parts) >= 3:
                try:
                    r = float(parts[0])
                    v = float(parts[1])
                    v_err = float(parts[2])
                    
                    # Basic data quality checks
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
                'n_points': len(radius)
            }
        else:
            return None
            
    except Exception as e:
        print(f"Error loading {galaxy_file}: {e}")
        return None


def load_sparc_database(data_directory):
    """
    Load all SPARC galaxy data from a directory.
    
    Parameters
    ----------
    data_directory : str
        Path to SPARC data directory
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing all galaxy data
    """
    if not os.path.exists(data_directory):
        raise ValueError(f"Data directory not found: {data_directory}")
    
    # Find all galaxy files
    galaxy_files = glob.glob(os.path.join(data_directory, "*.mrt"))
    if not galaxy_files:
        galaxy_files = glob.glob(os.path.join(data_directory, "*rotmod*.txt"))
    
    print(f"Found {len(galaxy_files)} galaxy files")
    
    # Load each galaxy
    galaxies = []
    for file_path in galaxy_files:
        galaxy_data = load_sparc_galaxy(file_path)
        if galaxy_data is not None:
            galaxies.append(galaxy_data)
    
    print(f"Successfully loaded {len(galaxies)} galaxies")
    
    return galaxies


def load_csp_supernova(sn_file):
    """
    Load individual CSP supernova data.
    
    Parameters
    ----------
    sn_file : str
        Path to supernova data file
        
    Returns
    -------
    dict
        Supernova data including redshift and photometry
    """
    try:
        # Read file
        with open(sn_file, 'r') as f:
            lines = f.readlines()
        
        # Extract SN name
        sn_name = os.path.basename(sn_file).split('_')[0]
        
        # Look for redshift
        redshift = None
        for line in lines[:50]:  # Check header area
            if 'z' in line.lower() or 'redshift' in line.lower():
                parts = line.split()
                for i, part in enumerate(parts):
                    if part.lower() in ['z', 'z:', 'redshift:', 'redshift']:
                        if i + 1 < len(parts):
                            try:
                                redshift = float(parts[i + 1])
                                break
                            except ValueError:
                                continue
        
        # Parse photometry data
        data_start = False
        time_data = []
        b_mag = []
        b_err = []
        
        for line in lines:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split()
                if len(parts) >= 6:  # Typical format has multiple columns
                    try:
                        mjd = float(parts[0])
                        # Look for B-band data (usually in specific columns)
                        if len(parts) > 3:
                            b = float(parts[2])  # Adjust based on actual format
                            b_e = float(parts[3])
                            if 10 < b < 30 and 0 < b_e < 2:  # Reasonable magnitude range
                                time_data.append(mjd)
                                b_mag.append(b)
                                b_err.append(b_e)
                    except ValueError:
                        continue
        
        if redshift and len(b_mag) > 0:
            # Find peak magnitude
            peak_idx = np.argmin(b_mag)
            return {
                'name': sn_name,
                'redshift': redshift,
                'B_peak_raw': b_mag[peak_idx],
                'B_error_raw': b_err[peak_idx],
                'n_observations': len(b_mag),
                'time_span': max(time_data) - min(time_data) if len(time_data) > 1 else 0
            }
        else:
            return None
            
    except Exception as e:
        print(f"Error loading {sn_file}: {e}")
        return None


def load_csp_database(data_directory):
    """
    Load all CSP supernova data from a directory.
    
    Parameters
    ----------
    data_directory : str
        Path to CSP data directory
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing all supernova data
    """
    if not os.path.exists(data_directory):
        raise ValueError(f"Data directory not found: {data_directory}")
    
    # Find all SN files
    sn_files = glob.glob(os.path.join(data_directory, "SN*_snpy.txt"))
    print(f"Found {len(sn_files)} supernova files")
    
    # Load each supernova
    supernovae = []
    for file_path in sn_files:
        sn_data = load_csp_supernova(file_path)
        if sn_data is not None:
            supernovae.append(sn_data)
    
    print(f"Successfully loaded {len(supernovae)} supernovae")
    
    # Convert to DataFrame
    if supernovae:
        df = pd.DataFrame(supernovae)
        return df
    else:
        return pd.DataFrame()


def load_pantheon_data(pantheon_file):
    """
    Load Pantheon+ supernova catalog.
    
    Parameters
    ----------
    pantheon_file : str
        Path to Pantheon+ data file
        
    Returns
    -------
    pd.DataFrame
        DataFrame with supernova data
    """
    try:
        # Read the file, skipping comment lines
        df = pd.read_csv(pantheon_file, delim_whitespace=True, comment='#',
                        names=['name', 'zcmb', 'zhel', 'mb', 'dmb', 'x1', 'dx1', 
                              'c', 'dc', 'mass', 'ra', 'dec', 'host', 'survey'])
        
        # Filter for quality
        df = df[(df['zcmb'] > 0.001) & (df['zcmb'] < 2.5)]
        df = df[df['dmb'] < 0.5]  # Reasonable error cut
        
        print(f"Loaded {len(df)} supernovae from Pantheon+")
        print(f"Redshift range: {df['zcmb'].min():.4f} - {df['zcmb'].max():.4f}")
        
        return df
        
    except Exception as e:
        print(f"Error loading Pantheon+ data: {e}")
        return pd.DataFrame()