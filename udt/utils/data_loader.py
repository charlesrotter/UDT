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
    list
        List of galaxy data dictionaries
    """
    if not os.path.exists(data_directory):
        raise ValueError(f"Data directory not found: {data_directory}")
    
    # Look for the MassModels file
    mass_models_file = None
    for filename in ['MassModels_Lelli2016c.mrt', 'MassModels_Lelli2016c.txt']:
        filepath = os.path.join(data_directory, filename)
        if os.path.exists(filepath):
            mass_models_file = filepath
            break
    
    if not mass_models_file:
        print("MassModels file not found. Looking for individual galaxy files...")
        # Fall back to individual files
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
    
    # Parse the MassModels file
    print(f"Loading from MassModels file: {mass_models_file}")
    
    try:
        # Read the file, skipping header lines
        with open(mass_models_file, 'r') as f:
            lines = f.readlines()
        
        # Find data start (after the header)
        data_start = 0
        for i, line in enumerate(lines):
            if line.strip().startswith('---'):
                data_start = i + 1
                break
        
        # Parse data lines
        galaxies_dict = {}
        for line in lines[data_start:]:
            if line.strip() and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 9:
                    galaxy_name = parts[0]
                    try:
                        r = float(parts[2])  # Galactocentric radius
                        v_obs = float(parts[3])  # Observed velocity
                        v_err = float(parts[4])  # Velocity error
                        
                        # Basic quality checks
                        if r > 0 and v_obs > 0 and v_err > 0:
                            if galaxy_name not in galaxies_dict:
                                galaxies_dict[galaxy_name] = {
                                    'name': galaxy_name,
                                    'radius': [],
                                    'velocity': [],
                                    'velocity_error': []
                                }
                            
                            galaxies_dict[galaxy_name]['radius'].append(r)
                            galaxies_dict[galaxy_name]['velocity'].append(v_obs)
                            galaxies_dict[galaxy_name]['velocity_error'].append(v_err)
                    except (ValueError, IndexError):
                        continue
        
        # Convert to proper format
        galaxies = []
        for galaxy_name, data in galaxies_dict.items():
            if len(data['radius']) > 3:  # Minimum points for fitting
                galaxies.append({
                    'name': galaxy_name,
                    'radius': np.array(data['radius']),
                    'velocity': np.array(data['velocity']),
                    'velocity_error': np.array(data['velocity_error']),
                    'n_points': len(data['radius'])
                })
        
        print(f"Successfully loaded {len(galaxies)} galaxies from MassModels file")
        return galaxies
        
    except Exception as e:
        print(f"Error parsing MassModels file: {e}")
        return []


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
        
        # Extract SN name and redshift from first line
        # Format: SN2004dt 0.019700 30.553208 -0.097639 (name, redshift, RA, Dec)
        if len(lines) == 0:
            return None
            
        header = lines[0].strip().split()
        if len(header) < 2:
            return None
            
        sn_name = header[0]
        try:
            redshift = float(header[1])
        except ValueError:
            return None
        
        # Parse B-band photometry data
        b_band_data = []
        in_b_filter = False
        
        for line in lines[1:]:  # Skip header line
            line = line.strip()
            if line.startswith('filter B') and not line.startswith('filter BV'):
                in_b_filter = True
                continue
            elif line.startswith('filter') and in_b_filter:
                break
            elif in_b_filter and line and not line.startswith('filter'):
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        mjd = float(parts[0])
                        mag = float(parts[1])
                        err = float(parts[2])
                        
                        # Quality check: reasonable magnitude range for Type Ia
                        if 10.0 < mag < 25.0 and 0.001 < err < 1.0:
                            b_band_data.append((mjd, mag, err))
                    except ValueError:
                        continue
        
        if redshift and len(b_band_data) >= 3:  # Need minimum observations
            # Convert to numpy array for easier manipulation
            b_band_data = np.array(b_band_data)
            
            # Find peak magnitude (minimum magnitude = maximum brightness)
            peak_idx = np.argmin(b_band_data[:, 1])
            
            return {
                'name': sn_name,
                'redshift': redshift,
                'B_peak_raw': b_band_data[peak_idx, 1],  # Peak magnitude
                'B_error_raw': b_band_data[peak_idx, 2],  # Peak error
                'n_observations': len(b_band_data),
                'time_span': b_band_data[:, 0].max() - b_band_data[:, 0].min(),
                'mag_range': b_band_data[:, 1].max() - b_band_data[:, 1].min()
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
        # Read the file with proper header handling
        df = pd.read_csv(pantheon_file, delim_whitespace=True, comment='#')
        
        # Use only the columns we need (including raw mB for contamination prevention)
        required_columns = ['zCMB', 'm_b_corr', 'm_b_corr_err_DIAG', 'mB', 'mBERR']
        
        # Check if required columns exist
        if not all(col in df.columns for col in required_columns):
            print(f"Available columns: {list(df.columns)}")
            print("Required columns not found. Using available columns...")
            # Try alternative column names
            if 'zcmb' in df.columns:
                df = df.rename(columns={'zcmb': 'zCMB'})
            if 'mb' in df.columns:
                df = df.rename(columns={'mb': 'm_b_corr'})
            if 'dmb' in df.columns:
                df = df.rename(columns={'dmb': 'm_b_corr_err_DIAG'})
        
        # Filter for available columns from required list
        available_columns = [col for col in required_columns if col in df.columns]
        df = df[available_columns].copy()
        
        # Convert to numeric and handle errors
        for col in available_columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # Remove rows with NaN values
        df = df.dropna()
        
        # Filter for quality
        df = df[(df['zCMB'] > 0.001) & (df['zCMB'] < 2.5)]
        df = df[df['m_b_corr_err_DIAG'] < 0.5]  # Reasonable error cut
        
        print(f"Loaded {len(df)} supernovae from Pantheon+")
        print(f"Redshift range: {df['zCMB'].min():.4f} - {df['zCMB'].max():.4f}")
        
        return df
        
    except Exception as e:
        print(f"Error loading Pantheon+ data: {e}")
        import traceback
        traceback.print_exc()
        return pd.DataFrame()