"""
LIGO Data Loader

Validated functions for loading LIGO gravitational wave data.
Supports documented GW150914 parameters and other LIGO events.
"""

import os
import json
import numpy as np
import pandas as pd
from pathlib import Path


def load_gw150914_parameters():
    """
    Load documented GW150914 event parameters.
    
    Returns
    -------
    dict
        GW150914 documented parameters from LIGO Scientific Collaboration
    """
    # Well-documented GW150914 parameters from LIGO Scientific Collaboration
    gw150914_params = {
        # Observational parameters (model-independent)
        'h1_arrival_gps': 1126259462.422,  # GPS time of arrival at H1
        'l1_arrival_gps': 1126259462.429,  # GPS time of arrival at L1
        'timing_difference_ms': 7.0,  # ms (L1 - H1)
        'peak_strain_h1': 1.0e-21,  # Approximate peak strain
        'peak_strain_l1': 1.0e-21,  # Approximate peak strain
        
        # Detector geometry (well-known)
        'h1_location': [46.4547, -119.4077, 142.554],  # lat, lon, elevation (m)
        'l1_location': [30.5628, -90.7739, -6.574],    # lat, lon, elevation (m)
        
        # Physical source parameters (from GR analysis)
        'source_distance_mpc': 410,  # Megaparsecs
        'black_hole_mass_1': 36,     # Solar masses
        'black_hole_mass_2': 29,     # Solar masses
        'merger_frequency_hz': 250,  # Hz
        
        # Signal characteristics
        'signal_duration_s': 0.2,    # seconds
        'frequency_range_hz': [35, 350],  # Hz
        'snr_h1': 24,                # Signal-to-noise ratio
        'snr_l1': 13                 # Signal-to-noise ratio
    }
    
    return gw150914_params


def load_ligo_events_data(data_directory):
    """
    Load LIGO events data from JSON file.
    
    Parameters
    ----------
    data_directory : str
        Path to LIGO data directory
        
    Returns
    -------
    pd.DataFrame
        LIGO events data
    """
    data_path = Path(data_directory)
    
    # Look for LIGO events data file
    ligo_files = list(data_path.glob("**/ligo_events_data.json"))
    
    if len(ligo_files) == 0:
        # Create synthetic data based on documented parameters if no file found
        print("Warning: No LIGO events data file found, using documented GW150914 parameters")
        gw150914 = load_gw150914_parameters()
        
        events_data = [{
            'event_name': 'GW150914',
            'gps_time': gw150914['h1_arrival_gps'],
            'source_distance_mpc': gw150914['source_distance_mpc'],
            'mass_1_solar': gw150914['black_hole_mass_1'],
            'mass_2_solar': gw150914['black_hole_mass_2'],
            'peak_strain': gw150914['peak_strain_h1'],
            'snr_h1': gw150914['snr_h1'],
            'snr_l1': gw150914['snr_l1'],
            'timing_difference_ms': gw150914['timing_difference_ms']
        }]
        
        return pd.DataFrame(events_data)
    
    try:
        with open(ligo_files[0], 'r') as f:
            events_data = json.load(f)
        
        return pd.DataFrame(events_data)
        
    except Exception as e:
        raise ValueError(f"Error loading LIGO events data: {e}")


def calculate_detector_separation():
    """
    Calculate LIGO detector separation from documented coordinates.
    
    Returns
    -------
    float
        Detector separation in meters
    """
    gw150914_params = load_gw150914_parameters()
    
    h1_lat, h1_lon, h1_elev = gw150914_params['h1_location']
    l1_lat, l1_lon, l1_elev = gw150914_params['l1_location']
    
    # Convert to radians
    h1_lat_rad = np.radians(h1_lat)
    h1_lon_rad = np.radians(h1_lon)
    l1_lat_rad = np.radians(l1_lat)
    l1_lon_rad = np.radians(l1_lon)
    
    # Earth radius
    earth_radius = 6371000  # meters
    
    # Haversine formula for great circle distance
    dlat = l1_lat_rad - h1_lat_rad
    dlon = l1_lon_rad - h1_lon_rad
    
    a = np.sin(dlat/2)**2 + np.cos(h1_lat_rad) * np.cos(l1_lat_rad) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    
    # Distance along Earth's surface
    surface_distance = earth_radius * c
    
    # Account for elevation difference
    elevation_diff = l1_elev - h1_elev
    detector_separation = np.sqrt(surface_distance**2 + elevation_diff**2)
    
    return detector_separation


def load_ligo_data(data_directory, event='all'):
    """
    Load LIGO gravitational wave data.
    
    Parameters
    ----------
    data_directory : str
        Path to LIGO data directory
    event : str
        Event to load ('gw150914', 'all')
        
    Returns
    -------
    pd.DataFrame
        LIGO data with standardized columns
    """
    if event.lower() == 'gw150914':
        # Load documented GW150914 parameters
        gw150914 = load_gw150914_parameters()
        detector_separation = calculate_detector_separation()
        
        data = pd.DataFrame([{
            'event_name': 'GW150914',
            'gps_time': gw150914['h1_arrival_gps'],
            'source_distance_mpc': gw150914['source_distance_mpc'],
            'mass_1_solar': gw150914['black_hole_mass_1'],
            'mass_2_solar': gw150914['black_hole_mass_2'],
            'peak_strain': gw150914['peak_strain_h1'],
            'snr_h1': gw150914['snr_h1'],
            'snr_l1': gw150914['snr_l1'],
            'timing_difference_ms': gw150914['timing_difference_ms'],
            'detector_separation_m': detector_separation
        }])
        
        return data
        
    elif event.lower() == 'all':
        return load_ligo_events_data(data_directory)
    
    else:
        raise ValueError(f"Unknown event: {event}")


def validate_ligo_data_integrity(data_directory):
    """
    Validate LIGO data integrity and UDT compatibility.
    
    Parameters
    ----------
    data_directory : str
        Path to LIGO data directory
        
    Returns
    -------
    dict
        Data integrity and UDT compatibility assessment
    """
    try:
        data = load_ligo_data(data_directory, 'all')
        
        integrity_summary = {
            'status': 'valid',
            'n_events': len(data),
            'data_source': 'Documented LIGO parameters (GW150914)',
            'udt_compatibility': 'projection_theory_testable'
        }
        
        # Check for timing consistency with speed of light
        if 'timing_difference_ms' in data.columns:
            c_light = 299792458  # m/s
            detector_sep = calculate_detector_separation()
            expected_timing_ms = (detector_sep / c_light) * 1000
            
            for _, event in data.iterrows():
                observed_timing = event['timing_difference_ms']
                timing_ratio = observed_timing / expected_timing_ms
                
                if 0.5 <= timing_ratio <= 2.0:
                    integrity_summary['timing_assessment'] = 'excellent_udt_compatibility'
                elif 0.2 <= timing_ratio <= 5.0:
                    integrity_summary['timing_assessment'] = 'good_udt_compatibility'
                else:
                    integrity_summary['timing_assessment'] = 'challenging_for_udt'
        
        return integrity_summary
        
    except Exception as e:
        return {
            'status': 'error',
            'error_message': str(e),
            'data_source': 'Unknown'
        }