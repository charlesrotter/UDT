#!/usr/bin/env python3
"""
CMB Raw Data Inspector
======================

Inspects raw CMB FITS files to understand their structure and content.
This replaces the need for healpy by using astropy for basic FITS operations.

Author: Charles Rotter
Date: 2025-01-17
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import sys
from pathlib import Path

class CMBDataInspector:
    """Inspect CMB FITS files without requiring healpy."""
    
    def __init__(self, fits_file):
        """Initialize with FITS file path."""
        self.fits_file = Path(fits_file)
        self.data = None
        self.header = None
        self.hdulist = None
        
    def load_fits_file(self):
        """Load and inspect FITS file structure."""
        print(f"Loading FITS file: {self.fits_file}")
        print("=" * 60)
        
        try:
            # Open FITS file
            self.hdulist = fits.open(self.fits_file)
            
            print(f"Number of HDUs: {len(self.hdulist)}")
            print()
            
            # Inspect each HDU
            for i, hdu in enumerate(self.hdulist):
                print(f"HDU {i}: {type(hdu).__name__}")
                print(f"  Name: {hdu.name}")
                if hasattr(hdu, 'data') and hdu.data is not None:
                    print(f"  Data shape: {hdu.data.shape}")
                    print(f"  Data type: {hdu.data.dtype}")
                    if hasattr(hdu.data, 'names') and hdu.data.names:
                        print(f"  Column names: {hdu.data.names}")
                if hasattr(hdu, 'header'):
                    print(f"  Header keys: {len(hdu.header)} items")
                print()
                
            return True
            
        except Exception as e:
            print(f"Error loading FITS file: {e}")
            return False
    
    def inspect_primary_data(self):
        """Inspect the primary data extension."""
        if self.hdulist is None:
            print("No FITS file loaded!")
            return
            
        print("INSPECTING PRIMARY DATA")
        print("=" * 40)
        
        # Find data extension (usually extension 1 for HEALPix)
        data_hdu = None
        for i, hdu in enumerate(self.hdulist):
            if hasattr(hdu, 'data') and hdu.data is not None:
                if hasattr(hdu.data, 'names') or len(hdu.data.shape) > 0:
                    data_hdu = hdu
                    print(f"Using HDU {i} for data analysis")
                    break
        
        if data_hdu is None:
            print("No data HDU found!")
            return
            
        self.data = data_hdu.data
        self.header = data_hdu.header
        
        print(f"Data type: {type(self.data)}")
        print(f"Data shape: {self.data.shape}")
        
        if hasattr(self.data, 'names') and self.data.names:
            print(f"Columns: {self.data.names}")
            print()
            
            # Inspect each column
            for name in self.data.names:
                col_data = self.data[name]
                print(f"Column '{name}':")
                print(f"  Shape: {col_data.shape}")
                print(f"  Type: {col_data.dtype}")
                print(f"  Min: {np.min(col_data):.6e}")
                print(f"  Max: {np.max(col_data):.6e}")
                print(f"  Mean: {np.mean(col_data):.6e}")
                print(f"  Std: {np.std(col_data):.6e}")
                print()
        else:
            # Array data
            print(f"Array data shape: {self.data.shape}")
            print(f"Min: {np.min(self.data):.6e}")
            print(f"Max: {np.max(self.data):.6e}")
            print(f"Mean: {np.mean(self.data):.6e}")
            print(f"Std: {np.std(self.data):.6e}")
            print()
    
    def inspect_header(self):
        """Inspect FITS header for important metadata."""
        if self.header is None:
            print("No header loaded!")
            return
            
        print("HEADER INFORMATION")
        print("=" * 40)
        
        # Key HEALPix parameters
        important_keys = [
            'COORDSYS', 'ORDERING', 'NSIDE', 'TUNIT1', 'EXTNAME',
            'CREATOR', 'DATE', 'COMMENT', 'HISTORY'
        ]
        
        for key in important_keys:
            if key in self.header:
                value = self.header[key]
                if isinstance(value, str) and len(value) > 50:
                    print(f"{key}: {value[:50]}...")
                else:
                    print(f"{key}: {value}")
        
        print()
        print("All header keys:")
        for i, key in enumerate(self.header.keys()):
            if i % 8 == 0:
                print()
            print(f"{key:8s}", end=" ")
        print("\n")
    
    def basic_statistics(self):
        """Calculate basic statistics on the data."""
        if self.data is None:
            print("No data loaded!")
            return
            
        print("BASIC STATISTICS")
        print("=" * 40)
        
        if hasattr(self.data, 'names') and self.data.names:
            # Table data - analyze temperature column
            temp_candidates = ['TEMPERATURE', 'I_STOKES', 'TEMP', 'I']
            temp_col = None
            
            for candidate in temp_candidates:
                if candidate in self.data.names:
                    temp_col = candidate
                    break
            
            if temp_col is None and len(self.data.names) > 0:
                temp_col = self.data.names[0]  # Use first column
                
            if temp_col:
                temp_data = self.data[temp_col]
                print(f"Analyzing column: {temp_col}")
                print(f"Number of pixels: {len(temp_data)}")
                print(f"Temperature range: {np.min(temp_data):.6f} to {np.max(temp_data):.6f}")
                print(f"Mean temperature: {np.mean(temp_data):.6f}")
                print(f"RMS: {np.std(temp_data):.6f}")
                print(f"Units: {self.header.get('TUNIT1', 'Unknown')}")
                
                # Check for suspicious values
                n_nan = np.sum(np.isnan(temp_data))
                n_inf = np.sum(np.isinf(temp_data))
                n_zero = np.sum(temp_data == 0)
                
                print(f"Data quality:")
                print(f"  NaN values: {n_nan}")
                print(f"  Infinite values: {n_inf}")
                print(f"  Zero values: {n_zero}")
                print()
                
                return temp_data
        else:
            # Array data
            print(f"Array shape: {self.data.shape}")
            print(f"Data range: {np.min(self.data):.6f} to {np.max(self.data):.6f}")
            print(f"Mean: {np.mean(self.data):.6f}")
            print(f"RMS: {np.std(self.data):.6f}")
            return self.data
    
    def create_basic_plots(self, data, output_dir="results/cmb_raw_analysis"):
        """Create basic plots of the data."""
        os.makedirs(output_dir, exist_ok=True)
        
        print("CREATING BASIC PLOTS")
        print("=" * 40)
        
        # Histogram
        plt.figure(figsize=(12, 8))
        
        plt.subplot(2, 2, 1)
        plt.hist(data, bins=100, alpha=0.7, edgecolor='black')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Number of pixels')
        plt.title('Temperature Distribution')
        plt.grid(True, alpha=0.3)
        
        # Statistics
        plt.subplot(2, 2, 2)
        stats_text = f"""
        Mean: {np.mean(data):.6f} K
        Std: {np.std(data):.6f} K
        Min: {np.min(data):.6f} K
        Max: {np.max(data):.6f} K
        Pixels: {len(data):,}
        """
        plt.text(0.1, 0.5, stats_text, transform=plt.gca().transAxes, 
                fontsize=12, verticalalignment='center',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray"))
        plt.axis('off')
        plt.title('Data Statistics')
        
        # Pixel index vs temperature (sample)
        plt.subplot(2, 2, 3)
        sample_size = min(10000, len(data))
        indices = np.random.choice(len(data), sample_size, replace=False)
        plt.scatter(indices, data[indices], alpha=0.5, s=0.5)
        plt.xlabel('Pixel Index')
        plt.ylabel('Temperature (K)')
        plt.title(f'Temperature vs Pixel Index (sample of {sample_size})')
        plt.grid(True, alpha=0.3)
        
        # Cumulative distribution
        plt.subplot(2, 2, 4)
        sorted_data = np.sort(data)
        cumulative = np.arange(len(sorted_data)) / len(sorted_data)
        plt.plot(sorted_data, cumulative, linewidth=1)
        plt.xlabel('Temperature (K)')
        plt.ylabel('Cumulative Fraction')
        plt.title('Cumulative Distribution')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plot_file = os.path.join(output_dir, 'cmb_raw_data_analysis.png')
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"Plots saved to: {plot_file}")
    
    def run_full_inspection(self):
        """Run complete inspection of CMB data."""
        print("CMB RAW DATA INSPECTION")
        print("=" * 60)
        print(f"File: {self.fits_file}")
        print()
        
        # Load file
        if not self.load_fits_file():
            return False
        
        # Inspect structure
        self.inspect_header()
        self.inspect_primary_data()
        
        # Basic analysis
        data = self.basic_statistics()
        
        if data is not None:
            self.create_basic_plots(data)
        
        # Summary
        print("INSPECTION SUMMARY")
        print("=" * 40)
        print("File successfully loaded and analyzed")
        print(f"Data contains {len(data) if data is not None else 'unknown'} pixels")
        print("Basic plots created")
        print()
        
        return True
    
    def close(self):
        """Close FITS file."""
        if self.hdulist:
            self.hdulist.close()


def main():
    """Main inspection routine."""
    # Default file location
    cmb_file = Path(__file__).parent.parent / "data" / "cmb_raw" / "COM_CMB_IQU-smica-nosz_2048_R3.00_full.fits"
    
    if not cmb_file.exists():
        print(f"CMB file not found: {cmb_file}")
        print("Please download the Planck SMICA file first")
        return False
    
    # Create inspector
    inspector = CMBDataInspector(cmb_file)
    
    try:
        # Run inspection
        success = inspector.run_full_inspection()
        return success
        
    finally:
        inspector.close()


if __name__ == "__main__":
    main()