"""
Download and process SPARC galaxy database
Real data for Information Curvature Theory testing
"""

import pandas as pd
import numpy as np
import os
import urllib.request

def download_sparc_database():
    """
    Download SPARC galaxy rotation curve database
    """
    print("📊 DOWNLOADING REAL SPARC DATABASE")
    print("=" * 40)
    
    # SPARC database URLs (check these are current)
    urls = {
        'table': 'http://astroweb.cwru.edu/SPARC/MasterTable.txt',
        'readme': 'http://astroweb.cwru.edu/SPARC/ReadMe.txt'
    }
    
    # Create data directory
    os.makedirs('data/sparc_raw', exist_ok=True)
    
    # Download files
    for name, url in urls.items():
        local_file = f'data/sparc_raw/sparc_{name}.txt'
        try:
            print(f"Downloading {name}...")
            urllib.request.urlretrieve(url, local_file)
            print(f"✅ Saved to {local_file}")
        except Exception as e:
            print(f"❌ Failed to download {name}: {e}")
            print("🔍 Check URL or download manually from:")
            print("   http://astroweb.cwru.edu/SPARC/")
    
    return True

def load_sparc_table():
    """
    Load and parse SPARC master table
    """
    print("\n📋 LOADING SPARC MASTER TABLE")
    print("-" * 30)
    
    file_path = 'data/sparc_raw/sparc_table.txt'
    
    if not os.path.exists(file_path):
        print("❌ SPARC table not found. Run download_sparc_database() first")
        return None
    
    try:
        # Read SPARC table (format may need adjustment)
        # Note: Real implementation would need proper column parsing
        df = pd.read_csv(file_path, sep='\s+', comment='#', 
                        names=['Galaxy', 'D', 'e_D', 'f_D', 'T', 'i', 'e_i', 
                              'L36', 'e_L36', 'Reff', 'SBeff', 'Rdisk', 'SBdisk',
                              'MHI', 'RHI', 'Vflat', 'e_Vflat', 'Q', 'Ref'])
        
        print(f"✅ Loaded {len(df)} galaxies from SPARC")
        print(f"Columns: {list(df.columns)}")
        print(f"Distance range: {df['D'].min():.1f} - {df['D'].max():.1f} Mpc")
        
        return df
        
    except Exception as e:
        print(f"❌ Error loading SPARC table: {e}")
        print("🔍 Check file format or download fresh copy")
        return None

def download_individual_rotation_curves(galaxy_list):
    """
    Download individual rotation curve files for selected galaxies
    """
    print(f"\n🌀 DOWNLOADING ROTATION CURVES")
    print("-" * 35)
    
    base_url = "http://astroweb.cwru.edu/SPARC/RotationCurves/"
    os.makedirs('data/sparc_curves', exist_ok=True)
    
    successful_downloads = []
    
    for galaxy in galaxy_list:
        try:
            # SPARC naming convention (may need adjustment)
            filename = f"{galaxy}_rotcur.dat"
            url = base_url + filename
            local_path = f"data/sparc_curves/{filename}"
            
            print(f"Downloading {galaxy}...")
            urllib.request.urlretrieve(url, local_path)
            successful_downloads.append(galaxy)
            print(f"✅ {galaxy}")
            
        except Exception as e:
            print(f"❌ {galaxy}: {e}")
    
    print(f"\n📊 Downloaded {len(successful_downloads)} rotation curves")
    return successful_downloads

def main():
    """
    Main function to set up SPARC data
    """
    print("🚀 SPARC DATABASE SETUP")
    print("Information Curvature Theory - Real Data Testing")
    print("=" * 50)
    
    # Step 1: Download master table
    download_sparc_database()
    
    # Step 2: Load and examine
    df = load_sparc_table()
    
    if df is not None:
        # Step 3: Select promising galaxies for analysis
        print(f"\n🎯 SELECTING GALAXIES FOR ANALYSIS")
        print("-" * 35)
        
        # Quality criteria (adjust as needed)
        good_galaxies = df[
            (df['Q'] >= 2) &  # Quality flag
            (df['D'] > 0) &   # Valid distance
            (df['Vflat'] > 50)  # Reasonable rotation velocity
        ]
        
        print(f"Quality-filtered sample: {len(good_galaxies)} galaxies")
        
        # Select top candidates
        candidates = good_galaxies.head(10)['Galaxy'].tolist()
        print(f"Selected candidates: {candidates}")
        
        # Step 4: Download rotation curves
        successful = download_individual_rotation_curves(candidates)
        
        print(f"\n✅ SETUP COMPLETE")
        print(f"Ready to test Information Curvature Theory on {len(successful)} real galaxies!")
        
        return successful
    
    else:
        print("❌ Setup failed. Check internet connection and SPARC URLs.")
        return []

if __name__ == "__main__":
    galaxies = main()
