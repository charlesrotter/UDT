#!/usr/bin/env python3
"""
Direct CMS Dataset Download
==========================

Download specific CMS datasets using known record IDs from CERN Open Data Portal.
Focus on educational and smaller datasets suitable for UDT analysis.

Author: Charles Rotter
Date: 2025-07-18
"""

import os
import requests
import json
import time
from urllib.parse import urljoin

class DirectCMSDownloader:
    def __init__(self):
        self.base_url = "https://opendata.cern.ch/"
        self.download_dir = "C:/UDT/data/cms_open_data/"
        
        os.makedirs(self.download_dir, exist_ok=True)
        
        print("DIRECT CMS DATASET DOWNLOADER")
        print("=" * 29)
        print(f"Target directory: {self.download_dir}")
        print()
    
    def download_record(self, record_id):
        """Download a specific record by ID."""
        print(f"DOWNLOADING RECORD {record_id}")
        print("-" * 30)
        
        try:
            # Get record metadata
            url = f"https://opendata.cern.ch/api/records/{record_id}"
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            record = response.json()
            metadata = record.get('metadata', {})
            
            title = metadata.get('title', f'Record_{record_id}')
            print(f"Title: {title}")
            
            # Check for files
            files = metadata.get('files', [])
            print(f"Files available: {len(files)}")
            
            if not files:
                print("No files found in this record")
                return False
            
            # Create directory for this record
            safe_title = "".join(c for c in title if c.isalnum() or c in (' ', '-', '_')).rstrip()[:50]
            record_dir = os.path.join(self.download_dir, f"record_{record_id}_{safe_title}")
            os.makedirs(record_dir, exist_ok=True)
            
            # Save metadata
            metadata_file = os.path.join(record_dir, 'metadata.json')
            with open(metadata_file, 'w') as f:
                json.dump(record, f, indent=2)
            
            # Download files (limit to first few)
            downloaded = 0
            total_size = 0
            max_files = 3
            max_size = 500e6  # 500 MB limit
            
            for i, file_info in enumerate(files):
                if downloaded >= max_files or total_size >= max_size:
                    break
                
                file_size = file_info.get('size', 0)
                if total_size + file_size > max_size:
                    continue
                
                if self.download_file(file_info, record_dir):
                    downloaded += 1
                    total_size += file_size
            
            print(f"Downloaded {downloaded} files to {record_dir}")
            return True
            
        except Exception as e:
            print(f"Error downloading record {record_id}: {e}")
            return False
    
    def download_file(self, file_info, record_dir):
        """Download a single file."""
        uri = file_info.get('uri', '')
        filename = file_info.get('filename', 'unknown_file')
        file_size = file_info.get('size', 0)
        
        if not uri:
            print(f"  No URI for {filename}")
            return False
        
        # Construct download URL
        if uri.startswith('root://'):
            print(f"  Skipping ROOT protocol file: {filename}")
            return False
        
        if uri.startswith('http'):
            download_url = uri
        else:
            download_url = urljoin(self.base_url, uri.lstrip('/'))
        
        local_path = os.path.join(record_dir, filename)
        
        print(f"  Downloading: {filename} ({file_size/1e6:.1f} MB)")
        
        try:
            # Check if already exists
            if os.path.exists(local_path) and os.path.getsize(local_path) == file_size:
                print(f"    Already exists, skipping")
                return True
            
            # Download
            response = requests.get(download_url, stream=True, timeout=120)
            response.raise_for_status()
            
            with open(local_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
            
            print(f"    Downloaded successfully")
            return True
            
        except Exception as e:
            print(f"    Error: {e}")
            if os.path.exists(local_path):
                os.remove(local_path)
            return False
    
    def download_known_good_records(self):
        """Download known good CMS records suitable for analysis."""
        print("DOWNLOADING KNOWN GOOD CMS RECORDS")
        print("=" * 35)
        
        # Known record IDs for educational/analysis datasets
        # These are from CMS open data documentation
        target_records = [
            5500,   # CMS 2010 collision data (small sample)
            5502,   # CMS 2011 collision data (small sample)
            12341,  # CMS MC data (simulated)
            12342,  # CMS MC data (simulated)
            14220,  # CMS educational datasets
            17111,  # CMS analysis examples
            30500,  # CMS open data examples
        ]
        
        successful_downloads = 0
        
        for record_id in target_records:
            print()
            success = self.download_record(record_id)
            if success:
                successful_downloads += 1
            time.sleep(2)  # Be respectful to server
        
        print(f"\n\nSUMMARY:")
        print(f"Successfully downloaded {successful_downloads}/{len(target_records)} records")
        
        return successful_downloads > 0
    
    def download_cms_education_data(self):
        """Download CMS educational datasets specifically."""
        print("DOWNLOADING CMS EDUCATIONAL DATASETS")
        print("=" * 33)
        
        # Try some educational record ranges
        education_records = [
            545,    # CMS educational examples
            5200,   # CMS 2010 data subset
            5503,   # CMS 2011 data subset
            5004,   # CMS 2012 data subset
            14210,  # Educational materials
            12102,  # Simulated data for education
        ]
        
        successful = 0
        
        for record_id in education_records:
            try:
                success = self.download_record(record_id)
                if success:
                    successful += 1
                time.sleep(1)
            except Exception as e:
                print(f"Record {record_id} failed: {e}")
                continue
        
        return successful > 0
    
    def run_download(self):
        """Run the complete download process."""
        print("STARTING DIRECT CMS DATA DOWNLOAD")
        print("=" * 34)
        print("Strategy: Download known record IDs for smaller datasets")
        print()
        
        # Try educational datasets first
        success1 = self.download_cms_education_data()
        
        print("\n" + "="*60)
        
        # Try known good records
        success2 = self.download_known_good_records()
        
        print("\n" + "="*60)
        print("DIRECT CMS DOWNLOAD COMPLETE")
        print("="*60)
        
        if success1 or success2:
            print(f"Downloaded datasets to: {self.download_dir}")
            print("\nThese files can be used for UDT validation:")
            print("- Test particle momentum distributions under UDT geometry")
            print("- Compare kinematic correlations with Standard Model predictions")
            print("- Analyze collision data for deviations from SM expectations")
        else:
            print("No datasets successfully downloaded.")
            print("May need to check CERN Open Data Portal for current record IDs.")

def main():
    """Main download routine."""
    downloader = DirectCMSDownloader()
    downloader.run_download()

if __name__ == "__main__":
    main()