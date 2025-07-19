#!/usr/bin/env python3
"""
Download Accessible CMS Files
=============================

Download specific CMS files that are web-accessible for UDT analysis.
Focus on HTTP downloadable files rather than ROOT protocol.

Author: Charles Rotter
Date: 2025-07-18
"""

import os
import requests
import json
import time
from urllib.parse import urljoin

class AccessibleCMSDownloader:
    def __init__(self):
        self.download_dir = "C:/UDT/data/cms_accessible/"
        os.makedirs(self.download_dir, exist_ok=True)
        
        print("ACCESSIBLE CMS FILE DOWNLOADER")
        print("=" * 30)
        print(f"Download directory: {self.download_dir}")
        print()
    
    def download_http_file(self, url, filename, max_size_gb=1):
        """Download a file via HTTP with size limit."""
        local_path = os.path.join(self.download_dir, filename)
        
        print(f"DOWNLOADING: {filename}")
        print(f"URL: {url}")
        
        try:
            # Get file info first
            response = requests.head(url, timeout=30)
            response.raise_for_status()
            
            file_size = int(response.headers.get('content-length', 0))
            print(f"Size: {file_size/1e9:.2f} GB")
            
            if file_size > max_size_gb * 1e9:
                print(f"File too large (>{max_size_gb} GB), skipping")
                return False
            
            # Check if already exists
            if os.path.exists(local_path) and os.path.getsize(local_path) == file_size:
                print("File already exists with correct size")
                return True
            
            # Download with progress
            print("Downloading...")
            response = requests.get(url, stream=True, timeout=300)
            response.raise_for_status()
            
            downloaded = 0
            with open(local_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        
                        if file_size > 0:
                            progress = (downloaded / file_size) * 100
                            print(f"Progress: {progress:.1f}%", end='\r')
            
            print(f"\nDownloaded successfully: {downloaded/1e6:.1f} MB")
            return True
            
        except Exception as e:
            print(f"Error downloading {filename}: {e}")
            if os.path.exists(local_path):
                os.remove(local_path)
            return False
    
    def try_cern_opendata_direct_links(self):
        """Try direct links to CERN open data files."""
        print("TRYING CERN OPEN DATA DIRECT LINKS")
        print("=" * 35)
        
        # Known smaller files that might be HTTP accessible
        target_files = [
            {
                "name": "CMS_validated_runs_2015.txt",
                "url": "https://opendata.cern.ch/record/14210/files/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt",
                "description": "CMS validated runs list for 2015"
            },
            {
                "name": "CMS_validated_runs_2016.txt", 
                "url": "https://opendata.cern.ch/record/14220/files/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
                "description": "CMS validated runs list for 2016"
            },
            {
                "name": "cms_higgs_analysis_example.py",
                "url": "https://opendata.cern.ch/record/5500/files/HiggsExample20112012/HiggsDemoAnalyzer.cc",
                "description": "CMS Higgs analysis example code"
            }
        ]
        
        successful = 0
        
        for file_info in target_files:
            print(f"\nTrying: {file_info['description']}")
            success = self.download_http_file(
                file_info['url'], 
                file_info['name'],
                max_size_gb=0.1  # Small files only
            )
            if success:
                successful += 1
            time.sleep(1)
        
        return successful
    
    def download_cms_educational_resources(self):
        """Download CMS educational and documentation resources."""
        print("DOWNLOADING CMS EDUCATIONAL RESOURCES")
        print("=" * 36)
        
        # Educational resources that might be directly accessible
        educational_urls = [
            {
                "name": "cms_open_data_guide.html",
                "url": "https://opendata.cern.ch/docs/cms-guide-for-research",
                "description": "CMS Open Data Research Guide"
            },
            {
                "name": "cms_physics_objects.html", 
                "url": "https://opendata.cern.ch/docs/cms-physics-objects-2016",
                "description": "CMS Physics Objects Documentation"
            }
        ]
        
        successful = 0
        
        for resource in educational_urls:
            print(f"\nDownloading: {resource['description']}")
            try:
                response = requests.get(resource['url'], timeout=30)
                response.raise_for_status()
                
                local_path = os.path.join(self.download_dir, resource['name'])
                with open(local_path, 'w', encoding='utf-8') as f:
                    f.write(response.text)
                
                print(f"Downloaded: {len(response.text)/1000:.1f} KB")
                successful += 1
                
            except Exception as e:
                print(f"Error: {e}")
        
        return successful
    
    def explore_alternative_cms_data_sources(self):
        """Explore alternative sources for CMS data."""
        print("EXPLORING ALTERNATIVE CMS DATA SOURCES")
        print("=" * 38)
        
        # Try HEPData for CMS results
        hepdata_urls = [
            {
                "name": "cms_higgs_data.json",
                "url": "https://www.hepdata.net/record/ins1696414",
                "description": "CMS Higgs measurement data from HEPData"
            }
        ]
        
        print("Checking HEPData for CMS results...")
        
        for data_source in hepdata_urls:
            try:
                response = requests.get(data_source['url'], timeout=30)
                if response.status_code == 200:
                    print(f"Found: {data_source['description']}")
                    # This would require parsing HEPData JSON structure
                else:
                    print(f"Not accessible: {data_source['description']}")
            except:
                print(f"Error accessing: {data_source['description']}")
        
        return 0
    
    def create_cms_data_summary(self):
        """Create a summary of what CMS data we have access to."""
        print("CREATING CMS DATA ACCESS SUMMARY")
        print("=" * 31)
        
        summary = {
            "download_directory": self.download_dir,
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "accessible_data_types": [
                "Metadata from CERN Open Data Portal records",
                "Run validation lists (JSON format)",
                "Analysis code examples", 
                "Documentation and guides"
            ],
            "data_access_challenges": [
                "Most collision data files use ROOT protocol (not HTTP)",
                "Large file sizes (GB range) for collision datasets",
                "Requires special ROOT/CERN software for full access"
            ],
            "recommendations_for_udt_analysis": [
                "Use metadata to understand data structure",
                "Focus on kinematic variables described in dataset semantics",
                "Consider requesting smaller samples from CERN collaboration",
                "Use LIGO data (already available) as primary validation"
            ],
            "files_downloaded": []
        }
        
        # List downloaded files
        if os.path.exists(self.download_dir):
            for item in os.listdir(self.download_dir):
                if os.path.isfile(os.path.join(self.download_dir, item)):
                    file_path = os.path.join(self.download_dir, item)
                    file_size = os.path.getsize(file_path)
                    summary["files_downloaded"].append({
                        "filename": item,
                        "size_bytes": file_size,
                        "size_mb": round(file_size/1e6, 2)
                    })
        
        # Save summary
        summary_file = os.path.join(self.download_dir, "cms_data_access_summary.json")
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        print(f"Summary saved to: {summary_file}")
        return summary
    
    def run_download(self):
        """Run the complete accessible CMS download process."""
        print("STARTING ACCESSIBLE CMS DATA DOWNLOAD")
        print("=" * 37)
        print("Strategy: Download HTTP-accessible CMS files and resources")
        print()
        
        total_successful = 0
        
        # Try direct links
        total_successful += self.try_cern_opendata_direct_links()
        
        print("\n" + "="*50)
        
        # Educational resources
        total_successful += self.download_cms_educational_resources()
        
        print("\n" + "="*50)
        
        # Alternative sources
        self.explore_alternative_cms_data_sources()
        
        print("\n" + "="*50)
        
        # Create summary
        summary = self.create_cms_data_summary()
        
        print("\n" + "="*60)
        print("ACCESSIBLE CMS DOWNLOAD COMPLETE")
        print("="*60)
        
        if total_successful > 0:
            print(f"Successfully downloaded {total_successful} files/resources")
            print(f"Files saved to: {self.download_dir}")
            print(f"Downloaded files: {len(summary['files_downloaded'])}")
        else:
            print("No files successfully downloaded via HTTP")
            print("CMS collision data requires ROOT protocol access")
        
        print("\nRecommendations for UDT validation:")
        print("1. Use existing LIGO data (excellent for geometric effects)")
        print("2. Focus on CMS metadata for understanding data structures")
        print("3. Consider contacting CMS collaboration for sample data")
        print("4. Use CMB data (already validated UDT 15.67x better than LCDM)")

def main():
    """Main accessible download routine."""
    downloader = AccessibleCMSDownloader()
    downloader.run_download()

if __name__ == "__main__":
    main()