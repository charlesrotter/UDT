#!/usr/bin/env python3
"""
CMS Open Data Download Script
============================

Downloads CMS open datasets for UDT kinematic analysis.
Focus on smaller, manageable datasets for testing geometric effects.

Author: Charles Rotter  
Date: 2025-07-18
"""

import os
import requests
import json
from urllib.parse import urljoin
import time

class CMSOpenDataDownloader:
    def __init__(self):
        self.base_url = "https://opendata.cern.ch/"
        self.api_url = "https://opendata.cern.ch/api/records/"
        self.download_dir = "C:/UDT/data/cms_open_data/"
        
        # Create download directory
        os.makedirs(self.download_dir, exist_ok=True)
        
        print("CMS OPEN DATA DOWNLOADER")
        print("=" * 24)
        print(f"Download directory: {self.download_dir}")
        print()
    
    def search_datasets(self, query_params):
        """Search for datasets using CERN Open Data API."""
        print(f"SEARCHING FOR DATASETS")
        print("-" * 20)
        
        # API endpoint for searching
        search_url = f"{self.api_url}?type=Dataset&experiment=CMS"
        
        # Add additional parameters
        for key, value in query_params.items():
            search_url += f"&{key}={value}"
        
        print(f"Search URL: {search_url}")
        
        try:
            response = requests.get(search_url, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            hits = data.get('hits', {}).get('hits', [])
            
            print(f"Found {len(hits)} datasets")
            print()
            
            return hits
            
        except Exception as e:
            print(f"Error searching datasets: {e}")
            return []
    
    def get_dataset_info(self, dataset_id):
        """Get detailed information about a specific dataset."""
        try:
            url = f"{self.api_url}{dataset_id}"
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            return response.json()
            
        except Exception as e:
            print(f"Error getting dataset {dataset_id}: {e}")
            return None
    
    def download_small_cms_samples(self):
        """Download small CMS samples suitable for kinematic analysis."""
        print("DOWNLOADING SMALL CMS SAMPLES FOR KINEMATIC ANALYSIS")
        print("=" * 52)
        
        # Look for specific known datasets that are smaller
        target_datasets = [
            # Search for educational or analysis-ready datasets
            {"q": "DoubleMuon Run2016", "size": "10"},
            {"q": "SingleMuon Run2016", "size": "10"}, 
            {"q": "DoubleEG Run2016", "size": "10"},
            {"q": "CMS MiniAOD", "size": "5"},
            {"q": "CMS NanoAOD", "size": "5"}
        ]
        
        all_datasets = []
        
        for search_params in target_datasets:
            datasets = self.search_datasets(search_params)
            all_datasets.extend(datasets)
            time.sleep(1)  # Be respectful to server
        
        # Remove duplicates
        unique_datasets = {}
        for dataset in all_datasets:
            dataset_id = dataset.get('id')
            if dataset_id not in unique_datasets:
                unique_datasets[dataset_id] = dataset
        
        print(f"Found {len(unique_datasets)} unique datasets")
        print()
        
        # Analyze datasets and select best candidates
        selected_datasets = self.select_best_datasets(list(unique_datasets.values()))
        
        # Download selected datasets
        for dataset in selected_datasets[:3]:  # Limit to 3 datasets
            self.download_dataset_sample(dataset)
    
    def select_best_datasets(self, datasets):
        """Select the best datasets for UDT kinematic analysis."""
        print("SELECTING BEST DATASETS FOR UDT ANALYSIS")
        print("-" * 36)
        
        scored_datasets = []
        
        for dataset in datasets:
            score = 0
            title = dataset.get('metadata', {}).get('title', '')
            
            # Scoring criteria for UDT analysis
            if 'Run2016' in title:
                score += 10  # Prefer 2016 data
            if '13TeV' in title:
                score += 8   # Prefer 13 TeV energy
            if 'AOD' in title:
                score += 6   # AOD format good for analysis
            if 'MiniAOD' in title:
                score += 7   # MiniAOD even better
            if 'NanoAOD' in title:
                score += 9   # NanoAOD most compact
            if 'Muon' in title:
                score += 5   # Muons good for kinematics
            if 'Electron' in title or 'EG' in title:
                score += 4   # Electrons also useful
            
            # Get file information
            files = dataset.get('metadata', {}).get('files', [])
            total_size = sum(f.get('size', 0) for f in files)
            
            # Prefer smaller datasets (in GB)
            if total_size < 1e9:  # < 1 GB
                score += 8
            elif total_size < 10e9:  # < 10 GB
                score += 6
            elif total_size < 100e9:  # < 100 GB
                score += 3
            
            scored_datasets.append((score, dataset, total_size))
            
            print(f"Dataset: {title[:60]}...")
            print(f"  Score: {score}, Size: {total_size/1e9:.2f} GB")
            print()
        
        # Sort by score (highest first)
        scored_datasets.sort(key=lambda x: x[0], reverse=True)
        
        print("TOP SELECTED DATASETS:")
        for i, (score, dataset, size) in enumerate(scored_datasets[:5]):
            title = dataset.get('metadata', {}).get('title', '')
            print(f"{i+1}. Score: {score}, Size: {size/1e9:.2f} GB")
            print(f"   {title}")
            print()
        
        return [item[1] for item in scored_datasets]
    
    def download_dataset_sample(self, dataset):
        """Download a sample of a dataset."""
        metadata = dataset.get('metadata', {})
        title = metadata.get('title', 'Unknown')
        files = metadata.get('files', [])
        
        print(f"DOWNLOADING DATASET SAMPLE: {title}")
        print("-" * 60)
        
        if not files:
            print("No files found in dataset")
            return
        
        # Select a small sample of files (max 3 files or 1 GB total)
        selected_files = []
        total_size = 0
        max_size = 1e9  # 1 GB limit
        
        for file_info in files[:10]:  # Check first 10 files
            file_size = file_info.get('size', 0)
            if total_size + file_size <= max_size:
                selected_files.append(file_info)
                total_size += file_size
            if len(selected_files) >= 3:  # Max 3 files
                break
        
        print(f"Selected {len(selected_files)} files, total size: {total_size/1e6:.1f} MB")
        
        # Create subdirectory for this dataset
        dataset_dir = os.path.join(self.download_dir, title.replace('/', '_')[:50])
        os.makedirs(dataset_dir, exist_ok=True)
        
        # Download files
        for i, file_info in enumerate(selected_files):
            self.download_file(file_info, dataset_dir, i+1, len(selected_files))
        
        # Save dataset metadata
        metadata_file = os.path.join(dataset_dir, 'dataset_metadata.json')
        with open(metadata_file, 'w') as f:
            json.dump(dataset, f, indent=2)
        
        print(f"Dataset sample downloaded to: {dataset_dir}")
        print()
    
    def download_file(self, file_info, dataset_dir, file_num, total_files):
        """Download a single file."""
        uri = file_info.get('uri', '')
        filename = file_info.get('filename', f'file_{file_num}')
        file_size = file_info.get('size', 0)
        
        if not uri:
            print(f"No URI for file {filename}")
            return
        
        # Full download URL
        if uri.startswith('http'):
            download_url = uri
        else:
            download_url = urljoin(self.base_url, uri)
        
        local_path = os.path.join(dataset_dir, filename)
        
        print(f"Downloading file {file_num}/{total_files}: {filename}")
        print(f"  Size: {file_size/1e6:.1f} MB")
        print(f"  URL: {download_url}")
        
        try:
            # Check if file already exists
            if os.path.exists(local_path) and os.path.getsize(local_path) == file_size:
                print(f"  File already exists, skipping")
                return
            
            # Download with progress
            response = requests.get(download_url, stream=True, timeout=60)
            response.raise_for_status()
            
            with open(local_path, 'wb') as f:
                downloaded = 0
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        
                        # Progress indicator
                        if file_size > 0:
                            progress = (downloaded / file_size) * 100
                            print(f"  Progress: {progress:.1f}%", end='\r')
            
            print(f"  Downloaded successfully")
            
        except Exception as e:
            print(f"  Error downloading file: {e}")
            if os.path.exists(local_path):
                os.remove(local_path)
    
    def run_download(self):
        """Run the complete download process."""
        print("STARTING CMS OPEN DATA DOWNLOAD")
        print("=" * 32)
        print("Goal: Download manageable CMS datasets for UDT kinematic analysis")
        print("Strategy: Focus on smaller files suitable for testing geometric effects")
        print()
        
        self.download_small_cms_samples()
        
        print("=" * 60)
        print("CMS OPEN DATA DOWNLOAD COMPLETE")
        print("=" * 60)
        print(f"Files downloaded to: {self.download_dir}")
        print("These datasets can be used to test UDT predictions for:")
        print("- Particle momentum distributions under different spacetime geometry")
        print("- Kinematic correlations in multi-particle events")
        print("- Deviations from Standard Model expectations in high-energy collisions")

def main():
    """Main download routine."""
    downloader = CMSOpenDataDownloader()
    downloader.run_download()

if __name__ == "__main__":
    main()