"""
Data Integrity Checker

Utilities for verifying data integrity and detecting synthetic data contamination.
Part of the comprehensive validation framework.
"""

import hashlib
import numpy as np
import pandas as pd
from pathlib import Path


class DataIntegrityChecker:
    """
    Data integrity verification and synthetic data detection.
    
    Ensures all analyses use real observational data only.
    """
    
    def __init__(self):
        """Initialize data integrity checker."""
        self.integrity_log = []
        
    def compute_file_hash(self, file_path):
        """
        Compute SHA-256 hash of a file for integrity verification.
        
        Parameters
        ----------
        file_path : str
            Path to file
            
        Returns
        -------
        str
            SHA-256 hash of file contents
        """
        try:
            with open(file_path, 'rb') as f:
                file_hash = hashlib.sha256(f.read()).hexdigest()
            return file_hash
        except Exception as e:
            return f"ERROR: {str(e)}"
            
    def verify_data_authenticity(self, data_array, data_source_description):
        """
        Check data for patterns indicating synthetic generation.
        
        Parameters
        ----------
        data_array : array-like
            Data to check for synthetic patterns
        data_source_description : str
            Description of data source
            
        Returns
        -------
        dict
            Authenticity assessment
        """
        data = np.array(data_array)
        
        # Statistical tests for synthetic data detection
        authenticity_score = 100  # Start with 100% authentic score
        warnings = []
        
        # Test 1: Perfect mathematical relationships (suspicious)
        if len(data) > 10:
            # Check for suspiciously perfect correlations with simple functions
            x = np.arange(len(data))
            
            # Linear correlation
            linear_corr = abs(np.corrcoef(x, data)[0,1])
            if linear_corr > 0.999:
                authenticity_score -= 30
                warnings.append("Suspiciously perfect linear correlation")
            
            # Polynomial correlation
            if len(data) > 20:
                poly_fit = np.polyfit(x, data, 3)
                poly_pred = np.polyval(poly_fit, x)
                poly_r2 = 1 - np.sum((data - poly_pred)**2) / np.sum((data - np.mean(data))**2)
                if poly_r2 > 0.9999:
                    authenticity_score -= 40
                    warnings.append("Suspiciously perfect polynomial fit")
        
        # Test 2: Unrealistic precision (too many significant figures)
        if len(data) > 5:
            # Check decimal precision
            decimal_places = []
            for val in data:
                if val != 0:
                    str_val = f"{val:.10f}".rstrip('0')
                    if '.' in str_val:
                        decimal_places.append(len(str_val.split('.')[1]))
                    else:
                        decimal_places.append(0)
            
            if len(decimal_places) > 0:
                avg_precision = np.mean(decimal_places)
                if avg_precision > 8:  # More than 8 decimal places is suspicious
                    authenticity_score -= 25
                    warnings.append("Unrealistic numerical precision")
        
        # Test 3: Missing observational uncertainties
        data_std = np.std(data)
        data_range = np.max(data) - np.min(data)
        if data_std / data_range < 0.001:  # Suspiciously low scatter
            authenticity_score -= 20
            warnings.append("Suspiciously low data scatter")
        
        # Test 4: Perfect noise patterns
        if len(data) > 50:
            # Check for artificial noise patterns
            diff_data = np.diff(data)
            diff_std = np.std(diff_data)
            if diff_std < 1e-10:  # No variation in differences
                authenticity_score -= 35
                warnings.append("No realistic observational noise")
        
        assessment = {
            'authenticity_score': max(0, authenticity_score),
            'classification': 'authentic' if authenticity_score > 70 else 'suspicious' if authenticity_score > 30 else 'likely_synthetic',
            'warnings': warnings,
            'data_source': data_source_description,
            'n_data_points': len(data)
        }
        
        return assessment
        
    def create_data_manifest(self, data_directory, output_file='data/manifest_sha256.txt'):
        """
        Create SHA-256 manifest of all data files for integrity checking.
        
        Parameters
        ----------
        data_directory : str
            Directory containing data files
        output_file : str
            Path to output manifest file
            
        Returns
        -------
        dict
            Manifest creation results
        """
        data_path = Path(data_directory)
        
        if not data_path.exists():
            raise ValueError(f"Data directory not found: {data_directory}")
        
        # Find all data files
        data_patterns = ['*.dat', '*.txt', '*.csv', '*.mrt']
        data_files = []
        
        for pattern in data_patterns:
            data_files.extend(data_path.rglob(pattern))
        
        # Compute hashes
        manifest_entries = []
        
        for file_path in sorted(data_files):
            file_hash = self.compute_file_hash(file_path)
            relative_path = file_path.relative_to(data_path)
            
            manifest_entries.append(f"{file_hash}  {relative_path}")
        
        # Write manifest
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            f.write("# SHA-256 Data Integrity Manifest\\n")
            f.write("# Generated by UDT Data Integrity Checker\\n")
            f.write("# Format: <hash>  <file_path>\\n")
            f.write("\\n")
            for entry in manifest_entries:
                f.write(entry + "\\n")
        
        manifest_results = {
            'manifest_file': str(output_path),
            'n_files_processed': len(data_files),
            'data_directory': str(data_path),
            'status': 'created'
        }
        
        print(f"Data manifest created: {output_path}")
        print(f"Files processed: {len(data_files)}")
        
        return manifest_results
        
    def verify_data_manifest(self, manifest_file='data/manifest_sha256.txt'):
        """
        Verify data integrity against SHA-256 manifest.
        
        Parameters
        ----------
        manifest_file : str
            Path to manifest file
            
        Returns
        -------
        dict
            Verification results
        """
        manifest_path = Path(manifest_file)
        
        if not manifest_path.exists():
            raise ValueError(f"Manifest file not found: {manifest_file}")
        
        # Read manifest
        with open(manifest_path, 'r') as f:
            lines = f.readlines()
        
        # Parse manifest entries
        verification_results = []
        
        for line in lines:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
                
            parts = line.split(None, 1)
            if len(parts) != 2:
                continue
                
            expected_hash, file_path = parts
            
            # Compute actual hash
            full_path = manifest_path.parent / file_path
            actual_hash = self.compute_file_hash(full_path)
            
            verification_results.append({
                'file_path': file_path,
                'expected_hash': expected_hash,
                'actual_hash': actual_hash,
                'verified': expected_hash == actual_hash
            })
        
        # Summary
        total_files = len(verification_results)
        verified_files = sum(1 for r in verification_results if r['verified'])
        
        summary = {
            'manifest_file': str(manifest_path),
            'total_files': total_files,
            'verified_files': verified_files,
            'failed_files': total_files - verified_files,
            'integrity_status': 'PASSED' if verified_files == total_files else 'FAILED',
            'detailed_results': verification_results
        }
        
        print(f"Data integrity verification: {summary['integrity_status']}")
        print(f"Files verified: {verified_files}/{total_files}")
        
        if summary['failed_files'] > 0:
            print("INTEGRITY FAILURE - Some files have been modified!")
            
        return summary
        
    def check_for_synthetic_data_in_script(self, script_path):
        """
        Scan a script for synthetic data generation functions.
        
        Parameters
        ----------
        script_path : str
            Path to Python script to check
            
        Returns
        -------
        dict
            Synthetic data detection results
        """
        try:
            with open(script_path, 'r', encoding='utf-8') as f:
                script_content = f.read()
        except Exception as e:
            return {'status': 'error', 'error': str(e)}
        
        # Patterns that indicate synthetic data generation
        synthetic_patterns = [
            'np.random',
            'random.normal',
            'generate_synthetic',
            'artificial_data',
            'fake_data',
            'simulate_',
            'mock_data',
            'toy_model'
        ]
        
        detected_patterns = []
        for pattern in synthetic_patterns:
            if pattern in script_content:
                detected_patterns.append(pattern)
        
        result = {
            'script_path': script_path,
            'synthetic_patterns_detected': detected_patterns,
            'contains_synthetic_data': len(detected_patterns) > 0,
            'status': 'clean' if len(detected_patterns) == 0 else 'contains_synthetic_code'
        }
        
        return result