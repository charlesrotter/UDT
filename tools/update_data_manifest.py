#!/usr/bin/env python3
"""
Update data integrity manifest (SHA-256 checksums) for all data files.

This script generates a deterministic manifest of all files under the data/ directory,
excluding the manifest itself. The output is sorted for reproducibility.

Usage:
    python tools/update_data_manifest.py
    
The manifest will be written to data/manifest_sha256.txt
"""

import hashlib
import os
from pathlib import Path
import sys


def calculate_sha256(filepath):
    """Calculate SHA-256 hash of a file."""
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        # Read in chunks to handle large files efficiently
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()


def generate_manifest(data_dir, output_file):
    """Generate SHA-256 manifest for all files in data directory."""
    data_path = Path(data_dir)
    manifest_path = Path(output_file)
    
    if not data_path.exists():
        print(f"Error: Data directory '{data_dir}' does not exist.")
        sys.exit(1)
    
    # Collect all files (excluding manifest itself and hidden files)
    files_to_hash = []
    for file_path in data_path.rglob("*"):
        if file_path.is_file() and file_path != manifest_path:
            # Skip hidden files and __pycache__
            if not any(part.startswith('.') for part in file_path.parts) and '__pycache__' not in str(file_path):
                files_to_hash.append(file_path)
    
    # Sort files for deterministic output
    files_to_hash.sort()
    
    # Generate manifest
    manifest_lines = []
    total_files = len(files_to_hash)
    
    print(f"Generating SHA-256 manifest for {total_files} files...")
    
    for i, file_path in enumerate(files_to_hash, 1):
        try:
            # Calculate hash
            file_hash = calculate_sha256(file_path)
            
            # Use forward slashes for cross-platform compatibility
            relative_path = file_path.relative_to(Path.cwd()).as_posix()
            
            manifest_lines.append(f"{file_hash}  {relative_path}")
            
            # Progress indicator
            if i % 50 == 0 or i == total_files:
                print(f"  Processed {i}/{total_files} files...")
                
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue
    
    # Write manifest
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(manifest_lines))
        if manifest_lines:  # Add final newline if there are entries
            f.write('\n')
    
    print(f"\nManifest written to: {manifest_path}")
    print(f"Total files: {len(manifest_lines)}")
    print(f"Manifest size: {manifest_path.stat().st_size / 1024:.1f} KB")


def verify_manifest(manifest_file):
    """Verify files against the manifest."""
    manifest_path = Path(manifest_file)
    
    if not manifest_path.exists():
        print(f"Error: Manifest file '{manifest_file}' does not exist.")
        sys.exit(1)
    
    print(f"Verifying files against manifest: {manifest_file}")
    
    errors = 0
    checked = 0
    
    with open(manifest_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            parts = line.split('  ', 1)
            if len(parts) != 2:
                print(f"Error: Invalid manifest line: {line}")
                errors += 1
                continue
            
            expected_hash, file_path = parts
            file_path = Path(file_path)
            
            if not file_path.exists():
                print(f"Missing: {file_path}")
                errors += 1
                continue
            
            actual_hash = calculate_sha256(file_path)
            if actual_hash != expected_hash:
                print(f"Mismatch: {file_path}")
                print(f"  Expected: {expected_hash}")
                print(f"  Actual:   {actual_hash}")
                errors += 1
            
            checked += 1
            if checked % 50 == 0:
                print(f"  Verified {checked} files...")
    
    print(f"\nVerification complete: {checked} files checked, {errors} errors")
    return errors == 0


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Update or verify data integrity manifest")
    parser.add_argument('--verify', action='store_true', 
                        help='Verify existing manifest instead of updating')
    parser.add_argument('--data-dir', default='data', 
                        help='Data directory to scan (default: data)')
    parser.add_argument('--output', default='data/manifest_sha256.txt',
                        help='Output manifest file (default: data/manifest_sha256.txt)')
    
    args = parser.parse_args()
    
    if args.verify:
        success = verify_manifest(args.output)
        sys.exit(0 if success else 1)
    else:
        generate_manifest(args.data_dir, args.output)