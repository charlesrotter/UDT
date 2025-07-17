#!/usr/bin/env python3
"""
Cleanup Script for Experimental Files
=====================================

This script helps clean up experimental files before publication.
Run with --dry-run first to see what would be deleted.

Usage:
    python scripts/cleanup_experimental.py --dry-run
    python scripts/cleanup_experimental.py --confirm

Author: Charles Rotter
Date: 2025-01-17
"""

import os
import sys
import glob
import shutil
from datetime import datetime
import argparse

def find_experimental_files():
    """Find all experimental files based on naming patterns."""
    patterns = [
        "experimental_*.py",
        "draft_*.py", 
        "test_idea_*.py",
        "failed_*.py",
        "musing_*.py",
        "**/*.tmp",
        "**/*.bak",
        "**/*~"
    ]
    
    experimental_files = []
    for pattern in patterns:
        experimental_files.extend(glob.glob(pattern, recursive=True))
        experimental_files.extend(glob.glob(f"**/{pattern}", recursive=True))
    
    # Remove duplicates and sort
    experimental_files = sorted(set(experimental_files))
    
    return experimental_files

def analyze_exploration_directory():
    """Analyze contents of exploration directory."""
    exploration_path = "exploration"
    
    if not os.path.exists(exploration_path):
        return None
    
    stats = {
        'total_files': 0,
        'total_size': 0,
        'by_subdirectory': {}
    }
    
    for root, dirs, files in os.walk(exploration_path):
        for file in files:
            filepath = os.path.join(root, file)
            size = os.path.getsize(filepath)
            stats['total_files'] += 1
            stats['total_size'] += size
            
            # Track by subdirectory
            subdir = root.replace(exploration_path, '').strip(os.sep).split(os.sep)[0]
            if subdir:
                if subdir not in stats['by_subdirectory']:
                    stats['by_subdirectory'][subdir] = {'files': 0, 'size': 0}
                stats['by_subdirectory'][subdir]['files'] += 1
                stats['by_subdirectory'][subdir]['size'] += size
    
    return stats

def cleanup_files(files_to_remove, dry_run=True):
    """Remove experimental files."""
    print(f"\n{'DRY RUN - ' if dry_run else ''}Cleaning up experimental files...")
    print("=" * 60)
    
    if not files_to_remove:
        print("No experimental files found to remove.")
        return
    
    print(f"Found {len(files_to_remove)} files to remove:\n")
    
    for file in files_to_remove:
        if os.path.exists(file):
            size = os.path.getsize(file)
            print(f"  {file} ({size:,} bytes)")
            
            if not dry_run:
                try:
                    os.remove(file)
                    print(f"    ✓ Removed")
                except Exception as e:
                    print(f"    ✗ Error: {e}")
    
    if dry_run:
        print(f"\nTotal files that would be removed: {len(files_to_remove)}")
        print("Run with --confirm to actually delete these files.")

def archive_exploration(dry_run=True):
    """Archive exploration directory."""
    exploration_path = "exploration"
    
    if not os.path.exists(exploration_path):
        print("\nNo exploration directory found to archive.")
        return
    
    stats = analyze_exploration_directory()
    if stats:
        print(f"\n{'DRY RUN - ' if dry_run else ''}Exploration directory analysis:")
        print("=" * 60)
        print(f"Total files: {stats['total_files']}")
        print(f"Total size: {stats['total_size']:,} bytes ({stats['total_size']/1024/1024:.2f} MB)")
        print("\nBy subdirectory:")
        for subdir, info in stats['by_subdirectory'].items():
            print(f"  {subdir}/: {info['files']} files, {info['size']:,} bytes")
    
    if not dry_run:
        # Create archive
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        archive_name = f"exploration_archive_{timestamp}.tar.gz"
        
        print(f"\nCreating archive: {archive_name}")
        try:
            import tarfile
            with tarfile.open(archive_name, "w:gz") as tar:
                tar.add(exploration_path)
            print(f"✓ Archive created successfully")
            
            # Ask for confirmation before removing directory
            response = input("\nRemove exploration directory? (yes/no): ")
            if response.lower() == 'yes':
                shutil.rmtree(exploration_path)
                print("✓ Exploration directory removed")
            else:
                print("✗ Exploration directory kept")
        except Exception as e:
            print(f"✗ Error creating archive: {e}")
    else:
        print("\nRun with --confirm to create archive and optionally remove directory.")

def main():
    parser = argparse.ArgumentParser(description="Clean up experimental files")
    parser.add_argument('--dry-run', action='store_true', 
                        help='Show what would be deleted without actually deleting')
    parser.add_argument('--confirm', action='store_true',
                        help='Actually perform the cleanup')
    parser.add_argument('--skip-exploration', action='store_true',
                        help='Skip exploration directory archiving')
    
    args = parser.parse_args()
    
    if not args.dry_run and not args.confirm:
        print("Error: Must specify either --dry-run or --confirm")
        sys.exit(1)
    
    if args.dry_run and args.confirm:
        print("Error: Cannot specify both --dry-run and --confirm")
        sys.exit(1)
    
    dry_run = args.dry_run
    
    print("UDT Experimental Files Cleanup Tool")
    print("===================================\n")
    
    # Find and clean experimental files
    experimental_files = find_experimental_files()
    cleanup_files(experimental_files, dry_run=dry_run)
    
    # Archive exploration directory
    if not args.skip_exploration:
        archive_exploration(dry_run=dry_run)
    
    print("\n✓ Cleanup complete!" if not dry_run else "\n✓ Dry run complete!")

if __name__ == "__main__":
    main()