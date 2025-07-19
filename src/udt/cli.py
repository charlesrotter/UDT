"""
UDT Command Line Interface

Entry point for UDT validation and analysis commands.
"""

import argparse
import sys
from pathlib import Path

from .diagnostics.validation import ValidationSuite
from .diagnostics.integrity import DataIntegrityChecker


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description='Universal Distance Dilation Theory Validation Suite',
        prog='udt-validate'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Validation command
    validate_parser = subparsers.add_parser('validate', help='Run UDT validation suite')
    validate_parser.add_argument('--sparc-data', type=str, help='Path to SPARC data directory')
    validate_parser.add_argument('--supernova-data', type=str, help='Path to supernova data directory')
    validate_parser.add_argument('--cmb-data', type=str, help='Path to CMB data directory')
    validate_parser.add_argument('--max-galaxies', type=int, default=50, help='Maximum galaxies to analyze')
    
    # Integrity command
    integrity_parser = subparsers.add_parser('integrity', help='Check data integrity')
    integrity_parser.add_argument('--data-dir', type=str, required=True, help='Data directory to check')
    integrity_parser.add_argument('--create-manifest', action='store_true', help='Create integrity manifest')
    integrity_parser.add_argument('--verify-manifest', action='store_true', help='Verify against manifest')
    
    args = parser.parse_args()
    
    if args.command == 'validate':
        # Run validation suite
        suite = ValidationSuite()
        
        data_directories = {}
        if args.sparc_data:
            data_directories['sparc'] = args.sparc_data
        if args.supernova_data:
            data_directories['supernova'] = args.supernova_data
        if args.cmb_data:
            data_directories['cmb'] = args.cmb_data
            
        if not data_directories:
            print("Error: No data directories specified")
            return 1
            
        results = suite.run_comprehensive_validation(data_directories, args.max_galaxies)
        
        # Print summary
        print("\\nValidation Summary:")
        for key, result in results.items():
            if isinstance(result, dict) and 'status' in result:
                print(f"  {key}: {result['status']}")
        
        return 0
        
    elif args.command == 'integrity':
        # Run integrity checks
        checker = DataIntegrityChecker()
        
        if args.create_manifest:
            results = checker.create_data_manifest(args.data_dir)
            print(f"Manifest created: {results['status']}")
            
        if args.verify_manifest:
            results = checker.verify_data_manifest()
            print(f"Integrity check: {results['integrity_status']}")
            return 0 if results['integrity_status'] == 'PASSED' else 1
            
        return 0
        
    else:
        parser.print_help()
        return 1


if __name__ == '__main__':
    sys.exit(main())