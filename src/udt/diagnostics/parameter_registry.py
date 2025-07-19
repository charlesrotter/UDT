"""
Parameter Registry System

Prevents parameter regression errors by maintaining validated parameter sets
for different analysis types and automatically checking parameter consistency.
"""

import json
import os
from pathlib import Path
from datetime import datetime
import warnings


class ParameterRegistry:
    """
    Central registry for validated UDT parameters across different scales and analyses.
    
    Prevents regression errors by:
    1. Storing validated parameter sets for each analysis type
    2. Warning when parameters deviate from validated values
    3. Providing scale-appropriate parameters automatically
    4. Tracking parameter evolution and validation history
    """
    
    def __init__(self, registry_file=None):
        """Initialize parameter registry."""
        if registry_file is None:
            # Default location in project
            self.registry_file = Path(__file__).parent.parent.parent.parent / "parameter_registry.json"
        else:
            self.registry_file = Path(registry_file)
        
        self.load_registry()
    
    def load_registry(self):
        """Load parameter registry from file."""
        if self.registry_file.exists():
            try:
                with open(self.registry_file, 'r') as f:
                    self.registry = json.load(f)
            except Exception as e:
                warnings.warn(f"Could not load parameter registry: {e}")
                self.registry = self._create_default_registry()
        else:
            self.registry = self._create_default_registry()
            self.save_registry()
    
    def save_registry(self):
        """Save parameter registry to file."""
        try:
            with open(self.registry_file, 'w') as f:
                json.dump(self.registry, f, indent=2)
        except Exception as e:
            warnings.warn(f"Could not save parameter registry: {e}")
    
    def _create_default_registry(self):
        """Create default parameter registry with validated values."""
        return {
            "metadata": {
                "created": datetime.now().isoformat(),
                "version": "1.0",
                "description": "UDT validated parameter registry"
            },
            "multi_scale_framework": {
                "galactic": {
                    "R0_mpc": 0.038,
                    "R0_kpc": 38.0,
                    "description": "Galaxy rotation curve scale",
                    "validated_analyses": ["sparc_rotation_curves"],
                    "typical_distances": "0.01-100 kpc",
                    "last_validated": None,
                    "chi2_dof_best": None
                },
                "cosmological": {
                    "R0_mpc": 3000.0,
                    "R0_gpc": 3.0,
                    "description": "Supernova distance scale",
                    "validated_analyses": ["pantheon_supernovae"],
                    "typical_distances": "10-10000 Mpc",
                    "last_validated": None,
                    "chi2_dof_best": None
                },
                "cmb": {
                    "R0_mpc": 13041.1,
                    "R0_gpc": 13.041,
                    "description": "CMB acoustic peak scale",
                    "validated_analyses": ["planck_cmb_analysis"],
                    "typical_distances": "14000 Mpc (recombination)",
                    "last_validated": None,
                    "chi2_dof_best": 27.32,
                    "alternative_method": {
                        "R0_mpc": 10316.4,
                        "description": "Optimized for l1=220",
                        "method": "angular_diameter_distance_modification"
                    }
                }
            },
            "analysis_specific_parameters": {
                "cmb_analysis": {
                    "standard_cosmology_baseline": {
                        "eta_rec_mpc_per_c": 288.0,
                        "r_s_mpc": 147.3,
                        "D_A_standard_mpc": 13975.0,
                        "l1_target": 220.0
                    },
                    "udt_modifications": {
                        "method": "angular_diameter_distance",
                        "formula": "D_A = eta_rec * tau(eta_rec)",
                        "tau_formula": "R0_cmb / (R0_cmb + eta_rec)"
                    }
                },
                "supernova_analysis": {
                    "artifact_correction": {
                        "systematic_uncertainty": 0.03,
                        "bias_correction_required": True,
                        "mean_correction_magnitude": 0.068
                    }
                },
                "galaxy_analysis": {
                    "enhancement_formula": "1/tau^2 = (1 + r/R0)^2",
                    "typical_enhancement_range": [1.1, 100.0]
                }
            },
            "validation_history": {
                "successful_analyses": [],
                "failed_analyses": [],
                "parameter_updates": []
            }
        }
    
    def get_parameters_for_analysis(self, analysis_type, scale=None):
        """
        Get validated parameters for a specific analysis type.
        
        Parameters
        ----------
        analysis_type : str
            Type of analysis ('cmb', 'supernova', 'galaxy', 'ligo', 'muon_g2', 'bao')
        scale : str, optional
            Specific scale ('galactic', 'cosmological', 'cmb')
            
        Returns
        -------
        dict
            Validated parameters for the analysis
        """
        
        # Map analysis types to scales
        scale_mapping = {
            'galaxy': 'galactic',
            'sparc': 'galactic',
            'supernova': 'cosmological',
            'pantheon': 'cosmological',
            'cmb': 'cmb',
            'planck': 'cmb',
            'ligo': 'galactic',  # Local scale for detector separation
            'muon_g2': 'galactic',  # Particle scale ~ galactic
            'bao': 'cosmological'  # Large scale structure
        }
        
        # Determine appropriate scale
        if scale is None:
            scale = scale_mapping.get(analysis_type.lower())
            if scale is None:
                raise ValueError(f"Unknown analysis type: {analysis_type}")
        
        # Get multi-scale parameters
        if scale in self.registry["multi_scale_framework"]:
            base_params = self.registry["multi_scale_framework"][scale].copy()
        else:
            raise ValueError(f"Unknown scale: {scale}")
        
        # Add analysis-specific parameters
        analysis_key = f"{analysis_type.lower()}_analysis"
        if analysis_key in self.registry["analysis_specific_parameters"]:
            base_params.update(self.registry["analysis_specific_parameters"][analysis_key])
        
        # Add metadata
        base_params["analysis_type"] = analysis_type
        base_params["scale"] = scale
        base_params["source"] = "validated_parameter_registry"
        
        return base_params
    
    def validate_parameters(self, analysis_type, provided_params):
        """
        Validate provided parameters against registry.
        
        Parameters
        ----------
        analysis_type : str
            Type of analysis
        provided_params : dict
            Parameters being used
            
        Returns
        -------
        dict
            Validation results with warnings/errors
        """
        validated_params = self.get_parameters_for_analysis(analysis_type)
        
        validation_results = {
            "status": "valid",
            "warnings": [],
            "errors": [],
            "recommendations": []
        }
        
        # Check R0 parameter
        if "R0_mpc" in provided_params and "R0_mpc" in validated_params:
            provided_R0 = provided_params["R0_mpc"]
            validated_R0 = validated_params["R0_mpc"]
            
            # Check for significant deviation
            if abs(provided_R0 - validated_R0) / validated_R0 > 0.1:  # >10% deviation
                validation_results["warnings"].append(
                    f"R0 parameter deviation: provided={provided_R0}, validated={validated_R0}"
                )
                validation_results["recommendations"].append(
                    f"Consider using validated R0={validated_R0} Mpc for {analysis_type}"
                )
            
            # Check for wrong scale usage
            if analysis_type.lower() == 'cmb' and provided_R0 < 1000:
                validation_results["errors"].append(
                    f"CRITICAL: Using galactic-scale R0={provided_R0} for CMB analysis"
                )
                validation_results["recommendations"].append(
                    f"Use CMB-scale R0={validated_R0} Mpc instead"
                )
                validation_results["status"] = "error"
            
            elif analysis_type.lower() in ['galaxy', 'sparc'] and provided_R0 > 1000:
                validation_results["warnings"].append(
                    f"Using large R0={provided_R0} for galactic analysis"
                )
                validation_results["recommendations"].append(
                    f"Consider galactic-scale R0={validated_R0} Mpc"
                )
        
        # Check for missing critical parameters
        critical_params = ["R0_mpc"]
        for param in critical_params:
            if param not in provided_params:
                validation_results["warnings"].append(f"Missing parameter: {param}")
        
        return validation_results
    
    def update_validation_results(self, analysis_type, parameters, results):
        """Update registry with new validation results."""
        
        # Create validation entry
        validation_entry = {
            "timestamp": datetime.now().isoformat(),
            "analysis_type": analysis_type,
            "parameters": parameters,
            "results": results,
            "chi2_dof": results.get("chi2_dof"),
            "success": results.get("status") == "success"
        }
        
        # Add to history
        if "validation_history" not in self.registry:
            self.registry["validation_history"] = {"successful_analyses": [], "failed_analyses": []}
        
        if validation_entry["success"]:
            self.registry["validation_history"]["successful_analyses"].append(validation_entry)
            
            # Update best result if better
            scale = self._get_scale_for_analysis(analysis_type)
            if scale and "chi2_dof" in results:
                current_best = self.registry["multi_scale_framework"][scale].get("chi2_dof_best")
                if current_best is None or results["chi2_dof"] < current_best:
                    self.registry["multi_scale_framework"][scale]["chi2_dof_best"] = results["chi2_dof"]
                    self.registry["multi_scale_framework"][scale]["last_validated"] = datetime.now().isoformat()
        else:
            self.registry["validation_history"]["failed_analyses"].append(validation_entry)
        
        self.save_registry()
    
    def _get_scale_for_analysis(self, analysis_type):
        """Get scale for analysis type."""
        scale_mapping = {
            'galaxy': 'galactic', 'sparc': 'galactic',
            'supernova': 'cosmological', 'pantheon': 'cosmological',
            'cmb': 'cmb', 'planck': 'cmb'
        }
        return scale_mapping.get(analysis_type.lower())
    
    def get_validation_summary(self):
        """Get summary of validation history."""
        history = self.registry.get("validation_history", {})
        
        summary = {
            "total_successful": len(history.get("successful_analyses", [])),
            "total_failed": len(history.get("failed_analyses", [])),
            "scales": {}
        }
        
        for scale, params in self.registry["multi_scale_framework"].items():
            summary["scales"][scale] = {
                "R0_mpc": params["R0_mpc"],
                "best_chi2_dof": params.get("chi2_dof_best"),
                "last_validated": params.get("last_validated"),
                "description": params["description"]
            }
        
        return summary


# Global parameter registry instance
parameter_registry = ParameterRegistry()


def get_validated_parameters(analysis_type, scale=None):
    """
    Get validated parameters for an analysis.
    
    Usage:
    params = get_validated_parameters('cmb')
    R0_cmb = params['R0_mpc']  # Returns 13041.1 Mpc
    """
    return parameter_registry.get_parameters_for_analysis(analysis_type, scale)


def validate_analysis_parameters(analysis_type, parameters):
    """
    Validate parameters before analysis.
    
    Usage:
    validation = validate_analysis_parameters('cmb', {'R0_mpc': 3582})
    if validation['status'] == 'error':
        print("Parameter errors:", validation['errors'])
    """
    return parameter_registry.validate_parameters(analysis_type, parameters)


def update_analysis_results(analysis_type, parameters, results):
    """Update registry with analysis results."""
    return parameter_registry.update_validation_results(analysis_type, parameters, results)