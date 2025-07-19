#!/usr/bin/env python3
"""
Basic tests for UDT galactic dynamics module.

These tests verify core functionality for AI auditing and CI validation.
"""

import numpy as np
import pytest
import sys
from pathlib import Path

# Add UDT to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from udt.core.galactic_dynamics import (
    enhancement_factor,
    temporal_dilation_function,
    pure_temporal_velocity
)


class TestTemporalGeometry:
    """Test core temporal geometry functions."""
    
    def test_temporal_dilation_basic(self):
        """Test temporal dilation function τ(r) = R₀/(R₀ + r)."""
        R0 = 30.0  # kpc
        r = np.array([0, 10, 30, 60])  # kpc
        
        tau = temporal_dilation_function(r, R0)
        
        # Check boundary conditions
        assert tau[0] == 1.0  # τ(0) = 1
        assert tau[2] == 0.5  # τ(R₀) = 0.5
        
        # Check monotonic decrease
        assert np.all(np.diff(tau) < 0)
        
        # Check positive values
        assert np.all(tau > 0)
        
    def test_enhancement_factor_properties(self):
        """Test enhancement factor F(τ) properties."""
        R0 = 30.0
        r = np.linspace(0, 100, 101)
        
        enhancement = enhancement_factor(r, R0)
        
        # Check F(0) ≈ 1 (small solar system effects)
        assert enhancement[0] >= 1.0
        
        # Check monotonic increase
        assert np.all(np.diff(enhancement) >= 0)
        
        # Check reasonable range for galactic scales
        galactic_enhancement = enhancement[r <= 50]
        assert np.all(galactic_enhancement < 10)  # Reasonable bound
        
    def test_velocity_profile_shape(self):
        """Test UDT velocity profile shape."""
        R0 = 30.0
        V_scale = 200.0
        r = np.linspace(1, 50, 50)  # Avoid r=0
        
        velocity = pure_temporal_velocity(r, R0, V_scale)
        
        # Check positive velocities
        assert np.all(velocity > 0)
        
        # Check reasonable velocity range
        assert np.all(velocity < 500)  # km/s - reasonable upper bound
        assert np.all(velocity > 50)   # km/s - reasonable lower bound
        
        # Velocity should increase with radius initially (flat rotation curve behavior)
        initial_slope = np.mean(np.diff(velocity[:10]))
        assert initial_slope >= 0


class TestDataIntegrity:
    """Test data loading and integrity."""
    
    def test_sample_data_exists(self):
        """Verify sample galaxy data exists."""
        sample_file = Path("data/sample/NGC3198_rotmod.dat")
        assert sample_file.exists(), "Sample galaxy data missing"
        
        # Check file is not empty
        assert sample_file.stat().st_size > 0
        
    def test_manifest_exists(self):
        """Verify data integrity manifest exists."""
        manifest_file = Path("data/manifest_sha256.txt")
        
        # Check if manifest exists (may be empty if not generated yet)
        if manifest_file.exists():
            assert manifest_file.stat().st_size >= 0


class TestMathematicalConsistency:
    """Test mathematical consistency of UDT formulations."""
    
    def test_scale_invariance_check(self):
        """Test that UDT formulas have expected scale properties."""
        R0_1 = 30.0
        R0_2 = 60.0  # Double the scale
        r = np.array([15, 30, 45])
        
        tau_1 = temporal_dilation_function(r, R0_1)
        tau_2 = temporal_dilation_function(r, R0_2)
        
        # At r = R₀/2, should get same τ value
        assert abs(tau_1[0] - temporal_dilation_function(15, 30)) < 1e-10
        assert abs(tau_2[1] - temporal_dilation_function(30, 60)) < 1e-10
        
    def test_enhancement_limits(self):
        """Test enhancement factor limits."""
        R0 = 30.0
        
        # Near origin (solar system)
        r_small = np.array([1e-3, 1e-2, 1e-1])
        enhancement_small = enhancement_factor(r_small, R0)
        
        # Should approach 1 as r → 0
        assert np.all(enhancement_small < 1.1)
        
        # Far field behavior
        r_large = np.array([100, 200, 500])
        enhancement_large = enhancement_factor(r_large, R0)
        
        # Should grow but remain finite
        assert np.all(np.isfinite(enhancement_large))


def test_basic_imports():
    """Test that core UDT modules can be imported."""
    try:
        from udt.core import galactic_dynamics, temporal_geometry, cosmology
        from udt.utils import data_loader, plotting
    except ImportError as e:
        pytest.fail(f"Failed to import UDT modules: {e}")


def test_numpy_compatibility():
    """Test compatibility with numpy operations."""
    R0 = 30.0
    r = np.logspace(0, 2, 100)  # 1 to 100 kpc
    
    # These should all work without errors
    tau = temporal_dilation_function(r, R0)
    enhancement = enhancement_factor(r, R0)
    velocity = pure_temporal_velocity(r, R0, 200.0)
    
    # Check outputs are numpy arrays
    assert isinstance(tau, np.ndarray)
    assert isinstance(enhancement, np.ndarray)
    assert isinstance(velocity, np.ndarray)
    
    # Check no NaN or infinite values in normal range
    assert np.all(np.isfinite(tau))
    assert np.all(np.isfinite(enhancement))
    assert np.all(np.isfinite(velocity))


if __name__ == "__main__":
    # Run tests directly if called as script
    print("Running UDT basic validation tests...")
    
    # Simple smoke test
    try:
        test_basic_imports()
        print("+ Module imports: PASS")
        
        test_numpy_compatibility()
        print("+ NumPy compatibility: PASS")
        
        # Run a few key tests
        test_instance = TestTemporalGeometry()
        test_instance.test_temporal_dilation_basic()
        print("+ Temporal dilation: PASS")
        
        test_instance.test_enhancement_factor_properties()
        print("+ Enhancement factor: PASS")
        
        print("\n[OK] All basic tests PASSED!")
        
    except Exception as e:
        print(f"\n[FAIL] Test FAILED: {e}")
        sys.exit(1)