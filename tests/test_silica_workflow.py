"""Comprehensive tests for silica fitting workflows."""

import numpy as np
import pytest
from numpy.typing import NDArray

from lib.state_manager import StateManager
from workflows.silica_workflow import (
    FittingResult,
    SilicaCAFittingParams,
    SilicaCAWorkflow,
)


@pytest.fixture
def state_manager():
    """Create fresh state manager for each test."""
    return StateManager()


@pytest.fixture
def sample_positions() -> NDArray:
    """Sample Z-scan positions in mm."""
    return np.linspace(-2, 2, 41)


@pytest.fixture
def sample_ca_data(sample_positions) -> NDArray:
    """Realistic CA data with Gaussian-like peak (closed aperture signature)."""
    # Peak-valley signature typical of CA
    z = sample_positions
    return 1.0 + 0.15 * np.exp(-(z**2) / 0.3)


class TestSilicaCAWorkflow:
    """Test silica closed aperture fitting workflow."""

    def test_manual_fit_returns_fitting_result(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that manual fit returns valid FittingResult."""
        workflow = SilicaCAWorkflow(state_manager, window_provider=None)

        params = SilicaCAFittingParams(
            dphi0=0.3,
            beamwaist=10e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-20,  # Silica's n2
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert isinstance(result, FittingResult)
        assert result.curve is not None
        assert len(result.curve) == len(sample_positions)
        assert result.converged is True
        assert result.message == "Manual CA fit successful"

    def test_manual_fit_curve_shape(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that fitted curve has reasonable shape."""
        workflow = SilicaCAWorkflow(state_manager, window_provider=None)

        params = SilicaCAFittingParams(
            dphi0=0.25,
            beamwaist=15e-6,
            zero_level=0.95,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert result.converged
        # Curve should be finite
        assert np.all(np.isfinite(result.curve))
        # Curve should have similar shape to input data (roughly)
        assert len(result.curve) == len(sample_ca_data)

    def test_invalid_dphi0_negative_fails_for_silica(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that negative DPhi0 is rejected for silica (n2 > 0 in VIS/NIR)."""
        workflow = SilicaCAWorkflow(state_manager, window_provider=None)

        params = SilicaCAFittingParams(
            dphi0=-0.3,  # Invalid for silica (positive n2 in VIS/NIR)
            beamwaist=10e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-20,  # Positive n2 for silica
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        # Should fail - negative DPhi0 is invalid for silica with positive n2
        assert result.converged is False
        assert "DPhi0" in result.message or "dphi0" in result.message.lower()

    def test_invalid_dphi0_unrealistic_magnitude_fails(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that unrealistically large DPhi0 magnitude is caught."""
        workflow = SilicaCAWorkflow(state_manager, window_provider=None)

        params = SilicaCAFittingParams(
            dphi0=500.0,  # Unrealistically large magnitude
            beamwaist=10e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert result.converged is False
        assert "DPhi0" in result.message or "dphi0" in result.message.lower()

    def test_invalid_beamwaist_too_large_fails(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that invalid (too large) beamwaist is caught."""
        workflow = SilicaCAWorkflow(state_manager, window_provider=None)

        params = SilicaCAFittingParams(
            dphi0=0.3,
            beamwaist=0.5,  # Way too large (500 mm)
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert result.converged is False
        assert "beamwaist" in result.message.lower()

    def test_invalid_beamwaist_too_small_fails(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that unrealistically small beamwaist is caught."""
        workflow = SilicaCAWorkflow(state_manager, window_provider=None)

        params = SilicaCAFittingParams(
            dphi0=0.3,
            beamwaist=1e-9,  # 1 nanometer - unrealistic
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert result.converged is False

    def test_invalid_zero_level_too_low_fails(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that invalid zero level is caught."""
        workflow = SilicaCAWorkflow(state_manager, window_provider=None)

        params = SilicaCAFittingParams(
            dphi0=0.3,
            beamwaist=10e-6,
            zero_level=0.2,  # Too low
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert result.converged is False

    def test_invalid_zero_level_too_high_fails(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that unrealistically high zero level is caught."""
        workflow = SilicaCAWorkflow(state_manager, window_provider=None)

        params = SilicaCAFittingParams(
            dphi0=0.3,
            beamwaist=10e-6,
            zero_level=10.0,  # Unrealistically high
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert result.converged is False

    def test_automatic_fit_returns_result_with_minimizer(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that automatic fit returns FittingResult with minimizer."""
        workflow = SilicaCAWorkflow(state_manager, window_provider=None)

        params = SilicaCAFittingParams(
            dphi0=0.25,
            beamwaist=12e-6,
            zero_level=0.98,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_automatic(params)

        assert isinstance(result, FittingResult)
        assert result.converged is True
        assert result.minimizer_result is not None
        assert "dphi0" in result.params
        assert "beamwaist" in result.params
        assert "center" in result.params
        assert "zero_level" in result.params

    def test_automatic_fit_returns_errors(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that automatic fit includes error estimation."""
        workflow = SilicaCAWorkflow(state_manager, window_provider=None)

        params = SilicaCAFittingParams(
            dphi0=0.25,
            beamwaist=12e-6,
            zero_level=0.98,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_automatic(params)

        assert result.converged
        assert "dphi0" in result.errors
        assert "beamwaist" in result.errors
        assert "center" in result.errors
        assert "zero_level" in result.errors

    def test_state_updated_on_manual_fit(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that state is updated correctly after manual fit."""
        workflow = SilicaCAWorkflow(state_manager, window_provider=None)

        # Initially should be False
        assert state_manager.silica_autofit_done is False

        params = SilicaCAFittingParams(
            dphi0=0.3,
            beamwaist=10e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert result.converged
        # Manual fit sets flag to False (not automatic)
        assert state_manager.silica_autofit_done is False

    def test_state_updated_on_automatic_fit(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that state is updated correctly after automatic fit."""
        workflow = SilicaCAWorkflow(state_manager, window_provider=None)

        assert state_manager.silica_autofit_done is False

        params = SilicaCAFittingParams(
            dphi0=0.25,
            beamwaist=12e-6,
            zero_level=0.98,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_automatic(params)

        assert result.converged
        # Automatic fit sets flag to True
        assert state_manager.silica_autofit_done is True

    def test_fit_with_window_provider_callback(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that workflow works with window_provider callback (UI mode)."""
        # Create a mock window
        mock_window = object()
        window_provider = lambda: mock_window

        workflow = SilicaCAWorkflow(
            state_manager, window_provider=window_provider
        )

        params = SilicaCAFittingParams(
            dphi0=0.3,
            beamwaist=10e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        # Should still work fine with window_provider
        assert result.converged

    def test_fit_gracefully_handles_broken_window_provider(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that workflow gracefully handles broken window_provider."""

        # Create a broken window provider
        def broken_provider():
            raise RuntimeError("Window not available")

        workflow = SilicaCAWorkflow(
            state_manager, window_provider=broken_provider
        )

        params = SilicaCAFittingParams(
            dphi0=0.3,
            beamwaist=10e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        # Should gracefully fall back to headless mode (window=None)
        result = workflow.fit_manual(params)

        # Fitting should still succeed
        assert result.converged

    def test_multiple_datasets_independent(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that multiple fits don't interfere with each other."""
        workflow = SilicaCAWorkflow(state_manager, window_provider=None)

        # First fit
        params1 = SilicaCAFittingParams(
            dphi0=0.2,
            beamwaist=8e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result1 = workflow.fit_manual(params1)
        curve1 = result1.curve.copy()

        # Second fit with different parameters
        params2 = SilicaCAFittingParams(
            dphi0=0.4,
            beamwaist=14e-6,
            zero_level=0.95,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result2 = workflow.fit_manual(params2)
        curve2 = result2.curve.copy()

        # Curves should be different
        assert not np.allclose(curve1, curve2)