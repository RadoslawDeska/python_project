"""Comprehensive tests for sample fitting workflows."""

import numpy as np
import pytest
from numpy.typing import NDArray

from lib.state_manager import StateManager
from workflows.sample_workflow import (
    FittingResult,
    SampleCAFittingParams,
    SampleCAWorkflow,
    SampleOAFittingParams,
    SampleOAWorkflow,
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
    """Realistic sample CA data with peak."""
    z = sample_positions
    return 1.0 + 0.08 * np.exp(-(z**2) / 0.4)


@pytest.fixture
def sample_oa_data(sample_positions) -> NDArray:
    """Realistic sample OA data with absorption dip."""
    z = sample_positions
    return 0.92 - 0.04 * np.exp(-(z**2) / 0.8)


class TestSampleCAWorkflow:
    """Test sample closed aperture fitting."""

    def test_manual_fit_returns_fitting_result(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that manual CA fit returns valid FittingResult."""
        workflow = SampleCAWorkflow(state_manager, window_provider=None)

        params = SampleCAFittingParams(
            dphi0=0.3,
            beamwaist=50e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-21,  # Sample n2 is typically smaller than solvent
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
        assert result.message == "Sample CA manual fit successful"

    def test_manual_fit_curve_is_finite(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that fitted curve contains only finite values."""
        workflow = SampleCAWorkflow(state_manager, window_provider=None)

        params = SampleCAFittingParams(
            dphi0=0.25,
            beamwaist=45e-6,
            zero_level=0.98,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-21,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert result.converged
        assert np.all(np.isfinite(result.curve))

    def test_automatic_fit_returns_result_with_minimizer(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that automatic CA fit returns result with minimizer."""
        workflow = SampleCAWorkflow(state_manager, window_provider=None)

        params = SampleCAFittingParams(
            dphi0=0.2,
            beamwaist=40e-6,
            zero_level=0.96,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-21,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_automatic(params)

        assert isinstance(result, FittingResult)
        assert result.converged is True
        assert result.minimizer_result is not None
        assert result.message == "Sample CA automatic fit successful"

    def test_automatic_fit_includes_all_parameters(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that automatic fit returns all required parameters."""
        workflow = SampleCAWorkflow(state_manager, window_provider=None)

        params = SampleCAFittingParams(
            dphi0=0.28,
            beamwaist=48e-6,
            zero_level=0.97,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-21,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_automatic(params)

        assert result.converged
        assert "dphi0" in result.params
        assert "beamwaist" in result.params
        assert "center" in result.params
        assert "zero_level" in result.params

        # All should be finite numbers
        for key, value in result.params.items():
            assert isinstance(value, (int, float))
            assert np.isfinite(value)

    def test_automatic_fit_returns_errors(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that automatic fit includes error estimation."""
        workflow = SampleCAWorkflow(state_manager, window_provider=None)

        params = SampleCAFittingParams(
            dphi0=0.3,
            beamwaist=50e-6,
            zero_level=0.99,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-21,
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

    def test_invalid_dphi0_unrealistic_magnitude_fails(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that unrealistically large dphi0 magnitude is caught."""
        workflow = SampleCAWorkflow(state_manager, window_provider=None)

        params = SampleCAFittingParams(
            dphi0=500.0,  # Unrealistically large magnitude
            beamwaist=50e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-21,
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
        """Test that invalid beamwaist is caught."""
        workflow = SampleCAWorkflow(state_manager, window_provider=None)

        params = SampleCAFittingParams(
            dphi0=0.3,
            beamwaist=1.0,  # Way too large (1 meter)
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-21,
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
        workflow = SampleCAWorkflow(state_manager, window_provider=None)

        params = SampleCAFittingParams(
            dphi0=0.3,
            beamwaist=1e-9,  # 1 nanometer - unrealistic
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-21,
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
        workflow = SampleCAWorkflow(state_manager, window_provider=None)

        params = SampleCAFittingParams(
            dphi0=0.3,
            beamwaist=50e-6,
            zero_level=0.1,  # Too low
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-21,
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
        workflow = SampleCAWorkflow(state_manager, window_provider=None)

        params = SampleCAFittingParams(
            dphi0=0.3,
            beamwaist=50e-6,
            zero_level=10.0,  # Unrealistically high
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-21,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert result.converged is False

    def test_state_updated_on_manual_fit(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that state is updated correctly after manual fit."""
        assert state_manager.sampleCA_autofit_done is False

        workflow = SampleCAWorkflow(state_manager, window_provider=None)

        params = SampleCAFittingParams(
            dphi0=0.3,
            beamwaist=50e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-21,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert result.converged
        # Manual fit keeps flag False
        assert state_manager.sampleCA_autofit_done is False

    def test_state_updated_on_automatic_fit(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that state is updated on automatic fit."""
        assert state_manager.sampleCA_autofit_done is False

        workflow = SampleCAWorkflow(state_manager, window_provider=None)

        params = SampleCAFittingParams(
            dphi0=0.2,
            beamwaist=40e-6,
            zero_level=0.96,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-21,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_automatic(params)

        assert result.converged
        # Automatic fit sets flag to True
        assert state_manager.sampleCA_autofit_done is True

    def test_fit_with_window_provider(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that workflow works with window_provider (UI mode)."""
        mock_window = object()
        window_provider = lambda: mock_window

        workflow = SampleCAWorkflow(
            state_manager, window_provider=window_provider
        )

        params = SampleCAFittingParams(
            dphi0=0.3,
            beamwaist=50e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-21,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert result.converged

    def test_graceful_fallback_on_broken_window_provider(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test graceful fallback when window_provider fails."""

        def broken_provider():
            raise RuntimeError("Window unavailable")

        workflow = SampleCAWorkflow(
            state_manager, window_provider=broken_provider
        )

        params = SampleCAFittingParams(
            dphi0=0.3,
            beamwaist=50e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            ca_data=sample_ca_data,
            n2=1e-21,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        # Should still work by falling back to headless mode
        result = workflow.fit_manual(params)

        assert result.converged


class TestSampleOAWorkflow:
    """Test sample open aperture fitting."""

    def test_oa_fit_returns_result(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test that OA fit returns valid FittingResult."""
        workflow = SampleOAWorkflow(state_manager, window_provider=None)

        params = SampleOAFittingParams(
            transmittance=0.6,
            beamwaist=50e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            beta=1e-10,  # Two-photon absorption coefficient
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert isinstance(result, FittingResult)
        assert result.converged is True
        assert result.curve is not None
        assert len(result.curve) == len(sample_positions)
        assert result.message == "Sample OA manual fit successful"

    def test_oa_fit_curve_is_finite(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test that OA fitted curve contains only finite values."""
        workflow = SampleOAWorkflow(state_manager, window_provider=None)

        params = SampleOAFittingParams(
            transmittance=0.55,
            beamwaist=50e-6,
            zero_level=0.95,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            beta=1e-10,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert result.converged
        assert np.all(np.isfinite(result.curve))

    def test_oa_automatic_fit_optimizes(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test that OA automatic fit optimizes parameters."""
        workflow = SampleOAWorkflow(state_manager, window_provider=None)

        params = SampleOAFittingParams(
            transmittance=0.5,
            beamwaist=50e-6,
            zero_level=0.92,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            beta=1e-10,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_automatic(params)

        assert result.converged
        assert result.minimizer_result is not None
        assert "transmittance" in result.params
        assert result.message == "Sample OA automatic fit successful"

    def test_oa_automatic_fit_includes_errors(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test that OA automatic fit includes error estimation."""
        workflow = SampleOAWorkflow(state_manager, window_provider=None)

        params = SampleOAFittingParams(
            transmittance=0.58,
            beamwaist=50e-6,
            zero_level=0.93,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            beta=1e-10,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_automatic(params)

        assert result.converged
        assert "transmittance" in result.errors
        assert "beamwaist" in result.errors
        assert "center" in result.errors
        assert "zero_level" in result.errors

    def test_oa_invalid_transmittance_magnitude_too_large_fails(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test that unrealistically large transmittance magnitude is caught."""
        workflow = SampleOAWorkflow(state_manager, window_provider=None)

        params = SampleOAFittingParams(
            transmittance=50.0,  # Unrealistically large magnitude
            beamwaist=50e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            beta=1e-10,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert result.converged is False
        assert "transmittance" in result.message.lower()

    def test_oa_invalid_zero_level_too_low_fails(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test that invalid zero level is caught."""
        workflow = SampleOAWorkflow(state_manager, window_provider=None)

        params = SampleOAFittingParams(
            transmittance=0.6,
            beamwaist=50e-6,
            zero_level=0.05,  # Too low
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            beta=1e-10,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert result.converged is False

    def test_oa_invalid_beamwaist_too_small_fails(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test that unrealistically small beamwaist is caught."""
        workflow = SampleOAWorkflow(state_manager, window_provider=None)

        params = SampleOAFittingParams(
            transmittance=0.6,
            beamwaist=1e-9,  # 1 nanometer - unrealistic
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            beta=1e-10,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert result.converged is False

    def test_oa_state_updated_on_manual_fit(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test that state is updated after OA manual fit."""
        assert state_manager.sampleOA_autofit_done is False

        workflow = SampleOAWorkflow(state_manager, window_provider=None)

        params = SampleOAFittingParams(
            transmittance=0.6,
            beamwaist=50e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            beta=1e-10,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert result.converged
        # Manual fit keeps flag False
        assert state_manager.sampleOA_autofit_done is False

    def test_oa_state_updated_on_automatic_fit(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test that state is updated on OA automatic fit."""
        assert state_manager.sampleOA_autofit_done is False

        workflow = SampleOAWorkflow(state_manager, window_provider=None)

        params = SampleOAFittingParams(
            transmittance=0.5,
            beamwaist=50e-6,
            zero_level=0.92,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            beta=1e-10,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_automatic(params)

        assert result.converged
        # Automatic fit sets flag to True
        assert state_manager.sampleOA_autofit_done is True

    def test_oa_with_window_provider(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test OA workflow with window_provider."""
        mock_window = object()
        window_provider = lambda: mock_window

        workflow = SampleOAWorkflow(
            state_manager, window_provider=window_provider
        )

        params = SampleOAFittingParams(
            transmittance=0.6,
            beamwaist=50e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            beta=1e-10,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_manual(params)

        assert result.converged

    def test_oa_graceful_fallback_on_broken_provider(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test OA graceful fallback on broken window_provider."""

        def broken_provider():
            raise RuntimeError("Window unavailable")

        workflow = SampleOAWorkflow(
            state_manager, window_provider=broken_provider
        )

        params = SampleOAFittingParams(
            transmittance=0.6,
            beamwaist=50e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            beta=1e-10,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        # Should still work
        result = workflow.fit_manual(params)

        assert result.converged

    def test_oa_beamwaist_kept_fixed_in_automatic(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test that beamwaist is kept fixed during OA automatic fitting."""
        workflow = SampleOAWorkflow(state_manager, window_provider=None)

        initial_beamwaist = 50e-6

        params = SampleOAFittingParams(
            transmittance=0.5,
            beamwaist=initial_beamwaist,
            zero_level=0.92,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            beta=1e-10,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
        )

        result = workflow.fit_automatic(params)

        assert result.converged
        # Beamwaist should remain fixed (not optimized for OA)
        fitted_beamwaist = result.params["beamwaist"]
        assert fitted_beamwaist == initial_beamwaist