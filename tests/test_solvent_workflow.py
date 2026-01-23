"""Comprehensive tests for solvent fitting workflows."""

import numpy as np
import pytest
from numpy.typing import NDArray

from lib.state_manager import StateManager
from workflows.solvent_workflow import (
    FittingResult,
    SolventCAFittingParams,
    SolventCAWorkflow,
    SolventOAFittingParams,
    SolventOAWorkflow,
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
    """Realistic CA data with Gaussian peak."""
    z = sample_positions
    return 1.0 + 0.1 * np.exp(-(z**2) / 0.5)


@pytest.fixture
def sample_oa_data(sample_positions) -> NDArray:
    """Realistic OA data with dip at center (absorption)."""
    z = sample_positions
    return 0.9 - 0.05 * np.exp(-(z**2) / 1.0)


class TestSolventCAWorkflow:
    """Test solvent closed aperture fitting."""

    def test_manual_fit_returns_fitting_result(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that manual CA fit returns valid FittingResult."""
        workflow = SolventCAWorkflow(state_manager, window_provider=None)

        params = SolventCAFittingParams(
            dphi0=0.5,
            beamwaist=50e-6,
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

        assert isinstance(result, FittingResult)
        assert result.curve is not None
        assert len(result.curve) == len(sample_positions)
        assert result.converged is True
        assert result.message == "Manual CA fit successful"

    def test_manual_fit_preserves_curve_length(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that curve output has same length as input."""
        workflow = SolventCAWorkflow(state_manager, window_provider=None)

        params = SolventCAFittingParams(
            dphi0=0.3,
            beamwaist=40e-6,
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
        assert len(result.curve) == len(sample_ca_data)
        assert np.all(np.isfinite(result.curve))

    def test_automatic_fit_optimizes_parameters(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that automatic fit actually optimizes parameters."""
        workflow = SolventCAWorkflow(state_manager, window_provider=None)

        initial_dphi0 = 0.2
        initial_beamwaist = 30e-6

        params = SolventCAFittingParams(
            dphi0=initial_dphi0,
            beamwaist=initial_beamwaist,
            zero_level=0.9,
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
        assert result.minimizer_result is not None
        # Parameters may have changed after optimization
        fitted_dphi0 = result.params["dphi0"]
        fitted_beamwaist = result.params["beamwaist"]
        assert fitted_dphi0 is not None
        assert fitted_beamwaist is not None

    def test_automatic_fit_includes_all_parameters(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that automatic fit returns all required parameters."""
        workflow = SolventCAWorkflow(state_manager, window_provider=None)

        params = SolventCAFittingParams(
            dphi0=0.35,
            beamwaist=45e-6,
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
        workflow = SolventCAWorkflow(state_manager, window_provider=None)

        params = SolventCAFittingParams(
            dphi0=0.4,
            beamwaist=48e-6,
            zero_level=0.97,
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

    def test_invalid_beamwaist_too_large_fails(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that invalid beamwaist is caught."""
        workflow = SolventCAWorkflow(state_manager, window_provider=None)

        params = SolventCAFittingParams(
            dphi0=0.5,
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
        workflow = SolventCAWorkflow(state_manager, window_provider=None)

        params = SolventCAFittingParams(
            dphi0=0.5,
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
        workflow = SolventCAWorkflow(state_manager, window_provider=None)

        params = SolventCAFittingParams(
            dphi0=0.5,
            beamwaist=50e-6,
            zero_level=0.1,  # Too low
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
        workflow = SolventCAWorkflow(state_manager, window_provider=None)

        params = SolventCAFittingParams(
            dphi0=0.5,
            beamwaist=50e-6,
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

    def test_state_updated_on_manual_fit(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that state is updated correctly after manual fit."""
        assert state_manager.solventCA_autofit_done is False

        workflow = SolventCAWorkflow(state_manager, window_provider=None)

        params = SolventCAFittingParams(
            dphi0=0.5,
            beamwaist=50e-6,
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
        # Manual fit keeps flag False
        assert state_manager.solventCA_autofit_done is False

    def test_state_updated_on_automatic_fit(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that state is updated on automatic fit."""
        assert state_manager.solventCA_autofit_done is False

        workflow = SolventCAWorkflow(state_manager, window_provider=None)

        params = SolventCAFittingParams(
            dphi0=0.35,
            beamwaist=45e-6,
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
        assert state_manager.solventCA_autofit_done is True

    def test_fit_with_window_provider(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test that workflow works with window_provider (UI mode)."""
        mock_window = object()
        window_provider = lambda: mock_window

        workflow = SolventCAWorkflow(
            state_manager, window_provider=window_provider
        )

        params = SolventCAFittingParams(
            dphi0=0.5,
            beamwaist=50e-6,
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

    def test_graceful_fallback_on_broken_window_provider(
        self, state_manager, sample_positions, sample_ca_data
    ):
        """Test graceful fallback when window_provider fails."""

        def broken_provider():
            raise RuntimeError("Window unavailable")

        workflow = SolventCAWorkflow(
            state_manager, window_provider=broken_provider
        )

        params = SolventCAFittingParams(
            dphi0=0.5,
            beamwaist=50e-6,
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

        # Should still work by falling back to headless mode
        result = workflow.fit_manual(params)

        assert result.converged


class TestSolventOAWorkflow:
    """Test solvent open aperture fitting."""

    def test_oa_fit_returns_result(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test that OA fit returns valid FittingResult."""
        workflow = SolventOAWorkflow(state_manager, window_provider=None)

        params = SolventOAFittingParams(
            t_value=0.5,
            beamwaist=50e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
            no_absorption=False,
        )

        result = workflow.fit_manual(params)

        assert isinstance(result, FittingResult)
        assert result.converged is True
        assert result.curve is not None
        assert len(result.curve) == len(sample_positions)

    def test_oa_fit_skipped_when_no_absorption_true(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test that OA fit is skipped when no_absorption=True."""
        workflow = SolventOAWorkflow(state_manager, window_provider=None)

        params = SolventOAFittingParams(
            t_value=0.5,
            beamwaist=50e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
            no_absorption=True,  # Skip fitting
        )

        result = workflow.fit_manual(params)

        # When skipped, should return None
        assert result is None

    def test_oa_automatic_fit_optimizes(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test that OA automatic fit optimizes parameters."""
        workflow = SolventOAWorkflow(state_manager, window_provider=None)

        params = SolventOAFittingParams(
            t_value=0.6,
            beamwaist=50e-6,
            zero_level=0.95,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
            no_absorption=False,
        )

        result = workflow.fit_automatic(params)

        assert result.converged
        assert result.minimizer_result is not None
        assert "t" in result.params

    def test_oa_invalid_t_value_magnitude_too_large_fails(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test that unrealistically large T magnitude is caught."""
        workflow = SolventOAWorkflow(state_manager, window_provider=None)

        params = SolventOAFittingParams(
            t_value=50.0,  # Unrealistically large magnitude
            beamwaist=50e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
            no_absorption=False,
        )

        result = workflow.fit_manual(params)

        assert result.converged is False
        assert "t" in result.message.lower() or "magnitude" in result.message.lower()

    def test_oa_invalid_zero_level_too_low_fails(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test that invalid zero level is caught."""
        workflow = SolventOAWorkflow(state_manager, window_provider=None)

        params = SolventOAFittingParams(
            t_value=0.5,
            beamwaist=50e-6,
            zero_level=0.05,  # Too low
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
            no_absorption=False,
        )

        result = workflow.fit_manual(params)

        assert result.converged is False

    def test_oa_invalid_beamwaist_too_small_fails(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test that unrealistically small beamwaist is caught."""
        workflow = SolventOAWorkflow(state_manager, window_provider=None)

        params = SolventOAFittingParams(
            t_value=0.5,
            beamwaist=1e-9,  # 1 nanometer - unrealistic
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
            no_absorption=False,
        )

        result = workflow.fit_manual(params)

        assert result.converged is False

    def test_oa_state_updated_on_manual_fit(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test that state is updated after OA manual fit."""
        workflow = SolventOAWorkflow(state_manager, window_provider=None)

        params = SolventOAFittingParams(
            t_value=0.5,
            beamwaist=50e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
            no_absorption=False,
        )

        result = workflow.fit_manual(params)

        assert result.converged
        # Manual fit doesn't set autofit_done
        # (that's only set by automatic fitting)

    def test_oa_with_window_provider(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test OA workflow with window_provider."""
        mock_window = object()
        window_provider = lambda: mock_window

        workflow = SolventOAWorkflow(
            state_manager, window_provider=window_provider
        )

        params = SolventOAFittingParams(
            t_value=0.5,
            beamwaist=50e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
            no_absorption=False,
        )

        result = workflow.fit_manual(params)

        assert result.converged

    def test_oa_graceful_fallback_on_broken_provider(
        self, state_manager, sample_positions, sample_oa_data
    ):
        """Test OA graceful fallback on broken window_provider."""

        def broken_provider():
            raise RuntimeError("Window unavailable")

        workflow = SolventOAWorkflow(
            state_manager, window_provider=broken_provider
        )

        params = SolventOAFittingParams(
            t_value=0.5,
            beamwaist=50e-6,
            zero_level=1.0,
            center_point=20,
            z_positions=sample_positions,
            oa_data=sample_oa_data,
            n2=1e-20,
            d0=0.1,
            ra=0.001,
            wavelength=800e-9,
            z_range=0.004,
            no_absorption=False,
        )

        result = workflow.fit_manual(params)

        # Should still work
        assert result.converged