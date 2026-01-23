"""
Pure business logic for silica fitting workflows.
Designed to be testable without PyQt5 (headless).
UI integration happens at the window level.

NOTE: Silica is a reference material with known physics:
- In visible and NIR range, n2 > 0, so DPhi0 MUST be positive
- This uses silica-specific validation, not the generic validator
"""

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
from numpy.typing import NDArray

from lib.config import (
    INTEGRATION_STEPS,
    N_COMPONENTS,
    SILICA_BETA,
)
from lib.fitting import Fitting
from lib.integration import Integration
from lib.state_manager import StateManager
from lib.validation import ParameterValidator


@dataclass
class FittingResult:
    """Result from a fitting operation."""

    curve: NDArray
    params: dict
    errors: dict
    converged: bool
    message: str = ""
    minimizer_result: Optional[object] = None


@dataclass
class SilicaCAFittingParams:
    """Parameters for silica CA fitting."""

    dphi0: float
    beamwaist: float
    zero_level: float
    center_point: float
    z_positions: NDArray
    ca_data: NDArray
    n2: float
    d0: float
    ra: float
    wavelength: float
    z_range: float


class SilicaCAWorkflow:
    """Handles silica closed aperture fitting logic.

    This workflow is completely independent of the UI. It accepts an optional
    window_provider callback that can inject the window object if needed, but
    can also work without it for headless testing.
    
    NOTE: Uses silica-specific validation since DPhi0 must be positive in VIS/NIR.
    """

    def __init__(
        self, state: StateManager, window_provider: Optional[Callable] = None
    ):
        """
        Args:
            state: StateManager instance
            window_provider: Optional callable that returns the window object.
                           If None, Fitting.automatic() will receive None (headless mode).
                           This allows testing without PyQt5 while supporting UI integration.
        """
        self.state = state
        self.window_provider = window_provider

    def _get_window(self):
        """Get window object if provider is available, otherwise None."""
        if self.window_provider is None:
            return None
        try:
            return self.window_provider()
        except Exception:
            return None

    def fit_manual(self, params: SilicaCAFittingParams) -> FittingResult:
        """Manually fit CA data.

        Args:
            params: SilicaCAFittingParams with all required data

        Returns:
            FittingResult with curve and parameters
        """
        try:
            # Validate using SILICA-SPECIFIC validator
            # (DPhi0 must be positive for silica in VIS/NIR)
            valid, error = ParameterValidator.validate_silica_fitting_parameters(
                params.dphi0,
                params.beamwaist,
                params.zero_level,
                params.center_point,
            )
            if not valid:
                return FittingResult(
                    curve=np.array([]),
                    params={},
                    errors={},
                    converged=False,
                    message=f"CA validation failed: {error}",
                )

            # Create integration object
            integration = Integration(
                beta=SILICA_BETA,
                n2=params.n2,
                DPhi0=params.dphi0,
                positions=params.z_positions,
                d0=params.d0,
                aperture_radius=params.ra,
                wavelength=params.wavelength,
                beamwaist=params.beamwaist,
                n_components=N_COMPONENTS,
                integration_steps=INTEGRATION_STEPS,
                stype="CA",
            )

            # Create fitter and calculate curve
            fitter = Fitting(
                integration,
                params.dphi0,
                params.beamwaist,
                params.zero_level,
                params.center_point,
                len(params.z_positions),
                params.ca_data,
            )

            result_curve = fitter.manual(
                params.zero_level,
                params.center_point,
                params.dphi0,
                params.beamwaist,
                params.z_range,
                window=self._get_window(),
                stype="CA",
            )

            # Update state
            self.state.silica_autofit_done = False

            return FittingResult(
                curve=result_curve,
                params={"dphi0": params.dphi0, "beamwaist": params.beamwaist},
                errors={},
                converged=True,
                message="Manual CA fit successful",
            )

        except Exception as e:
            return FittingResult(
                curve=np.array([]),
                params={},
                errors={},
                converged=False,
                message=f"CA fitting failed: {str(e)}",
            )

    def fit_automatic(self, params: SilicaCAFittingParams) -> FittingResult:
        """Automatically fit CA data using optimization.

        Args:
            params: SilicaCAFittingParams with all required data

        Returns:
            FittingResult with optimized curve and parameters
        """
        try:
            # Validate using SILICA-SPECIFIC validator
            valid, error = ParameterValidator.validate_silica_fitting_parameters(
                params.dphi0,
                params.beamwaist,
                params.zero_level,
                params.center_point,
            )
            if not valid:
                return FittingResult(
                    curve=np.array([]),
                    params={},
                    errors={},
                    converged=False,
                    message=f"Validation failed: {error}",
                )

            # Create integration
            integration = Integration(
                beta=SILICA_BETA,
                n2=params.n2,
                DPhi0=params.dphi0,
                positions=params.z_positions,
                d0=params.d0,
                aperture_radius=params.ra,
                wavelength=params.wavelength,
                beamwaist=params.beamwaist,
                n_components=N_COMPONENTS,
                integration_steps=INTEGRATION_STEPS,
                stype="CA",
            )

            # Create fitter
            fitter = Fitting(
                integration,
                params.dphi0,
                params.beamwaist,
                params.zero_level,
                params.center_point,
                len(params.z_positions),
                params.ca_data,
            )

            # Perform automatic fit - may receive window=None in headless mode
            minimizer_result, result_curve = fitter.automatic(
                z_range=params.z_range,
                ftype="Silica",
                stype="CA",
                line_xydata=(params.z_positions, params.ca_data),
                window=self._get_window(),
            )

            # Extract results
            fitted_params = {
                "dphi0": minimizer_result.params["DPhi0"].value,
                "beamwaist": minimizer_result.params["Beamwaist"].value,
                "center": minimizer_result.params["Center"].value,
                "zero_level": minimizer_result.params["Zero"].value,
            }

            fitted_errors = {
                "dphi0": minimizer_result.params["DPhi0"].stderr,
                "beamwaist": minimizer_result.params["Beamwaist"].stderr,
                "center": minimizer_result.params["Center"].stderr,
                "zero_level": minimizer_result.params["Zero"].stderr,
            }

            # Update state
            self.state.silica_autofit_done = True

            return FittingResult(
                curve=result_curve,
                params=fitted_params,
                errors=fitted_errors,
                converged=True,
                message="Automatic CA fit successful",
                minimizer_result=minimizer_result,
            )

        except Exception as e:
            return FittingResult(
                curve=np.array([]),
                params={},
                errors={},
                converged=False,
                message=f"CA automatic fitting failed: {str(e)}",
            )