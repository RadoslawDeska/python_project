"""
Pure business logic for sample fitting workflows.
No UI dependencies - can be tested without PyQt5.
"""

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
from numpy.typing import NDArray

from lib.config import (
    INTEGRATION_STEPS,
    N_COMPONENTS,
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
class SampleCAFittingParams:
    """Parameters for sample CA fitting."""

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


@dataclass
class SampleOAFittingParams:
    """Parameters for sample OA fitting."""

    transmittance: float
    beamwaist: float
    zero_level: float
    center_point: float
    z_positions: NDArray
    oa_data: NDArray
    beta: float
    d0: float
    ra: float
    wavelength: float
    z_range: float


class SampleCAWorkflow:
    """Handles sample closed aperture fitting logic.

    This workflow is completely independent of the UI. It accepts an optional
    window_provider callback that can inject the window object if needed, but
    can also work without it for headless testing.
    """

    def __init__(
        self, state: StateManager, window_provider: Optional[Callable] = None
    ):
        """
        Args:
            state: StateManager instance
            window_provider: Optional callable that returns the window object.
                           If None, workflows operate in headless mode.
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

    def fit_manual(self, params: SampleCAFittingParams) -> FittingResult:
        """Manually fit sample CA data."""
        try:
            # Validate individual parameters - numerical validity and bounds only
            valid, error = ParameterValidator.validate_beamwaist(params.beamwaist)
            if not valid:
                return FittingResult(
                    np.array([]), {}, {}, False, f"Invalid beamwaist: {error}"
                )
            
            valid, error = ParameterValidator.validate_zero_level(params.zero_level)
            if not valid:
                return FittingResult(
                    np.array([]), {}, {}, False, f"Invalid zero_level: {error}"
                )
            
            # DPhi0 can be positive or negative (material-dependent)
            # Just check it's numerically reasonable
            if not isinstance(params.dphi0, (int, float)):
                return FittingResult(
                    np.array([]), {}, {}, False, "DPhi0 must be a number"
                )
            if abs(params.dphi0) > 100:
                return FittingResult(
                    np.array([]), {}, {}, False, 
                    f"DPhi0 magnitude unrealistic: {abs(params.dphi0)}"
                )
            
            # Use the general validator for overall consistency
            valid, error = ParameterValidator.validate_fitting_parameters(
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
                    message=f"Sample CA validation failed: {error}",
                )

            # Create integration object for CA
            integration = Integration(
                beta=0,
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

            # Perform manual fit
            result_curve = fitter.manual(
                zero_level=params.zero_level,
                centerpoint=params.center_point,
                amplitude=params.dphi0,
                beamwaist=params.beamwaist,
                z_range=params.z_range,
                window=self._get_window(),
                stype="CA",
            )

            self.state.sampleCA_autofit_done = False

            return FittingResult(
                curve=result_curve,
                params={
                    "dphi0": params.dphi0,
                    "beamwaist": params.beamwaist,
                    "center": params.center_point,
                    "zero_level": params.zero_level,
                },
                errors={},
                converged=True,
                message="Sample CA manual fit successful",
            )

        except Exception as e:
            return FittingResult(
                curve=np.array([]),
                params={},
                errors={},
                converged=False,
                message=f"Sample CA manual fit failed: {str(e)}",
            )

    def fit_automatic(self, params: SampleCAFittingParams) -> FittingResult:
        """Automatically fit sample CA data."""
        try:
            # Validate
            valid, error = ParameterValidator.validate_fitting_parameters(
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

            # Create integration object for CA
            integration = Integration(
                beta=0,
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

            # Perform automatic fit
            minimizer_result, result_curve = fitter.automatic(
                z_range=params.z_range,
                ftype="Sample",
                stype="CA",
                line_xydata=(params.z_positions, params.ca_data),
                window=self._get_window(),
                vary_beamwaist=True,
                vary_centerpoint=True,
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

            self.state.sampleCA_autofit_done = True

            return FittingResult(
                curve=result_curve,
                params=fitted_params,
                errors=fitted_errors,
                converged=True,
                message="Sample CA automatic fit successful",
                minimizer_result=minimizer_result,
            )

        except Exception as e:
            return FittingResult(
                curve=np.array([]),
                params={},
                errors={},
                converged=False,
                message=f"Sample CA automatic fitting failed: {str(e)}",
            )


class SampleOAWorkflow:
    """Handles sample open aperture fitting logic.

    This workflow is completely independent of the UI. It accepts an optional
    window_provider callback that can inject the window object if needed, but
    can also work without it for headless testing.
    """

    def __init__(
        self, state: StateManager, window_provider: Optional[Callable] = None
    ):
        """
        Args:
            state: StateManager instance
            window_provider: Optional callable that returns the window object.
                           If None, workflows operate in headless mode.
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

    def fit_manual(self, params: SampleOAFittingParams) -> FittingResult:
        """Manually fit sample OA data."""
        try:
            # Validate transmittance parameter
            # T can be positive or negative depending on model - just check bounds
            valid, error = ParameterValidator.validate_oa_parameters(params.transmittance)
            if not valid:
                return FittingResult(
                    np.array([]), {}, {}, False, f"Invalid transmittance: {error}"
                )
            
            # Validate zero level
            valid, error = ParameterValidator.validate_zero_level(params.zero_level)
            if not valid:
                return FittingResult(
                    np.array([]), {}, {}, False, f"Invalid zero_level: {error}"
                )
            
            # Validate beamwaist
            valid, error = ParameterValidator.validate_beamwaist(params.beamwaist)
            if not valid:
                return FittingResult(
                    np.array([]), {}, {}, False, f"Invalid beamwaist: {error}"
                )

            # Create integration object for OA
            integration = Integration(
                beta=params.beta,
                n2=1e-15,  # n2 not used in OA fitting
                DPhi0=0,
                positions=params.z_positions,
                d0=params.d0,
                aperture_radius=params.ra,
                wavelength=params.wavelength,
                beamwaist=params.beamwaist,
                n_components=N_COMPONENTS,
                integration_steps=INTEGRATION_STEPS,
                stype="OA",
            )

            # Create fitter
            fitter = Fitting(
                integration,
                params.transmittance,
                params.beamwaist,
                params.zero_level,
                params.center_point,
                len(params.z_positions),
                params.oa_data,
            )

            # Perform manual fit
            result_curve = fitter.manual(
                zero_level=params.zero_level,
                centerpoint=params.center_point,
                amplitude=params.transmittance,
                beamwaist=params.beamwaist,
                z_range=params.z_range,
                window=self._get_window(),
                stype="OA",
            )

            self.state.sampleOA_autofit_done = False

            return FittingResult(
                curve=result_curve,
                params={
                    "transmittance": params.transmittance,
                    "beamwaist": params.beamwaist,
                    "center": params.center_point,
                    "zero_level": params.zero_level,
                },
                errors={},
                converged=True,
                message="Sample OA manual fit successful",
            )

        except Exception as e:
            return FittingResult(
                curve=np.array([]),
                params={},
                errors={},
                converged=False,
                message=f"Sample OA manual fit failed: {str(e)}",
            )

    def fit_automatic(self, params: SampleOAFittingParams) -> FittingResult:
        """Automatically fit sample OA data."""
        try:
            # Create integration object for OA
            integration = Integration(
                beta=params.beta,
                n2=1e-15,  # n2 not used in OA fitting
                DPhi0=0,
                positions=params.z_positions,
                d0=params.d0,
                aperture_radius=params.ra,
                wavelength=params.wavelength,
                beamwaist=params.beamwaist,
                n_components=N_COMPONENTS,
                integration_steps=INTEGRATION_STEPS,
                stype="OA",
            )

            # Create fitter
            fitter = Fitting(
                integration,
                params.transmittance,
                params.beamwaist,
                params.zero_level,
                params.center_point,
                len(params.z_positions),
                params.oa_data,
            )

            # Perform automatic fit
            minimizer_result, result_curve = fitter.automatic(
                z_range=params.z_range,
                ftype="Sample",
                stype="OA",
                line_xydata=(params.z_positions, params.oa_data),
                window=self._get_window(),
                vary_beamwaist=False,  # Keep beamwaist fixed for OA
                vary_centerpoint=True,
            )

            # Extract results
            fitted_params = {
                "transmittance": minimizer_result.params.get(
                    "T",
                    minimizer_result.params.get(
                        "Amplitude", params.transmittance
                    ),
                ),
                "beamwaist": minimizer_result.params["Beamwaist"].value,
                "center": minimizer_result.params["Center"].value,
                "zero_level": minimizer_result.params["Zero"].value,
            }

            fitted_errors = {
                "transmittance": minimizer_result.params.get(
                    "T", minimizer_result.params.get("Amplitude")
                ).stderr,
                "beamwaist": minimizer_result.params["Beamwaist"].stderr,
                "center": minimizer_result.params["Center"].stderr,
                "zero_level": minimizer_result.params["Zero"].stderr,
            }

            self.state.sampleOA_autofit_done = True

            return FittingResult(
                curve=result_curve,
                params=fitted_params,
                errors=fitted_errors,
                converged=True,
                message="Sample OA automatic fit successful",
                minimizer_result=minimizer_result,
            )

        except Exception as e:
            return FittingResult(
                curve=np.array([]),
                params={},
                errors={},
                converged=False,
                message=f"Sample OA automatic fitting failed: {str(e)}",
            )
