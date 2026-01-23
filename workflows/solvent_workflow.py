"""
Pure business logic for solvent fitting workflows.
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
class SolventCAFittingParams:
    """Parameters for solvent CA fitting."""

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
class SolventOAFittingParams:
    """Parameters for solvent OA fitting."""

    t_value: float  # Transmittance/amplitude
    beamwaist: float
    zero_level: float
    center_point: float
    z_positions: NDArray
    oa_data: NDArray
    n2: float
    d0: float
    ra: float
    wavelength: float
    z_range: float
    no_absorption: bool = False


class SolventCAWorkflow:
    """Handles solvent closed aperture fitting logic.

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

    def fit_manual(self, params: SolventCAFittingParams) -> FittingResult:
        """Manually fit CA data.

        Args:
            params: SolventCAFittingParams with all required data

        Returns:
            FittingResult with curve and parameters
        """
        try:
            # Validate if parameters meet the fitting bounds from configuration
            from lib.constants import CA_FITTING_PARAMS

            if not (
                CA_FITTING_PARAMS["Zero"]["min"]
                <= params.zero_level
                <= CA_FITTING_PARAMS["Zero"]["max"]
            ):
                return FittingResult(
                    np.array([]), {}, {}, False, "Invalid zero_level"
                )
            if not (
                CA_FITTING_PARAMS["DPhi0"]["min"]
                <= params.dphi0
                <= CA_FITTING_PARAMS["DPhi0"]["max"]
            ):
                return FittingResult(
                    np.array([]), {}, {}, False, "Invalid dphi0"
                )

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
                    message=f"CA validation failed: {error}",
                )

            # Create integration object
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
            self.state.solventCA_autofit_done = False

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

    def fit_automatic(self, params: SolventCAFittingParams) -> FittingResult:
        """Automatically fit CA data using optimization.

        Args:
            params: SolventCAFittingParams with all required data

        Returns:
            FittingResult with optimized curve and parameters
        """
        try:
            # Validate inputs
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

            # Create integration
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
                ftype="Solvent",
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

            # Update state
            self.state.solventCA_autofit_done = True

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


class SolventOAWorkflow:
    """Handles solvent open aperture fitting logic.

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

    def fit_manual(self, params: SolventOAFittingParams) -> Optional[FittingResult]:
        """Manually fit OA data.

        Returns:
            FittingResult or None if fitting should be skipped
        """
        # Skip if no absorption assumed
        if params.no_absorption:
            return None

        try:
            # Validate T parameter - can be positive or negative
            valid, error = ParameterValidator.validate_oa_parameters(params.t_value)
            if not valid:
                return FittingResult(
                    np.array([]), {}, {}, False, f"Invalid t_value: {error}"
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

            # Create integration for OA
            integration = Integration(
                beta=params.t_value,
                n2=params.n2,
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

            # Fit
            fitter = Fitting(
                integration,
                params.t_value,
                params.beamwaist,
                params.zero_level,
                params.center_point,
                len(params.z_positions),
                params.oa_data,
            )

            result_curve = fitter.manual(
                params.zero_level,
                params.center_point,
                params.t_value,
                params.beamwaist,
                params.z_range,
                window=self._get_window(),
                stype="OA",
            )

            return FittingResult(
                curve=result_curve,
                params={"t": params.t_value},
                errors={},
                converged=True,
                message="OA manual fit successful",
            )

        except Exception as e:
            return FittingResult(
                curve=np.array([]),
                params={},
                errors={},
                converged=False,
                message=f"OA manual fit failed: {str(e)}",
            )

    def fit_automatic(
        self, params: SolventOAFittingParams
    ) -> Optional[FittingResult]:
        """Automatically fit OA data using optimization.

        Returns:
            FittingResult or None if fitting should be skipped
        """
        # Skip if no absorption assumed
        if params.no_absorption:
            return None

        try:
            # Validate T parameter
            valid, error = ParameterValidator.validate_oa_parameters(
                params.t_value
            )
            if not valid:
                return FittingResult(
                    curve=np.array([]),
                    params={},
                    errors={},
                    converged=False,
                    message=f"OA validation failed: {error}",
                )

            # Create integration for OA
            integration = Integration(
                beta=params.t_value,
                n2=params.n2,
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
                params.t_value,
                params.beamwaist,
                params.zero_level,
                params.center_point,
                len(params.z_positions),
                params.oa_data,
            )

            # Perform automatic fit
            minimizer_result, result_curve = fitter.automatic(
                z_range=params.z_range,
                ftype="Solvent",
                stype="OA",
                line_xydata=(params.z_positions, params.oa_data),
                window=self._get_window(),
                vary_beamwaist=False,  # Keep beamwaist fixed for OA
                vary_centerpoint=True,
            )

            # Extract results
            fitted_params = {
                "t": minimizer_result.params.get(
                    "T",
                    minimizer_result.params.get("Amplitude", params.t_value),
                ),
                "beamwaist": minimizer_result.params["Beamwaist"].value,
                "center": minimizer_result.params["Center"].value,
                "zero_level": minimizer_result.params["Zero"].value,
            }

            fitted_errors = {
                "t": minimizer_result.params.get(
                    "T", minimizer_result.params.get("Amplitude")
                ).stderr,
                "beamwaist": minimizer_result.params["Beamwaist"].stderr,
                "center": minimizer_result.params["Center"].stderr,
                "zero_level": minimizer_result.params["Zero"].stderr,
            }

            return FittingResult(
                curve=result_curve,
                params=fitted_params,
                errors=fitted_errors,
                converged=True,
                message="OA automatic fit successful",
                minimizer_result=minimizer_result,
            )

        except Exception as e:
            return FittingResult(
                curve=np.array([]),
                params={},
                errors={},
                converged=False,
                message=f"OA automatic fitting failed: {str(e)}",
            )
