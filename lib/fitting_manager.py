"""Fitting Manager - Handles all fitting operations.

Manages curve fitting workflows for CA and OA measurements, including automatic
and manual fitting modes.
"""

from typing import Dict, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from lib.fitting import Fitting
from lib.integration import Integration
from lib.sample_fitting import SampleFitter


class FittingManager:
    """Manages fitting operations for Z-scan measurements."""

    def __init__(self) -> None:
        self.integration: Optional[Integration] = None
        self.fitting: Optional[Fitting] = None
        self.sample_fitter = SampleFitter()
        self.last_fit_result: Optional[dict] = None
        self.last_fit_curve: Optional[np.ndarray] = None

    def create_fitter(
        self,
        amplitude: float,
        beamwaist: float,
        zero_level: float,
        centerpoint: float,
        y_data: NDArray,
    ) -> Fitting:
        """Create a Fitting object for curve fitting.

        Assumes `self.integration` has already been created.
        """
        if self.integration is None:
            raise RuntimeError(
                "Integration must be created before creating fitter"
            )

        nop = int(len(y_data))

        self.fitting = Fitting(
            sample_type=self.integration,
            amplitude=amplitude,
            beamwaist=beamwaist,
            zero_level=zero_level,
            centerpoint=centerpoint,
            nop=nop,
            data=y_data,
        )
        return self.fitting

    def create_integration(
        self,
        beta: float,
        n2: float,
        dphi0: float,
        wavelength: float,
        beamwaist: float,
        positions: NDArray,
        d0: float,
        aperture_radius: float,
        n_components: int,
        integration_steps: int,
        stype: str = "CA",
    ) -> Integration:
        """Create an Integration object for physics calculations.

        Args:
            beta: Nonlinear absorption coefficient
            n2: Nonlinear refractive index
            dphi0: Phase shift
            amplitude: Field amplitude
            waist: Beam waist
            wavelength: Wavelength
            z_points: Z-position array

        Returns:
            Configured Integration object
        """
        self.integration = Integration(
            beta=beta,
            n2=n2,
            DPhi0=dphi0,
            positions=positions,
            d0=d0,
            aperture_radius=aperture_radius,
            wavelength=wavelength,
            beamwaist=beamwaist,
            n_components=n_components,
            integration_steps=integration_steps,
            stype=stype,
        )
        assert self.integration is not None
        return self.integration

    def fit_ca_automatic(
        self,
        z_range: float,
        ftype: str,
        x_data: NDArray,
        y_data: NDArray,
        window: object,
        vary_beamwaist: bool = True,
        vary_centerpoint: bool = True,
    ) -> Tuple[Dict, NDArray]:
        """Perform automatic CA (Closed Aperture) fitting.

        Args:
            x_data: X-axis data (positions)
            y_data: Y-axis data (normalized transmittance)
            initial_amplitude: Initial amplitude estimate
            initial_offset: Initial offset estimate
            initial_dphi0: Initial phase shift estimate
            roi_start: Region of interest start index
            roi_end: Region of interest end index

        Returns:
            Tuple of (fit_result dict, fitted curve array)
        """
        if self.fitting is None or self.integration is None:
            raise RuntimeError("Must create integration and fitter first")

        result, curve = self.fitting.automatic(
            z_range=z_range,
            ftype=ftype,
            stype="CA",
            line_xydata=(x_data, y_data),
            window=window,
            vary_beamwaist=vary_beamwaist,
            vary_centerpoint=vary_centerpoint,
        )

        self.last_fit_result = result
        self.last_fit_curve = curve
        return result, curve

    def fit_oa_automatic(
        self,
        z_range: float,
        ftype: str,
        x_data: NDArray,
        y_data: NDArray,
        window: object,
        vary_centerpoint: bool = True,
    ) -> Tuple[Dict, NDArray]:
        """Perform automatic OA (Open Aperture) fitting.

        Args:
            x_data: X-axis data
            y_data: Y-axis data
            initial_beta: Initial nonlinear absorption estimate
            initial_offset: Initial offset estimate
            absorption_model: Absorption model type
            roi_start: Region of interest start index
            roi_end: Region of interest end index

        Returns:
            Tuple of (fit_result dict, fitted curve array)
        """
        if self.fitting is None or self.integration is None:
            raise RuntimeError("Must create integration and fitter first")

        result, curve = self.fitting.automatic(
            z_range=z_range,
            ftype=ftype,
            stype="OA",
            line_xydata=(x_data, y_data),
            window=window,
            vary_beamwaist=False,  # OA code fixes beamwaist internally
            vary_centerpoint=vary_centerpoint,
        )

        self.last_fit_result = result
        self.last_fit_curve = curve
        return result, curve

    def fit_manual(
        self,
        zero_level: float,
        centerpoint: float,
        amplitude: float,
        beamwaist: float,
        z_range: float,
        window: object,
        stype: str,
    ) -> NDArray:
        """Calculate manual fit without optimization.

        Args:
            zero_level: Baseline transmittance level
            centerpoint: Center position of curve [data points]
            amplitude: Phase shift (CA) or transmittance (OA)
            beamwaist: Beam waist radius [m]
            z_range: Total z-scan range [m]
            window: Main window object (for hardware parameters)
            stype: "CA" or "OA"

        Returns:
            Fitted curve array
        """
        if self.fitting is None:
            raise RuntimeError("Fitter must be created before manual fit")

        curve = self.fitting.manual(
            zero_level=zero_level,
            centerpoint=centerpoint,
            amplitude=amplitude,
            beamwaist=beamwaist,
            z_range=z_range,
            window=window,
            stype=stype,
        )
        return curve

    def extract_ca_parameters(self, fit_result: Dict) -> Tuple[Dict, Dict]:
        """Extract CA parameters and errors from fit result.

        Args:
            fit_result: Result from automatic fitting

        Returns:
            Tuple of (parameters dict, errors dict)
        """
        params = {}
        errors = {}

        if "params" in fit_result:
            params = fit_result["params"]
        if "errors" in fit_result:
            errors = fit_result["errors"]

        return params, errors

    def extract_oa_parameters(self, fit_result: Dict) -> Tuple[Dict, Dict]:
        """Extract OA parameters and errors from fit result.

        Args:
            fit_result: Result from automatic fitting

        Returns:
            Tuple of (parameters dict, errors dict)
        """
        params = {}
        errors = {}

        if "params" in fit_result:
            params = fit_result["params"]
        if "errors" in fit_result:
            errors = fit_result["errors"]

        return params, errors

    def validate_sample_prerequisites(
        self, silica_fitted: bool, solvent_fitted: bool
    ) -> Tuple[bool, str]:
        """Validate that sample fitting prerequisites are met.

        Args:
            silica_fitted: Whether silica reference has been fitted
            solvent_fitted: Whether solvent baseline has been fitted

        Returns:
            Tuple of (is_valid, message)
        """
        if not silica_fitted:
            return False, "Silica reference must be fitted first"
        if not solvent_fitted:
            return False, "Solvent baseline must be fitted first"
        return True, "All prerequisites met"

    def get_last_fit_curve(self) -> Optional[np.ndarray]:
        """Get the most recently calculated fit curve.

        Returns:
            Last fit curve array or None
        """
        return self.last_fit_curve

    def get_last_fit_result(self) -> Optional[Dict]:
        """Get the most recent fit result.

        Returns:
            Last fit result dict or None
        """
        return self.last_fit_result

    def reset(self):
        """Reset fitting manager state."""
        self.integration = None
        self.fitting = None
        self.last_fit_result = None
        self.last_fit_curve = None
