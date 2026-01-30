"""
Fitting module for Z-Scan curve fitting.

This module contains the Fitting class for both manual and automatic
curve fitting of Z-scan data using parametric functions and optimization.

Classes:
    Fitting: Performs manual and automatic curve fitting for Z-scan measurements.
"""

from typing import Any, List, Optional, Tuple
import logging

import numpy as np
from lmfit import Minimizer, Parameters
from numpy.typing import NDArray

from lib.constants import CA_FITTING_PARAMS, OA_FITTING_PARAMS
from lib.integration import Integration


logger = logging.getLogger(__name__)


class Fitting:
    """
    Handles manual and automatic fitting of Z-scan curves.

    Provides methods for computing the parametrized response function (manual fitting)
    and for automatic fitting using least-squares optimization.

    Attributes:
        sample_type (Integration): Integration object for the sample
        amplitude (float): Phase shift amplitude or transmittance
        beamwaist (float): Beam waist radius [m]
        zero_level (float): Normalized zero level (baseline)
        centerpoint (float): Center position of the curve [data points]
        nop (int): Number of data points
        ydata (NDArray): Measured y-values to fit to
    """

    def __init__(
        self,
        sample_type: Integration,
        amplitude: float,
        beamwaist: float,
        zero_level: float,
        centerpoint: float,
        nop: int,
        data: NDArray,
    ):
        """
        Initialize the Fitting object.

        Args:
            sample_type (Integration): Integration instance for calculations
            amplitude (float): Initial phase shift (DPhi0) or transmittance (T)
            beamwaist (float): Beam waist radius [m]
            zero_level (float): Baseline transmittance level
            centerpoint (float): Initial center position [data points]
            nop (int): Number of data points
            data (NDArray): Y-values to fit
        """
        self.sample_type = sample_type
        self.amplitude = amplitude
        self.beamwaist = beamwaist
        self.zero_level = zero_level
        self.centerpoint = centerpoint
        self.nop = nop
        self.ydata = data
        self.params: Parameters | None = None

    def manual(
        self,
        zero_level: float,
        centerpoint: float,
        amplitude: float,
        beamwaist: float,
        z_range: float,
        window: Optional[Any] = None,
        stype: str = "CA",
    ) -> NDArray:
        """
        Calculate fitted curve for given parameters (manual fitting).

        Computes the integrated field response for the given parameters
        without optimization.

        Args:
            zero_level (float): Baseline transmittance level
            centerpoint (float): Center position of curve [data points]
            amplitude (float): Phase shift (DPhi0) or transmittance (T)
            beamwaist (float): Beam waist radius [m]
            z_range (float): Total z-scan range [m]
            window (Any): Main window object (for accessing hardware parameters)
            stype (str): Measurement type ("CA" or "OA")

        Returns:
            NDArray: Calculated transmittance curve
        """
        self.z_range = z_range
        self.sample_type.z = np.array(
            [
                self.z_range * (zz - centerpoint) / self.nop - self.z_range / 2
                for zz in range(self.nop)
            ],
            dtype=float,
        )

        # =========================================================================
        # GET HARDWARE PARAMETERS
        # =========================================================================
        if window is None:
            # In headless mode, use values already stored in Integration object
            d0 = self.sample_type.d0
            ra = self.sample_type.ra
        else:
            # In UI mode, get current values from window
            try:
                window.get_general_parameters()
                d0 = window.d0
                ra = window.ra
            except Exception:
                # Fallback to Integration values if window method fails
                d0 = self.sample_type.d0
                ra = self.sample_type.ra

        # =========================================================================
        # COMPUTE FITTED CURVE
        # =========================================================================
        if stype == "CA":
            self.sample_type.derive(amplitude, beamwaist, d0, ra, stype)
            cas = self.sample_type.closed_sum
            oas = self.sample_type.open_sum
            result = (cas / oas) + (zero_level - 1)
        elif stype == "OA":
            self.sample_type.derive(amplitude, beamwaist, d0, ra, stype)
            oas = self.sample_type.Tznorm
            result = oas + (zero_level - 1)
        else:
            raise ValueError(f"Unknown measurement type: {stype}")

        return np.asarray(result, dtype=float)

    def fcn2min(
        self, params: Parameters, weights: NDArray[np.floating]
    ) -> NDArray[np.floating]:
        """
        Function to minimize during automatic fitting.

        Computes the weighted squared error between the measured data
        and the model prediction for the current parameter set.

        Args:
            params (Parameters):
                lmfit Parameters object containing the current parameter values.
            weights (array-like):
                1D array of weights for each data point (typically 0 or 1),
                used to emphasize or suppress selected regions of the curve.

        Returns:
            NDArray:
                Array of weighted squared errors for each data point.
                This vector is passed to the least-squares optimizer.
        """
        vals = params.valuesdict()
        ynew = self.manual(
            vals["Zero"],
            vals["Center"],
            vals.get("DPhi0", vals.get("T", 0.0)),
            vals["Beamwaist"],
            vals["Zrange"],
            window=self._fitting_window,
            stype=self._fitting_stype,
        )

        weights = np.asarray(weights, dtype=float).ravel()
        ynew = np.asarray(ynew, dtype=float).ravel()
        ydata = np.asarray(self.ydata, dtype=float).ravel()

        # Compute residuals and guard against NaNs/Infs in model output
        residuals = ynew - ydata
        if not np.all(np.isfinite(residuals)):
            # Replace non-finite residuals with a large penalty so optimizer steers away
            residuals = np.where(np.isfinite(residuals), residuals, 1e6)

        return np.asarray(weights * (residuals) ** 2, dtype=float)

    def automatic(
        self,
        z_range: float,
        ftype: str,
        stype: str,
        line_xydata: Tuple[NDArray, NDArray],
        window: Optional[Any] = None,
        vary_beamwaist: bool = True,
        vary_centerpoint: bool = True,
        max_iterations: int = 1000,
    ) -> Tuple[Any, NDArray]:
        """
        Automatically fit curve using least-squares optimization.

        Performs non-linear least-squares fitting with optional weighting
        based on cursor positions for region of interest.

        Args:
            z_range (float): Total z-scan range [m]
            ftype (str): Sample type ("Silica", "Solvent", "Sample")
            stype (str): Measurement type ("CA" or "OA")
            line_xydata (Tuple[NDArray, NDArray]): (x_data, y_data) tuple
            window (Optional[Any]): Main window object (for hardware parameters and cursors).
                        If None, operates in headless/testing mode.
            vary_beamwaist (bool): Allow beamwaist to vary during fit
            vary_centerpoint (bool): Allow center point to vary during fit
            max_iterations (int): Maximum optimization iterations

        Returns:
            Tuple[Any, NDArray]: (fitted_parameters, result_curve)
        """

        # =========================================================================
        # HEADLESS MODE DETECTION
        # =========================================================================
        headless = window is None

        # Store window and stype for use in fcn2min during optimization
        self._fitting_window = window
        self._fitting_stype = stype

        self.z_range = z_range
        xs, ys = line_xydata
        weights: NDArray[np.floating] = np.ones(xs.shape, dtype=float)

        # =========================================================================
        # WEIGHTING FROM CURSORS (UI MODE ONLY)
        # =========================================================================
        # Only attempt to get cursor positions in UI mode
        if not headless:
            cursor_attribute = self._get_cursor_attribute(ftype, stype)
            fix_attr = f"{ftype.lower()}{stype}_fixROI_checkBox"
            try:
                if (
                    hasattr(window, fix_attr)
                    and getattr(window, fix_attr).isChecked()
                ):
                    if hasattr(window, cursor_attribute):
                        cursor_positions = getattr(window, cursor_attribute)
                        if len(cursor_positions) == 2:
                            weights = self._calculate_weights(
                                xs, cursor_positions
                            )
            except Exception:
                # On any UI/widget error, fall back to full weights
                weights = np.ones(np.shape(xs))

        # =========================================================================
        # SET UP FITTING PARAMETERS
        # =========================================================================
        self.params = Parameters()

        # Choose bounds from config where available to improve stability
        try:
            if stype == "CA":
                zero_min = CA_FITTING_PARAMS["Zero"]["min"]
                zero_max = CA_FITTING_PARAMS["Zero"]["max"]
                center_min = CA_FITTING_PARAMS["Center"]["min"]
                center_max = CA_FITTING_PARAMS["Center"]["max"]
                dphi0_min = CA_FITTING_PARAMS["DPhi0"]["min"]
                dphi0_max = CA_FITTING_PARAMS["DPhi0"]["max"]
                bw_min = CA_FITTING_PARAMS["Beamwaist"]["min"]
                bw_max = CA_FITTING_PARAMS["Beamwaist"]["max"]
            else:
                zero_min = OA_FITTING_PARAMS["Zero"]["min"]
                zero_max = OA_FITTING_PARAMS["Zero"]["max"]
                center_min = OA_FITTING_PARAMS["Center"]["min"]
                center_max = OA_FITTING_PARAMS["Center"]["max"]
                dphi0_min = OA_FITTING_PARAMS.get("DPhi0", {}).get("min", -2)
                dphi0_max = OA_FITTING_PARAMS.get("DPhi0", {}).get("max", 2)
                bw_min = OA_FITTING_PARAMS["Beamwaist"]["min"]
                bw_max = OA_FITTING_PARAMS["Beamwaist"]["max"]
        except Exception:
            # Fallback defaults
            zero_min, zero_max = 0.75, 1.25
            center_min, center_max = -50, 50
            dphi0_min, dphi0_max = -2, 2
            bw_min, bw_max = 15e-6, 150e-6

        # Common params
        self.params.add(
            "Zero", value=self.zero_level, min=zero_min, max=zero_max
        )
        self.params.add(
            "Center",
            value=self.centerpoint,
            min=center_min,
            max=center_max,
            vary=vary_centerpoint,
        )

        # Debug: log the Center parameter setup for OA fits
        if stype == "OA":
            logger.debug(
                f"OA Fitting: Center initial={self.centerpoint}, vary={vary_centerpoint}, bounds=[{center_min}, {center_max}]"
            )

        if stype == "CA":
            self.params.add(
                "DPhi0",
                value=self.amplitude if self.amplitude is not None else 0.0,
                min=dphi0_min,
                max=dphi0_max,
            )
            self.params.add(
                "Beamwaist",
                value=self.beamwaist,
                min=bw_min,
                max=bw_max,
                vary=vary_beamwaist,
            )
            self.params.add("Zrange", value=self.z_range, vary=False)
        elif stype == "OA":
            # T bounds from config
            try:
                t_min = OA_FITTING_PARAMS["T"]["min"]
                t_max = OA_FITTING_PARAMS["T"]["max"]
            except Exception:
                t_min, t_max = -2, 2

            # Use data-driven initial guess for OA amplitude if not provided
            try:
                ys_arr = np.asarray(ys, dtype=float).ravel()
                guess_T = (
                    self.amplitude
                    if (self.amplitude is not None and self.amplitude != 0)
                    else max(1.0 - np.min(ys_arr), 1e-6)
                )
            except Exception:
                guess_T = self.amplitude if self.amplitude is not None else 0.0

            # Clamp initial guess to bounds
            try:
                guess_T = max(min(guess_T, t_max), t_min)
            except Exception:
                pass

            self.params.add("T", value=float(guess_T), min=t_min, max=t_max)
            self.params.add(
                "Beamwaist",
                value=self.beamwaist,
                min=bw_min,
                max=bw_max,
                vary=False,
            )
            self.params.add("Zrange", value=self.z_range, vary=False)

        # =========================================================================
        # PERFORM FITTING
        # =========================================================================
        fitter = Minimizer(self.fcn2min, self.params, fcn_args=(weights,))
        
        # Try least_squares first (better for covariance calculation)
        result = fitter.minimize(
            method="least_squares", 
            max_nfev=max_iterations
        )
        
        # If least_squares fails to calculate covariance, try leastsq
        if result.covar is None:
            logger.warning(
                f"Covariance matrix not calculated with least_squares for {ftype} {stype}. "
                "Attempting with leastsq method..."
            )
            result = fitter.minimize(
                method="leastsq",
                max_nfev=max_iterations
            )
        
        # Final fallback: if still no covariance, use Nelder-Mead for final polish
        if result.covar is None:
            logger.warning(
                f"Still no covariance for {ftype} {stype}. Attempting Nelder-Mead polish..."
            )
            # Re-initialize with last best result
            self.params = result.params
            fitter = Minimizer(self.fcn2min, self.params, fcn_args=(weights,))
            result = fitter.minimize(
                method="nelder",
                max_nfev=max_iterations
            )

        # =========================================================================
        # CALCULATE RESULT CURVE WITH FITTED PARAMETERS
        # =========================================================================
        vals = result.params.valuesdict()

        result_line = self.manual(
            zero_level=vals["Zero"],
            centerpoint=vals["Center"],
            amplitude=vals.get("DPhi0", vals.get("T")),
            beamwaist=vals["Beamwaist"],
            z_range=vals["Zrange"],
            window=window,
            stype=stype,
        )

        # =========================================================================
        # PRINT RESULTS (UI MODE ONLY)
        # =========================================================================
        # Only print in UI mode to avoid spam in test output
        if not headless:
            result.params.pretty_print()
            
            # Log convergence info
            if result.covar is None:
                logger.warning(
                    f"{ftype} {stype}: Fit converged but covariance matrix is None. "
                    f"Errors will be unreliable. Success={result.success}, "
                    f"Method used last."
                )
            else:
                logger.info(
                    f"{ftype} {stype}: Fit successful with covariance matrix calculated."
                )

        return result, result_line

    @staticmethod
    def _get_cursor_attribute(ftype: str, stype: str) -> str:
        """
        Get the attribute name for cursor positions.

        Args:
            ftype (str): Sample type ("Silica", "Solvent", "Sample")
            stype (str): Measurement type ("CA" or "OA")

        Returns:
            str: Attribute name in main window object
        """
        if ftype == "Silica":
            return "silicaCA_cursorPositions"
        elif ftype == "Solvent":
            if stype == "CA":
                return "solventCA_cursorPositions"
            else:
                return "solventOA_cursorPositions"
        elif ftype == "Sample":
            if stype == "CA":
                return "sampleCA_cursorPositions"
            else:
                return "sampleOA_cursorPositions"
        else:
            raise ValueError(f"Unknown sample type: {ftype}")

    @staticmethod
    def _calculate_weights(
        xs: NDArray, cursor_positions: List[Tuple[float, float]]
    ) -> NDArray[np.floating]:
        """
        Calculate data point weights based on cursor positions.

        Sets weight to 1 for points between cursors, 0 elsewhere.

        Args:
            xs (NDArray): X-axis values (positions)
            cursor_positions (List[Tuple[float, float]]): Two cursor points

        Returns:
            NDArray[np.floating]: Weights for each point
        """
        x1 = cursor_positions[0][0]
        x2 = cursor_positions[1][0]

        x1_index = list(xs).index(min(xs, key=lambda x: abs(x - x1)))
        x2_index = list(xs).index(min(xs, key=lambda x: abs(x - x2)))

        x_sm, x_lg = sorted([x1_index, x2_index])
        weights = [
            0 if (xi < x_sm or xi > x_lg) else 1 for xi in range(len(xs))
        ]

        return np.asarray(weights, dtype=float)