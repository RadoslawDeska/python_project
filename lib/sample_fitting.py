"""
Sample fitting utilities for Z-Scan measurements.

This module provides specialized functionality for fitting Z-scan curves
of unknown samples, building on solvent baseline measurements.
"""

from typing import Any, Optional
import numpy as np
from numpy.typing import NDArray

from lib.integration import Integration
from lib.fitting import Fitting


class SampleFitter:
    """
    Handles fitting of sample Z-scan curves.
    
    Manages the fitting workflow for unknown samples, including:
    - Validation that silica reference and solvent baseline are fitted
    - Automatic and manual fitting procedures
    - Parameter extraction and error calculation
    
    Attributes:
        silica_fitted (bool): Whether silica reference has been fitted
        solvent_fitted (bool): Whether solvent baseline has been fitted
    """
    
    def __init__(self):
        """Initialize the SampleFitter."""
        self.silica_fitted = False
        self.solvent_fitted = False
    
    def check_prerequisites(
        self,
        silica_fitted: bool,
        solvent_fitted: bool
    ) -> tuple[bool, str]:
        """
        Verify that prerequisite fittings are complete.
        
        Checks that both silica reference and solvent baseline have
        been successfully fitted before proceeding with sample fitting.
        
        Args:
            silica_fitted (bool): Whether silica CA fitting is complete
            solvent_fitted (bool): Whether solvent CA fitting is complete
        
        Returns:
            tuple[bool, str]: (is_valid, message)
                - is_valid: True if all prerequisites are met
                - message: Error or status message
        """
        if not silica_fitted:
            return False, "Please fit silica reference first."
        if not solvent_fitted:
            return False, "Please fit solvent baseline first."
        return True, "Prerequisites met."
    
    def fit_automatically_ca(
        self,
        sample_curves: Integration,
        initial_dphi0: float,
        initial_beamwaist: float,
        initial_zero_level: float,
        initial_centerpoint: float,
        nop: int,
        z_range: float,
        line_data: tuple[NDArray, NDArray],
        window: Any,
        solvent_dphi0: Optional[float] = None
    ) -> tuple[Any, NDArray]:
        """
        Automatically fit sample closed aperture (CA) curve.
        
        Performs non-linear least-squares fitting of the sample CA data.
        Can optionally fix the DPhi0 parameter to solvent value for consistency.
        
        Args:
            sample_curves (Integration): Integration object for sample
            initial_dphi0 (float): Initial phase shift estimate [rad]
            initial_beamwaist (float): Initial beam waist [m]
            initial_zero_level (float): Initial baseline transmittance
            initial_centerpoint (float): Initial center position [points]
            nop (int): Number of data points
            z_range (float): Total z-scan range [m]
            line_data (tuple): (x_data, y_data) from plot
            window (Any): Main window object (for parameters and cursors)
            solvent_dphi0 (float, optional): Solvent DPhi0 to use as constraint
        
        Returns:
            tuple: (fitted_parameters, result_curve)
        
        Raises:
            ValueError: If solvent DPhi0 is required but not provided
        """
        # Create fitter object
        sample_calculation = Fitting(
            sample_curves,
            initial_dphi0,
            initial_beamwaist,
            initial_zero_level,
            initial_centerpoint,
            nop,
            line_data[1]  # y-data
        )
        
        # Determine if beamwaist should vary
        vary_beamwaist = window.sampleCA_customBeamwaist_checkBox.isChecked()
        
        # Perform automatic fitting
        minimizer_result, result = sample_calculation.automatic(
            z_range=z_range,
            ftype="Sample",
            stype="CA",
            line_xydata=line_data,
            window=window,
            vary_beamwaist=vary_beamwaist,
            vary_centerpoint=True
        )
        
        return minimizer_result, result
    
    def fit_automatically_oa(
        self,
        sample_curves: Integration,
        initial_transmittance: float,
        initial_beamwaist: float,
        initial_zero_level: float,
        initial_centerpoint: float,
        nop: int,
        z_range: float,
        line_data: tuple[NDArray, NDArray],
        window: Any,
        absorption_model: str = "2PA"
    ) -> tuple[Any, NDArray]:
        """
        Automatically fit sample open aperture (OA) curve.
        
        Performs non-linear least-squares fitting of the sample OA data
        for nonlinear absorption characterization.
        
        Args:
            sample_curves (Integration): Integration object for sample
            initial_transmittance (float): Initial transmittance estimate
            initial_beamwaist (float): Beam waist from CA fitting [m]
            initial_zero_level (float): Initial baseline transmittance
            initial_centerpoint (float): Initial center position [points]
            nop (int): Number of data points
            z_range (float): Total z-scan range [m]
            line_data (tuple): (x_data, y_data) from plot
            window (Any): Main window object (for parameters)
            absorption_model (str): Model for nonlinear absorption
        
        Returns:
            tuple: (fitted_parameters, result_curve)
        """
        # Create fitter object
        sample_calculation = Fitting(
            sample_curves,
            initial_transmittance,
            initial_beamwaist,
            initial_zero_level,
            initial_centerpoint,
            nop,
            line_data[1]  # y-data
        )
        
        # Determine center point variation
        vary_centerpoint = window.sampleOA_customCenterPoint_checkBox.isChecked()
        
        # Perform automatic fitting
        minimizer_result, result = sample_calculation.automatic(
            z_range=z_range,
            ftype="Sample",
            stype="OA",
            line_xydata=line_data,
            window=window,
            vary_beamwaist=False,  # OA doesn't fit beamwaist
            vary_centerpoint=vary_centerpoint
        )
        
        return minimizer_result, result
    
    def fit_manually_ca(
        self,
        sample_curves: Integration,
        dphi0: float,
        beamwaist: float,
        zero_level: float,
        centerpoint: float,
        nop: int,
        z_range: float,
        window: Any
    ) -> NDArray:
        """
        Calculate sample CA curve with manual parameters (no optimization).
        
        Used during interactive fitting with sliders to preview curve
        with current parameter values.
        
        Args:
            sample_curves (Integration): Integration object for sample
            dphi0 (float): Phase shift [rad]
            beamwaist (float): Beam waist [m]
            zero_level (float): Baseline transmittance
            centerpoint (float): Center position [points]
            nop (int): Number of data points
            z_range (float): Total z-scan range [m]
            window (Any): Main window object
        
        Returns:
            NDArray: Calculated curve
        """
        fitter = Fitting(
            sample_curves,
            dphi0,
            beamwaist,
            zero_level,
            centerpoint,
            nop,
            np.zeros(nop)  # Dummy data
        )
        
        result = fitter.manual(
            zero_level=zero_level,
            centerpoint=centerpoint,
            amplitude=dphi0,
            beamwaist=beamwaist,
            z_range=z_range,
            window=window,
            stype="CA"
        )
        
        return result
    
    def fit_manually_oa(
        self,
        sample_curves: Integration,
        transmittance: float,
        beamwaist: float,
        zero_level: float,
        centerpoint: float,
        nop: int,
        z_range: float,
        window: Any,
        absorption_model: str = "2PA"
    ) -> NDArray:
        """
        Calculate sample OA curve with manual parameters (no optimization).
        
        Used during interactive fitting with sliders to preview OA curve
        with current parameter values.
        
        Args:
            sample_curves (Integration): Integration object for sample
            transmittance (float): Nonlinear transmittance
            beamwaist (float): Beam waist [m]
            zero_level (float): Baseline transmittance
            centerpoint (float): Center position [points]
            nop (int): Number of data points
            z_range (float): Total z-scan range [m]
            window (Any): Main window object
            absorption_model (str): Nonlinear absorption model
        
        Returns:
            NDArray: Calculated curve
        """
        # Update absorption model if needed
        if hasattr(sample_curves, 'calculate_Tz_for_OA'):
            sample_curves.calculate_Tz_for_OA(model=absorption_model)
        
        fitter = Fitting(
            sample_curves,
            transmittance,
            beamwaist,
            zero_level,
            centerpoint,
            nop,
            np.zeros(nop)  # Dummy data
        )
        
        result = fitter.manual(
            zero_level=zero_level,
            centerpoint=centerpoint,
            amplitude=transmittance,
            beamwaist=beamwaist,
            z_range=z_range,
            window=window,
            stype="OA"
        )
        
        return result
    
    @staticmethod
    def extract_ca_parameters(minimizer_result: Any) -> dict:
        """
        Extract and organize CA fitting parameters.
        
        Args:
            minimizer_result: lmfit Minimizer result object
        
        Returns:
            dict: Extracted parameters with values and uncertainties
        """
        params = minimizer_result.params
        return {
            'dphi0': {
                'value': params['DPhi0'].value,
                'stderr': params['DPhi0'].stderr
            },
            'beamwaist': {
                'value': params['Beamwaist'].value,
                'stderr': params['Beamwaist'].stderr
            },
            'zero_level': {
                'value': params['Zero'].value,
                'stderr': params['Zero'].stderr
            },
            'centerpoint': {
                'value': params['Center'].value,
                'stderr': params['Center'].stderr
            }
        }
    
    @staticmethod
    def extract_oa_parameters(minimizer_result: Any) -> dict:
        """
        Extract and organize OA fitting parameters.
        
        Args:
            minimizer_result: lmfit Minimizer result object
        
        Returns:
            dict: Extracted parameters with values and uncertainties
        """
        params = minimizer_result.params
        return {
            'transmittance': {
                'value': params['T'].value,
                'stderr': params['T'].stderr
            },
            'zero_level': {
                'value': params['Zero'].value,
                'stderr': params['Zero'].stderr
            },
            'centerpoint': {
                'value': params['Center'].value,
                'stderr': params['Center'].stderr
            }
        }
