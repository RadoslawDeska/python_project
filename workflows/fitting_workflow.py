# workflows/fitting_workflow.py
"""
Pure business logic for fitting workflows.
No UI dependencies - testable and reusable.
"""

from typing import Tuple, Optional
import numpy as np
from numpy.typing import NDArray
from dataclasses import dataclass

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


class SolventOAFittingWorkflow:
    """Encapsulates solvent OA fitting logic."""
    
    def __init__(self, fitting_manager, integration_manager, state: StateManager):
        """
        Args:
            fitting_manager: FittingManager instance
            integration_manager: Integration instance
            state: StateManager instance
        """
        self.fitting = fitting_manager
        self.integration = integration_manager
        self.state = state
    
    def fit_manual(
        self,
        z_positions: NDArray,
        oa_data: NDArray,
        beamwaist: float,
        zero_level: float,
        center_point: float,
        t_value: float,
        z_range: float,
        **general_params
    ) -> Optional[FittingResult]:
        """Manually fit OA data.
        
        Args:
            z_positions: Z-scan position array
            oa_data: Open aperture intensity data
            beamwaist: Beam waist radius [m]
            zero_level: Baseline transmittance
            center_point: Center position [data points]
            t_value: Transmittance/amplitude parameter
            z_range: Total z-scan range [m]
            **general_params: Additional parameters (wavelength, etc.)
            
        Returns:
            FittingResult if successful, None if skipped
        """
        # Validate inputs
        valid, error = ParameterValidator.validate_oa_parameters(t_value)
        if not valid:
            raise ValueError(f"Invalid OA parameters: {error}")
        
        # Check if fitting should be skipped
        if self.state.solventOA_no_absorption:
            return None  # Fitting skipped
        
        try:
            # Create Integration object for OA
            integration = Integration(
                beta=t_value,
                n2=general_params.get('n2', 0),
                DPhi0=0,
                positions=z_positions,
                d0=general_params.get('d0', 0),
                aperture_radius=general_params.get('ra', 0),
                wavelength=general_params.get('wavelength', 1e-6),
                beamwaist=beamwaist,
                n_components=general_params.get('n_components', 8),
                integration_steps=general_params.get('integration_steps', 30),
                stype="OA"
            )
            
            # Fit
            fitter = Fitting(
                integration, t_value, beamwaist, zero_level, 
                center_point, len(z_positions), oa_data
            )
            
            result_curve = fitter.manual(
                zero_level,
                center_point,
                t_value,
                beamwaist,
                z_range,
                window=None,   # headless mode
                stype="OA",
            )

            # Update state
            self.state.solventOA_autofit_done = False

            return FittingResult(
                curve=result_curve,
                params={
                    "t_value": t_value,
                    "beamwaist": beamwaist,
                    "zero_level": zero_level,
                    "center_point": center_point,
                },
                errors={},
                converged=True,
            )

        except Exception as e:
            return FittingResult(
                curve=np.array([]),
                params={},
                errors={},
                converged=False,
                message=f"OA manual fitting failed: {str(e)}",
            )
