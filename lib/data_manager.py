"""Data Manager - Handles data storage, calculations, and processing.

Manages measurement data, parameter calculations, and error propagation.
"""

from typing import Dict, List, TypedDict

import numpy as np
from numpy.typing import NDArray


class MeasurementData(TypedDict):
    positions: List[float]
    absolute: List[List[float]]
    relative: List[List[float]]


class DataManager:
    """Manages application data storage and calculations."""

    def __init__(self) -> None:
        """Initialize data manager."""

        self.data: MeasurementData = {
            "positions": [],
            "absolute": [[], [], [], []],
            "relative": [[], [], [], []],
        }

        self.rms_value: float = 0.0

        # Fitting parameters
        self.silica_params: Dict[str, float] = {}
        self.solvent_params: Dict[str, float] = {}
        self.solvent_ca_params: Dict[str, float] = {}
        self.solvent_oa_params: Dict[str, float] = {}
        self.sample_ca_params: Dict[str, float] = {}
        self.sample_oa_params: Dict[str, float] = {}

        # Fitting parameter errors
        self.silica_params_errors: Dict[str, float] = {}
        self.solvent_params_errors: Dict[str, float] = {}
        self.solvent_ca_params_errors: Dict[str, float] = {}
        self.solvent_oa_params_errors: Dict[str, float] = {}
        self.sample_ca_params_errors: Dict[str, float] = {}
        self.sample_oa_params_errors: Dict[str, float] = {}

    def clear_measurement_data(self) -> None:
        """Clear all measurement data."""
        self.data = {
            "positions": [],
            "absolute": [[], [], [], []],
            "relative": [[], [], [], []],
        }

        self.rms_value = 0.0

    def add_measurement_point(self, position: float, values: NDArray) -> None:
        """Add a measurement point to the data.

        Args:
            position: Motor position
            values: Array of channel values
        """
        self.data["positions"].append(position)
        for chan_no in range(len(values)):
            self.data["absolute"][chan_no].append(values[chan_no])

    def add_relative_values(self, num_channels: int) -> None:
        """Calculate and add relative values (normalized by reference channel).

        Args:
            num_channels: Number of channels to process
        """
        step = len(self.data["positions"]) - 1
        for chan_no in range(num_channels):
            relative_value = (
                self.data["absolute"][chan_no][step]
                / self.data["absolute"][1][step]
            )
            self.data["relative"][chan_no].append(relative_value)

    def calculate_rms_noise(self, num_channels: int) -> float:
        """Calculate RMS noise from reference channel data."""
        y: List[float] = self.data["absolute"][1]

        if len(y) < 2:
            return 0.0

        mean_sq: float = float(np.mean([yi**2 for yi in y]))
        rms: float = float(np.abs(np.sqrt(mean_sq) - y[0]) / y[0])
        self.rms_value = rms
        return rms

    def store_fitting_parameters(
        self, sample_type: str, ca_or_oa: str, params: Dict
    ) -> None:
        """Store fitting parameters for a sample.

        Args:
            sample_type: "silica", "solvent", or "sample"
            ca_or_oa: "CA" or "OA"
            params: Dictionary of parameters
        """
        if sample_type == "silica":
            self.silica_params = params
        elif sample_type == "solvent":
            if ca_or_oa == "CA":
                self.solvent_ca_params = params
            else:
                self.solvent_oa_params = params
        elif sample_type == "sample":
            if ca_or_oa == "CA":
                self.sample_ca_params = params
            else:
                self.sample_oa_params = params

    def store_fitting_parameter_errors(
        self, sample_type: str, ca_or_oa: str, params_errors: Dict
    ) -> None:
        """Store fitting parameter errors for a sample.

        Args:
            sample_type: "silica", "solvent", or "sample"
            ca_or_oa: "CA" or "OA"
            params_errors: Dictionary of parameter errors
        """
        if sample_type == "silica":
            self.silica_params_errors = params_errors
        elif sample_type == "solvent":
            if ca_or_oa == "CA":
                self.solvent_ca_params_errors = params_errors
            else:
                self.solvent_oa_params_errors = params_errors
        elif sample_type == "sample":
            if ca_or_oa == "CA":
                self.sample_ca_params_errors = params_errors
            else:
                self.sample_oa_params_errors = params_errors

    def get_fitting_parameters(self, sample_type: str, ca_or_oa: str) -> Dict:
        """Get stored fitting parameters for a sample.

        Args:
            sample_type: "silica", "solvent", or "sample"
            ca_or_oa: "CA" or "OA"

        Returns:
            Dictionary of parameters
        """
        if sample_type == "silica":
            return self.silica_params
        elif sample_type == "solvent":
            if ca_or_oa == "CA":
                return self.solvent_ca_params
            else:
                return self.solvent_oa_params
        elif sample_type == "sample":
            if ca_or_oa == "CA":
                return self.sample_ca_params
            else:
                return self.sample_oa_params
        return {}

    def get_fitting_parameter_errors(
        self, sample_type: str, ca_or_oa: str
    ) -> Dict:
        """Get stored fitting parameter errors for a sample.

        Args:
            sample_type: "silica", "solvent", or "sample"
            ca_or_oa: "CA" or "OA"

        Returns:
            Dictionary of parameter errors
        """
        if sample_type == "silica":
            return self.silica_params_errors
        elif sample_type == "solvent":
            if ca_or_oa == "CA":
                return self.solvent_ca_params_errors
            else:
                return self.solvent_oa_params_errors
        elif sample_type == "sample":
            if ca_or_oa == "CA":
                return self.sample_ca_params_errors
            else:
                return self.sample_oa_params_errors
        return {}

    def get_measurement_data(
        self, data_type: str = "relative"
    ) -> Dict[int, List[float]]:
        """Get measurement data by type.

        Args:
            data_type: "absolute" or "relative"

        Returns:
            Dictionary mapping channel numbers to data lists
        """
        if data_type == "absolute":
            source = self.data["absolute"]
        elif data_type == "relative":
            source = self.data["relative"]
        else:
            raise ValueError(f"Unknown data_type: {data_type}")

        result: Dict[int, List[float]] = {}
        for chan_no, values in enumerate(source):
            result[chan_no] = values
        return result

    def reverse_measurement_data(self) -> None:
        """Reverse all measurement data (for backward scans)."""
        self.data["positions"].reverse()
        for channel_data in self.data["absolute"]:
            channel_data.reverse()
        for channel_data in self.data["relative"]:
            channel_data.reverse()
