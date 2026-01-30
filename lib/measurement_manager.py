"""Measurement Manager - Handles hardware measurement and data acquisition.

Manages motor control, DAQ hardware interaction, and measurement workflow.
"""

import time
from typing import TYPE_CHECKING, Any, Optional

import numpy as np
from numpy.typing import NDArray

winsound: Any
try:
    import winsound
except ImportError:
    # winsound is Windows-only
    winsound = None

try:
    import nidaqmx
    from nidaqmx.constants import AcquisitionType, Edge
except ImportError:
    # nidaqmx may not be available in all environments
    nidaqmx = None


if TYPE_CHECKING:
    from lib.mgmotor import MG17Motor
else:
    MG17Motor = Any


class MeasurementManager:
    """Manages hardware measurement and motor positioning."""

    def __init__(self, motor: Optional["MG17Motor"] = None):
        """Initialize measurement manager.

        Args:
            motor: Optional MG17Motor instance for Z-positioning
        """
        self.motor = motor
        self.where_to_start: Optional[str] = None
        self.measurement_stopped = False

    def setup_motor(self, motor: "MG17Motor"):
        """Set or replace the motor instance.

        Args:
            motor: MG17Motor instance
        """
        self.motor = motor

    def move_to_start(self, start_pos: float, end_pos: float):
        """Move motor to starting position.

        Args:
            start_pos: Starting position
            end_pos: Ending position
        """
        if self.motor is None:
            raise RuntimeError("Motor not configured")

        mid_point = (start_pos + end_pos) / 2
        if self.motor.position <= mid_point:
            self.motor.move_to(start_pos)
            self.where_to_start = "start"
        else:
            self.motor.move_to(end_pos)
            self.where_to_start = "end"

        if self.motor.is_in_motion:
            print("Moving to starting position")

        while self.motor.is_in_motion:
            time.sleep(0.01)

    def move_to_custom_position(self, target: float):
        """Move motor to a custom position.

        Args:
            target: Target position in mm

        Raises:
            ValueError: If target > 100 mm
        """
        if self.motor is None:
            raise RuntimeError("Motor not configured")

        if target > 100:
            raise ValueError("Target position exceeds maximum (100 mm)")

        self.motor.move_to(target)
        while self.motor.is_in_motion:
            time.sleep(0.01)

    def move_by_step(self, start_pos: float, end_pos: float, num_steps: int):
        """Move motor by a single step.

        Args:
            start_pos: Starting position
            end_pos: Ending position
            num_steps: Total number of steps
        """
        if self.motor is None:
            raise RuntimeError("Motor not configured")

        step_size = (end_pos - start_pos) / num_steps

        if self.where_to_start == "end":
            self.motor.move_by(-step_size, blocking=True)
        else:
            self.motor.move_by(step_size, blocking=True)

        time.sleep(0.05)  # Allow motor to settle

    def home_motor(self):
        """Perform homing operation on motor.

        Returns:
            Status message
        """
        if self.motor is None:
            raise RuntimeError("Motor not configured")

        if not self.motor.has_homing_been_completed:
            self.motor.move_home()
            time.sleep(0.2)

            if self.motor.is_in_motion:
                print("Homing now")

            while not self.motor.has_homing_been_completed:
                time.sleep(0.01)

            time.sleep(0.2)

        return "Homing performed!"

    def get_motor_position(self) -> float:
        """Get current motor position.

        Returns:
            Current position in mm
        """
        if self.motor is None:
            return 0.0
        return float(self.motor.position)

    def acquire_data(
        self,
        num_channels: int,
        samples_per_step: int,
        detector_core_name: str = "Dev1/ai",
        trigger_src: str = "/Dev1/PFI0",
    ) -> NDArray:
        """Acquire data from DAQ hardware.

        Args:
            num_channels: Number of channels to read
            samples_per_step: Samples per channel
            detector_core_name: DAQ device channel prefix
            trigger_src: Digital trigger source

        Returns:
            Array of acquired data (num_channels x samples_per_step)
        """
        if nidaqmx is None:
            raise RuntimeError(
                "nidaqmx not available - install with: pip install nidaqmx"
            )

        with nidaqmx.Task() as task:
            # Create multichannel input
            task.ai_channels.add_ai_voltage_chan(
                f"{detector_core_name}0:{num_channels}"
            )

            # Setup digital trigger
            task.triggers.start_trigger.dig_edge_src = trigger_src

            # Setup sample clock (1 kHz acquisition rate)
            rate = 1000
            task.timing.cfg_samp_clk_timing(
                rate,
                source=trigger_src,
                active_edge=Edge.FALLING,
                sample_mode=AcquisitionType.FINITE,
                samps_per_chan=samples_per_step,
            )

            # Create reader and data buffer
            reader = nidaqmx.stream_readers.AnalogMultiChannelReader(
                task.in_stream
            )
            values_read = np.zeros(
                (num_channels + 1, samples_per_step), dtype=np.float64
            )

            # Start task and acquire data
            task.start()
            time.sleep(0.2)  # Allow settling time
            reader.read_many_sample(
                values_read, number_of_samples_per_channel=samples_per_step
            )
            task.stop()

        return values_read

    def beep_completion(
        self, frequency: int = 800, duration: int = 500, count: int = 3
    ):
        """Emit audio signal for measurement completion.

        Args:
            frequency: Beep frequency in Hz
            duration: Beep duration in milliseconds
            count: Number of beeps
        """
        if winsound is None:
            # Not on Windows - print message instead
            print(f"[BEEP] {frequency} Hz for {duration} ms (x{count})")
            return

        for _ in range(count):
            try:
                winsound.Beep(frequency, duration)
            except Exception:
                print(f"Beep: {frequency} Hz for {duration} ms")
            time.sleep(0.05)

    def stop_measurement(self):
        """Signal to stop ongoing measurement."""
        self.measurement_stopped = True

    def reset_stop_flag(self):
        """Reset the stop measurement flag."""
        self.measurement_stopped = False

    def is_measurement_stopped(self) -> bool:
        """Check if measurement stop has been requested.

        Returns:
            True if stop has been requested
        """
        return self.measurement_stopped
