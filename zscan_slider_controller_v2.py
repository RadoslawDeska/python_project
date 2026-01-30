"""
Slider controller for Z-scan UI
Maps slider movements to fitting parameters
Connects to UI elements by name
"""

from dataclasses import dataclass
from typing import Dict, Callable, Optional, Tuple
from enum import Enum

import numpy as np


class SampleType(Enum):
    SILICA = "silica"
    SOLVENT = "solvent"
    SAMPLE = "sample"


class ApertureType(Enum):
    CA = "CA"
    OA = "OA"


@dataclass
class SliderConfig:
    """Configuration for a single slider"""

    slider_name: str  # UI widget name
    spinbox_name: Optional[str]  # Corresponding display spinbox
    param_name: str  # Parameter: amplitude, center, zero_level, etc
    min_val: float  # Physical value in non-prefixed SI base unit
    max_val: float  # Physical value in non-prefixed SI base unit
    step: float
    display_format: str  # Format string for display


class SliderController:
    """
    Manages slider-to-parameter mapping
    Maps actual slider ranges to physical quantities
    """

    # Define slider configs for each sample type + aperture
    SLIDER_CONFIGS = {
        # SILICA CA
        ("silica", "CA"): {
            "zero_level": SliderConfig(
                slider_name="silicaCA_zeroLevel_slider",
                spinbox_name="silicaCA_zeroLevel_doubleSpinBox",
                param_name="zero_level",
                min_val=0.75,  # no units
                max_val=1.25,  # no units
                step=0.01,
                display_format="{:.4f}",
            ),
            "amplitude": SliderConfig(
                slider_name="silicaCA_DPhi0_slider",
                spinbox_name="silicaCA_deltaPhi0Summary_doubleSpinBox",
                param_name="DPhi0",  # in radians
                min_val=0.0,
                max_val=2.0,
                step=0.01,
                display_format="{:.4f}",
            ),
            "centerpoint": SliderConfig(
                slider_name="silicaCA_centerPoint_slider",
                spinbox_name=None,
                param_name="centerpoint",  # in meters
                min_val=-0.01,
                max_val=0.01,
                step=0.001,
                display_format="{:.1f}",
            ),
            "beamwaist": SliderConfig(
                slider_name="silicaCA_Beamwaist_slider",
                spinbox_name="silicaCA_beamwaistSummary_doubleSpinBox",
                param_name="beamwaist",  # in meters
                min_val=5e-6,
                max_val=150e-6,
                step=1,
                display_format="{:.2f}",
            ),
        },
        # SOLVENT CA
        ("solvent", "CA"): {
            "zero_level": SliderConfig(
                slider_name="solventCA_zeroLevel_slider",
                spinbox_name=None,
                param_name="zero_level",
                min_val=0.75,
                max_val=1.25,
                step=0.01,
                display_format="{:.4f}",
            ),
            "amplitude": SliderConfig(
                slider_name="solventCA_DPhi0_slider",
                spinbox_name="solventCA_deltaPhi0Summary_doubleSpinBox",
                param_name="DPhi0",
                min_val=-2.0,
                max_val=2.0,
                step=0.01,
                display_format="{:.4f}",
            ),
            "centerpoint": SliderConfig(
                slider_name="solventCA_centerPoint_slider",
                spinbox_name=None,
                param_name="centerpoint",
                min_val=-50,
                max_val=50,
                step=1.0,
                display_format="{:.1f}",
            ),
            "beamwaist": SliderConfig(
                slider_name="solventCA_Beamwaist_slider",
                spinbox_name="solventCA_beamwaistSummary_doubleSpinBox",
                param_name="beamwaist",
                min_val=0.0,
                max_val=1.0,
                step=0.1,
                display_format="{:.2f}",
            ),
        },
        # SAMPLE CA
        ("sample", "CA"): {
            "zero_level": SliderConfig(
                slider_name="sampleCA_zeroLevel_slider",
                spinbox_name=None,
                param_name="zero_level",
                min_val=0.75,
                max_val=1.25,
                step=0.01,
                display_format="{:.4f}",
            ),
            "amplitude": SliderConfig(
                slider_name="sampleCA_DPhi0_slider",
                spinbox_name="sampleCA_deltaPhi0Summary_doubleSpinBox",
                param_name="DPhi0",
                min_val=-2.0,
                max_val=2.0,
                step=0.01,
                display_format="{:.4f}",
            ),
            "centerpoint": SliderConfig(
                slider_name="sampleCA_centerPoint_slider",
                spinbox_name=None,
                param_name="centerpoint",
                min_val=-50,
                max_val=50,
                step=1.0,
                display_format="{:.1f}",
            ),
            "beamwaist": SliderConfig(
                slider_name="sampleCA_Beamwaist_slider",
                spinbox_name="sampleCA_beamwaistSummary_doubleSpinBox",
                param_name="beamwaist",
                min_val=0.0,
                max_val=1.0,
                step=0.1,
                display_format="{:.2f}",
            ),
        },
        # SILICA OA
        ("silica", "OA"): {},  # No sliders for OA in UI yet
        # SOLVENT OA
        ("solvent", "OA"): {
            "zero_level": SliderConfig(
                slider_name="solventOA_zeroLevel_slider",
                spinbox_name="solventOA_zeroLevel_doubleSpinBox",
                param_name="zero_level",
                min_val=0.75,
                max_val=1.25,
                step=0.01,
                display_format="{:.4f}",
            ),
            "amplitude": SliderConfig(
                slider_name="solventOA_T_slider",
                spinbox_name="solventOA_T_doubleSpinBox",
                param_name="T",
                min_val=-2.0,
                max_val=2.0,
                step=0.01,
                display_format="{:.4f}",
            ),
            "centerpoint": SliderConfig(
                slider_name="solventOA_centerPoint_slider",
                spinbox_name="solventOA_centerPoint_doubleSpinBox",
                param_name="centerpoint",
                min_val=-50,
                max_val=50,
                step=1.0,
                display_format="{:.1f}",
            ),
        },
        # SAMPLE OA
        ("sample", "OA"): {
            "zero_level": SliderConfig(
                slider_name="sampleOA_zeroLevel_slider",
                spinbox_name="sampleOA_zeroLevel_doubleSpinBox",
                param_name="zero_level",
                min_val=0.75,
                max_val=1.25,
                step=0.01,
                display_format="{:.4f}",
            ),
            "amplitude": SliderConfig(
                slider_name="sampleOA_T_slider",
                spinbox_name="sampleOA_T_doubleSpinBox",
                param_name="T",
                min_val=-2.0,
                max_val=2.0,
                step=0.01,
                display_format="{:.4f}",
            ),
            "centerpoint": SliderConfig(
                slider_name="sampleOA_centerPoint_slider",
                spinbox_name="sampleOA_centerPoint_doubleSpinBox",
                param_name="centerpoint",
                min_val=-50,
                max_val=50,
                step=1.0,
                display_format="{:.1f}",
            ),
        },
    }

    def __init__(self, ui_window) -> None:
        """
        Initialize slider controller with reference to UI window

        Args:
            ui_window: Main window with loaded UI widgets
        """
        self.ui = ui_window
        self.slider_to_config: Dict[str, SliderConfig] = {}
        self.config_by_sample_aperture: Dict[Tuple, Dict] = self.SLIDER_CONFIGS

        # Build lookup table
        for (sample, aperture), configs in self.SLIDER_CONFIGS.items():
            for param, config in configs.items():
                self.slider_to_config[config.slider_name] = config

    def connect_sliders(
        self,
        sample_type: SampleType,
        aperture: ApertureType,
        on_slider_change: Callable,
    ) -> None:
        """
        Connect all sliders for a given sample/aperture to callback

        Args:
            sample_type: SampleType.SILICA, .SOLVENT, or .SAMPLE
            aperture: ApertureType.CA or .OA
            on_slider_change: Callback function(sample, aperture, param_name, physical_value, config)
        """
        key = (sample_type.value, aperture.value)
        configs = self.SLIDER_CONFIGS.get(key, {})

        for param_name, config in configs.items():
            slider = getattr(self.ui, config.slider_name, None)
            if slider is None:
                print(f"Warning: Slider {config.slider_name} not found in UI")
                continue

            # Get actual slider range from widget
            slider_min = slider.minimum()
            slider_max = slider.maximum()

            # Create lambda with proper closure - PASS SLIDER WIDGET
            def make_callback(s, a, p, cfg, sldr, s_min, s_max):
                def cb(value):
                    # Convert slider value to physical using ACTUAL slider range
                    physical_value = (
                        SliderController.slider_to_physical_with_range(
                            cfg, value, s_min, s_max
                        )
                    )
                    on_slider_change(s, a, p, physical_value, cfg)

                return cb

            callback = make_callback(
                sample_type,
                aperture,
                param_name,
                config,
                slider,
                slider_min,
                slider_max,
            )
            slider.valueChanged.connect(callback)

            print(
                f"✓ Connected {config.slider_name} ({param_name}) range=[{slider_min}, {slider_max}]"
            )

    @staticmethod
    def slider_to_physical_with_range(
        config: SliderConfig,
        slider_value: int,
        slider_min: int,
        slider_max: int,
    ) -> float:
        """
        Convert slider integer value to physical quantity

        Args:
            config: SliderConfig with min/max physical values
            slider_value: Current slider value
            slider_min: Actual slider minimum from widget
            slider_max: Actual slider maximum from widget
        """
        if slider_max == slider_min:
            return config.min_val

        # Map slider range [slider_min, slider_max] to physical range [min_val, max_val]
        ratio = (slider_value - slider_min) / (slider_max - slider_min)
        physical_value = config.min_val + ratio * (
            config.max_val - config.min_val
        )

        return physical_value

    @staticmethod
    def slider_to_physical(
        config: SliderConfig, slider_value: int, slider_widget=None
    ) -> float:
        """
        Convert slider integer value to physical quantity

        Args:
            config: SliderConfig with min/max physical values
            slider_value: Current slider value (0 to slider_max)
            slider_widget: Optional actual slider widget to get real max value
        """
        # Get actual slider range from widget if provided
        if slider_widget is not None:
            slider_min = slider_widget.minimum()
            slider_max = slider_widget.maximum()
        else:
            # Fallback: assume standard Qt range
            slider_min = 0
            slider_max = 100

        return SliderController.slider_to_physical_with_range(
            config, slider_value, slider_min, slider_max
        )

    @staticmethod
    def physical_to_slider(
        config: SliderConfig, physical_value: float, slider_widget=None
    ) -> int:
        """Convert physical quantity back to slider value"""
        if slider_widget is not None:
            slider_min = slider_widget.minimum()
            slider_max = slider_widget.maximum()
        else:
            slider_min = 0
            slider_max = 100

        # Map physical range to slider range
        ratio = (physical_value - config.min_val) / (
            config.max_val - config.min_val
        )
        ratio = np.clip(ratio, 0, 1)  # Clamp to valid range
        slider_value = int(slider_min + ratio * (slider_max - slider_min))

        return slider_value

    def set_slider_value(
        self,
        sample_type: SampleType,
        aperture: ApertureType,
        param_name: str,
        physical_value: float,
    ):
        """Set slider to a physical value"""
        key = (sample_type.value, aperture.value)
        config = self.SLIDER_CONFIGS.get(key, {}).get(param_name)

        if not config:
            print(f"No config for {key}/{param_name}")
            return

        slider = getattr(self.ui, config.slider_name, None)
        if slider is None:
            return

        slider_value = self.physical_to_slider(config, physical_value)
        slider.blockSignals(True)
        slider.setValue(slider_value)
        slider.blockSignals(False)

        # Update spinbox if exists
        if config.spinbox_name:
            spinbox = getattr(self.ui, config.spinbox_name, None)
            if spinbox:
                spinbox.blockSignals(True)
                spinbox.setValue(physical_value)
                spinbox.blockSignals(False)

    def set_slider_values(
        self,
        sample_type: SampleType,
        aperture: ApertureType,
        param_values: Dict[str, float],
    ) -> None:
        """
        Update multiple sliders at once with a dictionary of param_name: physical_value pairs

        Args:
            sample_type: SampleType.SILICA, SOLVENT, or SAMPLE
            aperture: ApertureType.CA or OA
            param_values: Dict mapping param_name to physical_value
                         e.g., {"amplitude": -0.5, "zero_level": 1.0, "centerpoint": 0.0}

        Example:
            self.slider_controller.set_slider_values(
                sample_type=SampleType("solvent"),
                aperture=ApertureType("CA"),
                param_values={
                    "amplitude": -0.482,
                    "zero_level": 1.001,
                    "centerpoint": 0.5,
                    "beamwaist": 2.55,
                }
            )
        """
        for param_name, physical_value in param_values.items():
            self.set_slider_value(
                sample_type=sample_type,
                aperture=aperture,
                param_name=param_name,
                physical_value=physical_value,
            )

    def update_spinbox_display(
        self, spinbox_name: Optional[str], value: float, format_str: str
    ):
        """Update a spinbox display value"""
        spinbox = getattr(self.ui, spinbox_name or "", None)
        if spinbox:
            spinbox.blockSignals(True)
            spinbox.setValue(value)
            spinbox.blockSignals(False)


if __name__ == "__main__":
    # Test configuration
    print("=" * 70)
    print("SLIDER CONFIGURATION TEST")
    print("=" * 70)

    for (sample, aperture), configs in SliderController.SLIDER_CONFIGS.items():
        print(f"\n{sample.upper()} - {aperture}:")
        for param, config in configs.items():
            print(
                f"  {param:<15} [{config.min_val:>6.2f}, {config.max_val:>6.2f}]"
            )
            print(f"    → {config.slider_name}")
            if config.spinbox_name:
                print(f"    → {config.spinbox_name}")

    print("\n" + "=" * 70)
    print("✓ Configuration loaded successfully")
