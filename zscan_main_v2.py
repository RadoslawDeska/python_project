"""
Main Z-scan application with complete parameter management
- Auto-reads all params from UI spinboxes each fit
- Progress feedback + background threading
- Invalidation warnings only when params actually change from stored values
- Set custom checkbox enable/disable spinboxes
- ROI selection with two vertical lines (click order irrelevant)
- Conditional absorption/saturation visibility
- Sliders enabled per sample type when data loaded
"""

import sys
from pathlib import Path
from typing import Any, Optional, Dict, Tuple
import numpy as np

from PyQt5.QtWidgets import QMainWindow, QApplication, QFileDialog, QMessageBox
from PyQt5.QtCore import pyqtSignal, QObject, QThread
from PyQt5.uic import loadUi
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from zscan_tab_manager_v2 import TabManager

from zscan_data_parser_v2 import (
    RawZScanData,
    ZScanFileParser,
    ZScanProcessor,
    ProcessedZScanData,
)
from zscan_config_manager_v2 import FittingConfig
from zscan_physics_v2 import (
    ClosedAperturePhysics,
    OpenAperturePhysics,
    FittingResult,
)
from zscan_slider_controller_v2 import (
    ApertureType,
    SampleType,
    SliderController,
)


class FitProgressEmitter(QObject):
    """Emits progress signals from fitting thread"""

    progress = pyqtSignal(str)
    finished = pyqtSignal(str, object)  # (sample_aperture, result)
    error = pyqtSignal(str, str)  # (sample_aperture, error_msg)


class FittingThread(QThread):
    """Background thread for running fits without blocking UI"""

    def __init__(self, fit_func, sample_aperture: str, **kwargs):
        super().__init__()
        self.fit_func = fit_func
        self.sample_aperture = sample_aperture
        self.kwargs = kwargs
        self.emitter = FitProgressEmitter()

    def run(self):
        try:
            self.emitter.progress.emit(f"Fitting {self.sample_aperture}...")
            result = self.fit_func(**self.kwargs)
            if result:
                self.emitter.finished.emit(self.sample_aperture, result)
            else:
                self.emitter.error.emit(
                    self.sample_aperture, "Fitting returned None"
                )
        except Exception as e:
            self.emitter.error.emit(self.sample_aperture, str(e))
            import traceback

            traceback.print_exc()


class CanvasROISelector:
    """Manages ROI (Region of Interest) selection with two vertical lines"""

    def __init__(self, ax, on_roi_selected=None):
        self.ax = ax
        self.on_roi_selected = on_roi_selected
        self.roi_limits = None
        self.click_count = 0
        self.line1 = None
        self.line2 = None
        self.cid = None  # Connection ID for click event
        self.enabled = False

        # Crosshair lines
        self.crosshair_v = None
        self.crosshair_h = None
        self.cid_motion = None
        self.cid_leave = None

    def _safe_remove(self, artist):
        """Remove an artist safely across all Matplotlib backends."""
        if artist is None:
            return

        # Try direct removal
        try:
            artist.remove()
            return
        except Exception:
            pass

        # Try removing from axes containers
        try:
            if artist in self.ax.lines:
                self.ax.lines[:] = [l for l in self.ax.lines if l is not artist]
                return
        except Exception:
            pass

        try:
            if artist in self.ax.artists:
                self.ax.artists[:] = [
                    a for a in self.ax.artists if a is not artist
                ]
                return
        except Exception:
            pass

        # Last resort: hide it
        try:
            artist.set_visible(False)
        except Exception:
            print("Lines hidden!!")
            pass

    def enable(self):
        """Enable ROI selection mode"""
        self.enabled = True
        self.click_count = 0
        self.roi_limits = None
        self._remove_lines()
        # Connect mouse click event
        self.cid = self.ax.figure.canvas.mpl_connect(
            "button_press_event", self._on_click
        )

        # Crosshair events
        self.cid_motion = self.ax.figure.canvas.mpl_connect(
            "motion_notify_event", self._on_motion
        )
        self.cid_leave = self.ax.figure.canvas.mpl_connect(
            "axes_leave_event", self._on_leave
        )

    def disable(self):
        """Disable ROI selection mode"""
        self.enabled = False
        self.roi_limits = None

        if self.cid is not None:
            self.ax.figure.canvas.mpl_disconnect(self.cid)
            self.cid = None

        if self.cid_motion is not None:
            self.ax.figure.canvas.mpl_disconnect(self.cid_motion)
            self.cid_motion = None

        if self.cid_leave is not None:
            self.ax.figure.canvas.mpl_disconnect(self.cid_leave)
            self.cid_leave = None

        self._remove_lines()
        self._remove_crosshair()

    def _remove_lines(self):
        """Remove vertical line markers"""
        self._safe_remove(self.line1)
        self._safe_remove(self.line2)
        self.line1 = None
        self.line2 = None
        self.ax.figure.canvas.draw_idle()

    def _remove_crosshair(self):
        self._safe_remove(self.crosshair_v)
        self._safe_remove(self.crosshair_h)
        self.crosshair_v = None
        self.crosshair_h = None
        self.ax.figure.canvas.draw_idle()

    def _on_motion(self, event):
        """Update crosshair on mouse movement"""
        if not self.enabled or event.inaxes != self.ax:
            return

        x, y = event.xdata, event.ydata
        if x is None or y is None:
            return

        # Create crosshair lines if missing
        if self.crosshair_v is None:
            self.crosshair_v = self.ax.axvline(
                x, color="gray", linestyle=":", linewidth=1
            )
        if self.crosshair_h is None:
            self.crosshair_h = self.ax.axhline(
                y, color="gray", linestyle=":", linewidth=1
            )

        # Update positions
        self.crosshair_v.set_xdata([x, x])
        self.crosshair_h.set_ydata([y, y])

        self.ax.figure.canvas.draw_idle()

    def _on_leave(self, event):
        """Hide crosshair when cursor leaves axes"""
        self._remove_crosshair()

    def _on_click(self, event):
        """Handle mouse click on canvas"""
        if not self.enabled or event.inaxes != self.ax or event.xdata is None:
            return

        self.click_count += 1

        if self.click_count == 1:
            # First click - draw first line
            self.line1 = self.ax.axvline(
                event.xdata,
                color="red",
                linestyle="--",
                linewidth=1.5,
                alpha=0.7,
            )
            self.ax.figure.canvas.draw_idle()

        elif self.click_count == 2:
            # Second click - draw second line and store limits
            self.line2 = self.ax.axvline(
                event.xdata,
                color="red",
                linestyle="--",
                linewidth=1.5,
                alpha=0.7,
            )

            # Store limits (order doesn't matter
            x_min = min(self.line1.get_xdata()[0], event.xdata)
            x_max = max(self.line1.get_xdata()[0], event.xdata)
            self.roi_limits = (x_min, x_max)

            if self.on_roi_selected:
                self.on_roi_selected(self.roi_limits)

            self.ax.figure.canvas.draw_idle()

        elif self.click_count >= 3:
            # Third click - reset selection
            self.click_count = 0
            self.roi_limits = None
            self._remove_lines()

            if self.on_roi_selected:
                self.on_roi_selected(None)

    def get_mask(self, x_data):
        """Get boolean mask for data points within ROI"""
        if self.roi_limits is None:
            return np.ones(len(x_data), dtype=bool)
        x_min, x_max = self.roi_limits
        return (x_data >= x_min) & (x_data <= x_max)


class UIParameterManager:
    """Manages reading and writing parameters from/to UI spinboxes"""

    def __init__(self, window):
        self.window = window

    def get_all_params(self) -> Dict[str, float]:
        """Read ALL parameters from UI spinboxes"""
        return {
            "wavelength_nm": self.window.wavelength_dataFittingTab_doubleSpinBox.value(),
            "zscan_range_mm": self.window.zscanRange_doubleSpinBox.value(),
            "aperture_diameter_mm": self.window.apertureDiameter_doubleSpinBox.value(),
            "d0_mm": self.window.apertureToFocusDistance_doubleSpinBox.value(),
            "silica_thickness_mm": self.window.silicaThickness_dataFittingTab_doubleSpinBox.value(),
            "concentration_pct": self.window.concentration_dataFittingTab_doubleSpinBox.value(),
        }

    def set_params(self, params: Dict[str, float]):
        """Write parameters to UI (blocking signals to avoid loops)"""
        mappings = {
            "wavelength_nm": "wavelength_dataFittingTab_doubleSpinBox",
            "zscan_range_mm": "zscanRange_doubleSpinBox",
            "aperture_diameter_mm": "apertureDiameter_doubleSpinBox",
            "d0_mm": "apertureToFocusDistance_doubleSpinBox",
            "silica_thickness_mm": "silicaThickness_dataFittingTab_doubleSpinBox",
            "concentration_pct": "concentration_dataFittingTab_doubleSpinBox",
        }

        for param_key, widget_name in mappings.items():
            if param_key in params and hasattr(self.window, widget_name):
                widget = getattr(self.window, widget_name)
                widget.blockSignals(True)
                widget.setValue(params[param_key])
                widget.blockSignals(False)

    def enable_spinbox(self, spinbox_name: Optional[str], enable: bool):
        """Enable/disable spinbox and toggle readonly state"""
        if spinbox_name is not None and hasattr(self.window, spinbox_name):
            widget = getattr(self.window, spinbox_name)
            widget.setEnabled(True)  # Allows programmatic update
            widget.setReadOnly(not enable)  # Toggle user input

    def enable_sliders(self, sample_type: str, enable: bool):
        """Enable/disable all sliders for a sample type"""
        slider_names = {
            "silica": [
                "silicaCA_zeroLevel_slider",
                "silicaCA_DPhi0_slider",
                "silicaCA_centerPoint_slider",
                "silicaCA_RayleighLength_slider",
                "silicaCA_filterSize_slider",
            ],
            "solvent": [
                "solventCA_zeroLevel_slider",
                "solventCA_DPhi0_slider",
                "solventCA_centerPoint_slider",
                "solventCA_RayleighLength_slider",
                "solventCA_filterSize_slider",
                "solventOA_zeroLevel_slider",
                "solventOA_T_slider",
                "solventOA_centerPoint_slider",
                "solventOA_filterSize_slider",
            ],
            "sample": [
                "sampleCA_zeroLevel_slider",
                "sampleCA_DPhi0_slider",
                "sampleCA_centerPoint_slider",
                "sampleCA_RayleighLength_slider",
                "sampleCA_filterSize_slider",
                "sampleOA_zeroLevel_slider",
                "sampleOA_T_slider",
                "sampleOA_centerPoint_slider",
                "sampleOA_filterSize_slider",
            ],
        }

        for slider_name in slider_names.get(sample_type, []):
            if hasattr(self.window, slider_name):
                slider = getattr(self.window, slider_name)
                slider.setEnabled(enable)


class ZScanMainWindow(QMainWindow):
    """Main application window with full parameter management"""

    def __init__(self, ui_path: str = "window.ui"):
        """Modified __init__ to add curve cache"""
        super().__init__()

        # Load UI
        try:
            self.ui = loadUi(ui_path, self)
            print(f"✓ Loaded UI from {ui_path}")
        except Exception as e:
            print(f"✗ Failed to load UI: {e}")
            sys.exit(1)

        # Create mapping of fitting tabs to switch easily
        self.tab_mapping = {
            self.fittingTabs.tabText(i).lower(): i
            for i in range(self.fittingTabs.count())
        }

        # Initialize managers
        self.param_manager = UIParameterManager(self)
        self.tab_manager = TabManager(self.fittingTabs)
        self.config = FittingConfig()

        # Data storage
        self.data: dict[str, ProcessedZScanData] = {}
        self.file_data: dict[str, RawZScanData] = {}
        self.fit_results: dict[tuple[str, str], FittingResult] = {}
        self.roi_selectors: dict[str, CanvasROISelector] = {}
        self.canvases: dict[str, dict[str, Any]] = {}
        self.fitted_params_snapshot: dict[tuple[str, str], Any] = {}
        self.fit_thread: Optional[FittingThread] = None

        self.fitted_curves_original: dict[tuple[str, str], np.ndarray] = {}
        self.fitted_curves_extended: dict[tuple[str, str], np.ndarray] = {}

        # Initialize slider controller
        self.slider_controller = SliderController(self)
        print("✓ Slider controller initialized")

        # Setup
        self._setup_canvases()
        self._setup_custom_checkboxes()
        self._setup_roi_selection()
        self._setup_absorption_initial_state()
        self._setup_absorption_visibility()
        self._connect_file_buttons()
        self._connect_fit_buttons()
        self._connect_absorption_changes()

        self.statusbar.showMessage("Ready - Load data to begin")

    def _setup_canvases(self):
        """Create matplotlib canvases for all CA/OA tabs"""
        configs = [
            ("silica", "CA", "silicaCA_layout"),
            ("silica", "OA", "silicaOA_layout"),
            ("solvent", "CA", "solventCA_layout"),
            ("solvent", "OA", "solventOA_layout"),
            ("sample", "CA", "sampleCA_layout"),
            ("sample", "OA", "sampleOA_layout"),
        ]

        for sample_type, aperture, layout_name in configs:
            if hasattr(self, layout_name):
                fig = Figure(figsize=(5, 4), dpi=100)
                ax = fig.add_subplot(111)
                canvas = FigureCanvas(fig)

                getattr(self, layout_name).addWidget(canvas)

                key = f"{sample_type}_{aperture}"
                self.canvases[key] = {"fig": fig, "ax": ax, "canvas": canvas}
                print(f"✓ Canvas: {key}")

    def _setup_custom_checkboxes(self):
        """Connect 'Set custom' checkboxes to enable/disable spinboxes"""
        custom_map = {
            "customWavelength_checkBox": (
                "wavelength_dataFittingTab_doubleSpinBox",
            ),
            "customZscanRange_checkBox": ("zscanRange_doubleSpinBox",),
            "customApertureDiameter_checkBox": (
                "apertureDiameter_doubleSpinBox",
            ),
            "customApertureToFocusDistance_checkBox": (
                "apertureToFocusDistance_doubleSpinBox",
            ),
            "customSilicaThickness_checkBox": (
                "silicaThickness_dataFittingTab_doubleSpinBox",
            ),
            "customConcentration_checkBox": (
                "concentration_dataFittingTab_doubleSpinBox",
            ),
        }

        for checkbox_name, (spinbox_name,) in custom_map.items():
            if hasattr(self, checkbox_name):
                checkbox = getattr(self, checkbox_name)
                checkbox.stateChanged.connect(
                    lambda state,
                    sn=spinbox_name: self._on_custom_checkbox_changed(
                        sn, state == 2
                    )
                )
                print(f"✓ Connected {checkbox_name}")

    def _on_custom_checkbox_changed(
        self, spinbox_name: Optional[str], enabled: bool
    ):
        """Handle custom checkbox - enable/disable spinbox editing without blocking value updates

        Args:
            spinbox_name: Name of the spinbox to control
            enabled: True = user can edit, False = readonly (but values still update)
        """
        if spinbox_name is not None and hasattr(self, spinbox_name):
            widget = getattr(self, spinbox_name)

            # CRITICAL FIX: Use setReadOnly instead of setEnabled
            # This allows values to be updated programmatically even when user can't edit
            widget.setReadOnly(not enabled)

            # Keep widget visible and enabled so updates work
            widget.setEnabled(True)

            status = "custom" if enabled else "locked"
            print(f"✓ {spinbox_name} set to {status} (readonly={not enabled})")

    def _check_invalidated_fits_on_param_change(
        self, old_params: Dict, new_params: Dict
    ):
        """Check if any fits are invalidated by parameter changes"""
        param_map = {
            "wavelength_nm": "wavelength_nm",
            "zscan_range_mm": "zscan_range_mm",
            "aperture_diameter_mm": "aperture_diameter_mm",
            "d0_mm": "d0_mm",
            "silica_thickness_mm": "silica_thickness_mm",
            "concentration_pct": "concentration_pct",
        }

        invalidated = set()

        for param_key, config_key in param_map.items():
            # Check if parameter actually changed
            if param_key in old_params and param_key in new_params:
                if abs(old_params[param_key] - new_params[param_key]) > 1e-6:
                    # Parameter changed - check what's invalid
                    newly_invalid = self.config.invalidate_dependents(
                        config_key
                    )
                    invalidated.update(newly_invalid)

        if invalidated:
            self.config.clear_fitted(invalidated)
            samples = ", ".join(sorted(invalidated))
            self.statusbar.showMessage(
                f"⚠ Parameter changed - please refit: {samples}"
            )

            QMessageBox.warning(
                self,
                "Refit Required",
                f"Parameter changed. The following fits are now invalid:\n\n{samples}\n\n"
                f"Please run Fit again to update the results.",
                QMessageBox.Ok,
            )

    def _setup_roi_selection(self):
        """Setup ROI selection checkboxes for CA tabs"""
        ca_configs = [
            ("silica", "CA", "silicaCA_fixROI_checkBox"),
            ("solvent", "CA", "solventCA_fixROI_checkBox"),
            ("sample", "CA", "sampleCA_fixROI_checkBox"),
        ]

        for sample_type, aperture, checkbox_name in ca_configs:
            if hasattr(self, checkbox_name):
                checkbox = getattr(self, checkbox_name)
                key = f"{sample_type}_{aperture}"
                checkbox.stateChanged.connect(
                    lambda state, k=key: self._toggle_roi(k, state == 2)
                )
                print(f"✓ ROI checkbox connected: {checkbox_name}")

    def _init_roi_selector(self, canvas_key: str):
        """Initialize ROI selector for a canvas after it's created"""
        if canvas_key not in self.canvases:
            return

        if canvas_key in self.roi_selectors:
            return  # Already initialized

        ax = self.canvases[canvas_key]["ax"]
        roi = CanvasROISelector(
            ax, lambda lim, k=canvas_key: self._on_roi_selected(k, lim)
        )
        self.roi_selectors[canvas_key] = roi
        print(f"✓ ROI selector initialized: {canvas_key}")

    def _toggle_roi(self, canvas_key: str, enabled: bool):
        """Enable/disable ROI selection on a canvas"""
        self._init_roi_selector(canvas_key)

        if canvas_key in self.roi_selectors:
            roi_selector = self.roi_selectors[canvas_key]

            if enabled:
                # Enable ROI selection mode
                roi_selector.enable()
                self.statusbar.showMessage(
                    "ROI selection ON - click twice to set limits (3rd click resets)"
                )
            else:
                # CRITICAL FIX: Clear ROI limits when checkbox is unchecked
                roi_selector.roi_limits = None
                roi_selector.disable()

                # Redraw plot without ROI highlighting (full range)
                sample_type = canvas_key.split("_")[0]
                aperture = canvas_key.split("_")[1]

                data_key = {
                    "silica": "silica_ca",
                    "solvent": "solvent_ca",
                    "sample": "sample_ca",
                }.get(sample_type)

                if data_key in self.data:
                    result = self.fit_results.get((sample_type, aperture))
                    self.plot_data(
                        sample_type, aperture, self.data[data_key], result
                    )
                    print(
                        f"✓ ROI cleared for {canvas_key} - showing full range"
                    )

                self.statusbar.showMessage(
                    "ROI selection OFF - full range active"
                )

    def _on_roi_selected(
        self, canvas_key: str, limits: Optional[Tuple[float, float]]
    ):
        """Callback when ROI is selected"""
        if limits:
            self.statusbar.showMessage(
                f"ROI: {limits[0]:.2f} to {limits[1]:.2f} mm (3rd click to reset)"
            )
        else:
            self.statusbar.showMessage("ROI selection reset")

    def _setup_absorption_initial_state(self):
        """Setup initial state for absorption checkboxes - CORRECTED LOGIC

        "Assume no absorption" checkbox:
        - CHECKED (True): No absorption mode - HIDE absorption/saturation controls
        - UNCHECKED (False): Absorption exists - SHOW absorption/saturation controls
        """
        # Solvent OA: Default to NO absorption (checked = True)
        if hasattr(self, "solventOA_isAbsorption_checkBox"):
            self.solventOA_isAbsorption_checkBox.setChecked(True)
            self.solventOA_isAbsorption_checkBox.setEnabled(True)
            print("✓ Solvent OA: Assume NO absorption (checked)")

        # Sample OA: Default to WITH absorption (unchecked = False)
        if hasattr(self, "sampleOA_isAbsorption_checkBox"):
            self.sampleOA_isAbsorption_checkBox.setChecked(False)
            self.sampleOA_isAbsorption_checkBox.setEnabled(True)
            print("✓ Sample OA: Assume absorption PRESENT (unchecked)")

    def _setup_absorption_visibility(self):
        """Setup conditional visibility for absorption and saturation models

        Logic:
        - isAbsorption checkbox CHECKED = "Assume no absorption" = HIDE models
        - isAbsorption checkbox UNCHECKED = "Assume absorption" = SHOW models
        """
        for sample_type in ["solvent", "sample"]:
            checkbox_name = f"{sample_type}OA_isAbsorption_checkBox"

            if hasattr(self, checkbox_name):
                checkbox = getattr(self, checkbox_name)
                checkbox.stateChanged.connect(
                    lambda state, st=sample_type: self._update_absorption_ui(
                        st,
                        state == 2,  # state==2 means CHECKED
                    )
                )

                # Set initial visibility based on checkbox state
                # If checked (True), hide absorption controls
                self._update_absorption_ui(sample_type, checkbox.isChecked())
                print(f"✓ Absorption visibility setup: {sample_type}")

    def _update_absorption_ui(self, sample_type: str, hide_absorption: bool):
        """Update visibility of absorption and saturation model controls

        Args:
            sample_type: "solvent" or "sample"
            hide_absorption: True when "Assume no absorption" is CHECKED
                            False when "Assume no absorption" is UNCHECKED
        """
        abs_label = f"{sample_type}OA_absorptionModel_label"
        abs_combo = f"{sample_type}OA_absorptionModel_comboBox"
        sat_label = f"{sample_type}OA_saturationModel_label"
        sat_combo = f"{sample_type}OA_saturationModel_comboBox"

        # Show/hide absorption model controls
        # Hide when checkbox is checked (no absorption mode)
        # Show when checkbox is unchecked (absorption present)
        for name in [abs_label, abs_combo]:
            if hasattr(self, name):
                getattr(self, name).setVisible(not hide_absorption)

        # Show saturation model only if:
        # 1. Absorption controls are visible (not hidden)
        # 2. Saturation model is selected in the absorption combo
        show_saturation = False
        if not hide_absorption and hasattr(self, abs_combo):
            text = getattr(self, abs_combo).currentText()
            show_saturation = "SA" in text or "RSA" in text

        for name in [sat_label, sat_combo]:
            if hasattr(self, name):
                getattr(self, name).setVisible(show_saturation)

        status = "hidden" if hide_absorption else "visible"
        print(f"✓ {sample_type} OA absorption controls: {status}")

    def _connect_absorption_changes(self):
        """Connect absorption model combo changes to update saturation visibility"""
        for sample_type in ["solvent", "sample"]:
            combo_name = f"{sample_type}OA_absorptionModel_comboBox"
            if hasattr(self, combo_name):
                combo = getattr(self, combo_name)
                combo.currentTextChanged.connect(
                    lambda text,
                    st=sample_type: self._on_absorption_model_changed(st)
                )
                print(f"✓ Absorption combo connected: {combo_name}")

    def _on_absorption_model_changed(self, sample_type: str):
        """Update saturation visibility when absorption model changes"""
        checkbox_name = f"{sample_type}OA_isAbsorption_checkBox"
        if hasattr(self, checkbox_name):
            # Get checkbox state: True = "Assume no absorption" (hide)
            is_checked = getattr(self, checkbox_name).isChecked()
            self._update_absorption_ui(sample_type, is_checked)

    def _connect_file_buttons(self):
        """Connect data file loading buttons"""
        buttons = {
            "customSilicaFile_pushButton": ("silica", "silica_ca"),
            "customSolventFile_pushButton": ("solvent", "solvent_ca"),
            "customSampleFile_pushButton": ("sample", "sample_ca"),
        }

        for btn_name, (st, dk) in buttons.items():
            if hasattr(self, btn_name):
                getattr(self, btn_name).clicked.connect(
                    lambda checked=False, s=st, d=dk: self.load_file(s, d)
                )
                print(f"✓ File button connected: {btn_name}")

    def _connect_fit_buttons(self):
        """Connect fit buttons"""
        buttons = [
            ("silicaCA_fit_pushButton", "silica", "CA"),
            ("silicaOA_fit_pushButton", "silica", "OA"),
            ("solventCA_fit_pushButton", "solvent", "CA"),
            ("solventOA_fit_pushButton", "solvent", "OA"),
            ("sampleCA_fit_pushButton", "sample", "CA"),
            ("sampleOA_fit_pushButton", "sample", "OA"),
        ]

        for btn_name, st, ap in buttons:
            if hasattr(self, btn_name):
                getattr(self, btn_name).clicked.connect(
                    lambda checked=False, s=st, a=ap: self.on_fit_clicked(s, a)
                )
                print(f"✓ Fit button connected: {btn_name}")

    def load_file(self, sample_type: str, data_key: str):
        """Load Z-scan data file"""
        filepath, _ = QFileDialog.getOpenFileName(
            self,
            f"Load {sample_type.title()} Data",
            "",
            "Text Files (*.txt);;CSV Files (*.csv)",
        )

        if not filepath:
            return

        raw = ZScanFileParser.parse(filepath)
        if not raw:
            QMessageBox.critical(
                self,
                "Parse Error",
                f"Cannot parse {Path(filepath).name}\n\n"
                f"This file may be corrupted or in unsupported format.",
            )
            return

        # Check for missing header information
        warnings = []
        if raw.wavelength_nm == 800.0:
            warnings.append(
                "Wavelength not found in file - using default 800 nm"
            )

        if warnings:
            msg = "File loaded with warnings:\n\n" + "\n".join(
                f"• {w}" for w in warnings
            )
            msg += "\n\nPlease verify parameters in the 'Data fitting' tab"
            QMessageBox.warning(self, "Load Warnings", msg)

            # Auto-enable custom checkboxes for missing params
            if "Wavelength not found" in "\n".join(warnings):
                if hasattr(self, "customWavelength_checkBox"):
                    self.customWavelength_checkBox.setChecked(True)

        # Process data
        processed = ZScanProcessor.normalize(raw)
        self.data[data_key] = processed
        self.file_data[sample_type] = raw

        # Update filename display
        fn_widget = f"{sample_type}Filename_lineEdit"
        if hasattr(self, fn_widget):
            getattr(self, fn_widget).setText(Path(filepath).name)

        # Update UI parameters from file
        self.param_manager.set_params(
            {
                "wavelength_nm": raw.wavelength_nm,
                "zscan_range_mm": processed.z_range_mm,
            }
        )

        # Enable sliders for this sample type
        self.param_manager.enable_sliders(sample_type, True)

        # Switch to currently loaded sample type tab
        self.tab_manager.switch_to(sample_type, "CA")

        # Enable fit buttons
        for btn in [
            f"{sample_type}CA_fit_pushButton",
            f"{sample_type}OA_fit_pushButton",
        ]:
            if hasattr(self, btn):
                getattr(self, btn).setEnabled(True)

        # Plot initial raw data (no fit yet)
        for aperture in ["CA", "OA"]:
            self.plot_data(sample_type, aperture, processed, None)

        # Connect sliders for this sample type
        sample_enum = SampleType(sample_type)
        self.slider_controller.connect_sliders(
            sample_enum, ApertureType.CA, self.on_slider_changed
        )

        self.statusbar.showMessage(f"✓ Loaded {raw.sample_code}")

    def on_slider_changed(
        self, sample_type, aperture, param_name, physical_value, config
    ):
        """Callback when slider changes - modify curve exactly like minimal test"""
        import numpy as np
        
        # Get data and fit
        data_map = {
            "silica": "silica_ca",
            "solvent": "solvent_ca",
            "sample": "sample_ca",
        }
        data_key = data_map.get(sample_type.value)
        data = self.data.get(data_key)
        
        if data is None:
            return
        
        sample_val = sample_type.value
        aperture_val = aperture.value
        
        fit_result = self.fit_results.get((sample_val, aperture_val))
        if fit_result is None:
            return
        
        # Update spinbox display
        if config.spinbox_name:
            spinbox = getattr(self, config.spinbox_name, None)
            if spinbox:
                spinbox.blockSignals(True)
                spinbox.setValue(physical_value)
                spinbox.blockSignals(False)
        
        # Status message
        status_msg = (
            f"{sample_val.upper()} {aperture_val} - {param_name}: {physical_value:.4f}"
        )
        self.statusbar.showMessage(status_msg)
        
        # Get canvas
        key = f"{sample_val}_{aperture_val}"
        if key not in self.canvases:
            return
        
        canvas_dict = self.canvases[key]
        ax = canvas_dict["ax"]
        fig = canvas_dict["fig"]
        canvas_widget = canvas_dict["canvas"]
        
        x = data.position_centered
        y_raw = data.ca_antisym if aperture_val == "CA" else data.oa
                
        if param_name == "zero_level":
            # Shift curve: y_new = y_current + (new_zero - old_zero)
            old_zero = fit_result.params.zero_level
            shift = physical_value - old_zero
            fit_result.y_fit = fit_result.y_fit + shift
            fit_result.params.zero_level = physical_value
        
        elif param_name == "DPhi0":
            # Scale amplitude: y_new = zero + (y_current - zero) * (new_amp / old_amp)
            old_amp = fit_result.params.amplitude
            old_zero = fit_result.params.zero_level
            
            if old_amp != 0:
                scale = physical_value / old_amp
                fit_result.y_fit = old_zero + (fit_result.y_fit - old_zero) * scale
            
            fit_result.params.amplitude = physical_value
        
        elif param_name == "z0":
            # Scale peak width
            wavelength = data.wavelength_nm * 1e-9
            z0_m = physical_value * 1e-3
            new_beamwaist = np.sqrt(z0_m * wavelength / np.pi)
            
            old_beamwaist = fit_result.params.beamwaist
            old_zero = fit_result.params.zero_level
            
            if old_beamwaist > 0 and new_beamwaist > 0:
                scale = old_beamwaist / new_beamwaist
                fit_result.y_fit = old_zero + (fit_result.y_fit - old_zero) * scale
            
            fit_result.params.beamwaist = new_beamwaist
        
        elif param_name == "centerpoint":
            # Use pre-calculated extended curve and shift to new center
            old_center = fit_result.params.centerpoint
            center_shift_points = int(physical_value - old_center)
            
            if center_shift_points != 0:
                # Get the pre-calculated extended curve (calculated once during fit)
                y_extended = self.fitted_curves_extended.get((sample_val, aperture_val))
                
                if y_extended is not None:
                    # Roll the extended curve to the new center position
                    y_fit_shifted = np.roll(y_extended, center_shift_points)
                    # Trim to original length
                    fit_result.y_fit = y_fit_shifted[:len(fit_result.y_fit)]
                else:
                    # Fallback: just roll current curve if extended not available
                    fit_result.y_fit = np.roll(fit_result.y_fit, center_shift_points)
            
            fit_result.params.centerpoint = physical_value
        
        # Clear and redraw (exactly like minimal test)
        ax.clear()
        
        # Check for ROI
        roi_selector = self.roi_selectors.get(key)
        if roi_selector and roi_selector.roi_limits:
            mask = roi_selector.get_mask(x)
            ax.scatter(
                x[~mask],
                y_raw[~mask],
                color="gray",
                s=15,
                alpha=0.3,
                label="Outside ROI",
            )
            ax.scatter(
                x[mask],
                y_raw[mask],
                color="#3b82f6",
                s=30,
                alpha=0.7,
                label="Inside ROI",
            )
        else:
            ax.scatter(x, y_raw, color="#3b82f6", s=30, alpha=0.7, label="Data")
        
        # Plot fit curve
        ax.plot(x, fit_result.y_fit, color="#ef4444", linewidth=2, label="Fit")
        
        # Labels
        ax.set_xlabel("Position (mm)")
        ylabel = "ΔT/T₀" if aperture_val == "CA" else "1 - T"
        ax.set_ylabel(ylabel)
        ax.set_title(f"{sample_val.upper()} {aperture_val} Z-scan")
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        fig.tight_layout()
        canvas_widget.draw()
        
        # Update parameter spinboxes
        prefix = f"{sample_val}{aperture_val}"
        
        if aperture_val == "CA":
            amp_widget = f"{prefix}_deltaPhi0Summary_doubleSpinBox"
        else:
            amp_widget = f"{prefix}_TSummary_doubleSpinBox"
        
        if hasattr(self, amp_widget):
            widget = getattr(self, amp_widget)
            widget.blockSignals(True)
            widget.setValue(fit_result.params.amplitude)
            widget.blockSignals(False)
        
        if aperture_val == "CA":
            bw_widget = f"{prefix}_beamwaistSummary_doubleSpinBox"
            if hasattr(self, bw_widget):
                widget = getattr(self, bw_widget)
                widget.blockSignals(True)
                widget.setValue(fit_result.params.beamwaist * 1e6)
                widget.blockSignals(False)
        
        print(f"✓ {sample_val} {aperture_val} - {param_name}: {physical_value:.4f}")

    def plot_data(
        self,
        sample_type: str,
        aperture: str,
        data: ProcessedZScanData,
        result: Optional[FittingResult],
    ):
        """Plot data with optional fit curve"""
        try:
            key = f"{sample_type}_{aperture}"
            if key not in self.canvases:
                return

            # Initialize ROI selector if needed
            self._init_roi_selector(key)

            ax = self.canvases[key]["ax"]
            fig = self.canvases[key]["fig"]
            canvas = self.canvases[key]["canvas"]

            ax.clear()

            # Select data based on aperture type
            x = data.position_centered
            y = data.ca_antisym if aperture == "CA" else data.oa
            ylabel = "ΔT/T₀" if aperture == "CA" else "1 - T"

            # Check if ROI is selected
            roi_selector = self.roi_selectors.get(key)
            if roi_selector and roi_selector.roi_limits:
                mask = roi_selector.get_mask(x)
                ax.scatter(
                    x[~mask],
                    y[~mask],
                    color="gray",
                    s=15,
                    alpha=0.3,
                    label="Outside ROI",
                )
                ax.scatter(
                    x[mask],
                    y[mask],
                    color="#3b82f6",
                    s=30,
                    alpha=0.7,
                    label="Inside ROI",
                )
            else:
                ax.scatter(x, y, color="#3b82f6", s=30, alpha=0.7, label="Data")

            # Plot fit if available
            if result:
                ax.plot(
                    x, result.y_fit, color="#ef4444", linewidth=2, label="Fit"
                )

            ax.set_xlabel("Position (mm)")
            ax.set_ylabel(ylabel)
            ax.set_title(f"{sample_type.upper()} {aperture} Z-scan")
            ax.grid(True, alpha=0.3)
            ax.legend()

            fig.tight_layout()
            canvas.draw()
            print(f"✓ Plot updated: {key}")

        except Exception as e:
            print(f"! Plot error: {e}")

    def on_fit_clicked(self, sample_type: str, aperture: str):
        """Handle fit button click"""
        data_map = {
            "silica": "silica_ca",
            "solvent": "solvent_ca",
            "sample": "sample_ca",
        }
        data_key = data_map[sample_type]

        data = self.data.get(data_key)
        file_data = self.file_data.get(sample_type)

        if data is None or file_data is None:
            QMessageBox.warning(
                self, "No Data", f"Please load {sample_type} data first"
            )
            return

        # READ ALL PARAMETERS FROM UI
        new_params = self.param_manager.get_all_params()

        # Check if parameters changed since last fit
        old_params = self.fitted_params_snapshot.get((sample_type, aperture))
        if old_params:
            self._check_invalidated_fits_on_param_change(old_params, new_params)

        self.statusbar.showMessage(f"Fitting {sample_type} {aperture}...")

        if aperture == "CA":
            # Infer initial parameters from data
            initial = ClosedAperturePhysics.infer_initial_params(
                data.ca_antisym,
                data.position_centered,
                new_params["wavelength_nm"],
            )

            # Create and start fitting thread
            fit_thread = FittingThread(
                ClosedAperturePhysics.fit_ca,
                f"{sample_type}_{aperture}",
                ca_antisym=data.ca_antisym,
                position_centered_mm=data.position_centered,
                wavelength_nm=new_params["wavelength_nm"],
                params=initial,
            )
        else:
            # OA requires CA reference - CRITICAL VALIDATION
            ca_result = self.fit_results.get((sample_type, "CA"))
            if not ca_result:
                QMessageBox.warning(
                    self,
                    "Need CA Reference",
                    f"Please fit {sample_type} CA first to determine beamwaist.\n\n"
                    f"The closed aperture measurement provides the beam profile needed "
                    f"for open aperture fitting.",
                )
                self.statusbar.showMessage(
                    "OA fit cancelled - need CA reference"
                )
                return

            # Validate beamwaist is reasonable
            if (
                ca_result.params.beamwaist <= 0
                or ca_result.params.beamwaist > 1e-3
            ):
                QMessageBox.warning(
                    self,
                    "Invalid Beamwaist",
                    f"Beamwaist from CA fit is invalid: {ca_result.params.beamwaist * 1e6:.2f} µm\n\n"
                    f"Please refit the CA data with valid parameters.",
                )
                self.statusbar.showMessage(
                    "OA fit cancelled - invalid CA beamwaist"
                )
                return

            fit_thread = FittingThread(
                OpenAperturePhysics.fit_oa,
                f"{sample_type}_{aperture}",
                oa=data.oa,
                position_centered_mm=data.position_centered,
                wavelength_nm=new_params["wavelength_nm"],
                beamwaist_m=ca_result.params.beamwaist,
                beta=0.01,
                zero_level=1.0,
                centerpoint=0.0,
                d0_m=new_params["d0_mm"] * 1e-3,
                ra_m=new_params["aperture_diameter_mm"] * 1e-3 / 2,
            )

        # Connect signals
        self.fit_thread = fit_thread
        fit_thread.emitter.progress.connect(self.statusbar.showMessage)
        fit_thread.emitter.finished.connect(self._on_fit_done)
        fit_thread.emitter.error.connect(self._on_fit_error)
        fit_thread.start()

    def _on_fit_done(self, sample_aperture: str, result: FittingResult):
        """Fit completed successfully - STORE ORIGINAL CURVE"""
        sample_type, aperture = sample_aperture.rsplit("_", 1)
        data_map = {
            "silica": "silica_ca",
            "solvent": "solvent_ca",
            "sample": "sample_ca",
        }

        # Store result and mark as fitted
        self.fit_results[(sample_type, aperture)] = result
        self.config.mark_fitted(sample_type, aperture)
        self.fitted_params_snapshot[(sample_type, aperture)] = (
            self.param_manager.get_all_params()
        )

        # ADD THIS: Cache the original fitted curve
        self.fitted_curves_original[(sample_type, aperture)] = (
            result.y_fit.copy()
        )

        # Update plot
        data = self.data[data_map[sample_type]]
        self.plot_data(sample_type, aperture, data, result)

        # Update results display
        self._update_result_display(sample_type, aperture, result)

        # Update status
        self.statusbar.showMessage(
            f"✓ {sample_aperture} fitted | R²={result.r_squared:.4f}"
        )
        print(f"✓ Fit complete: {sample_aperture}")

    def _update_result_display(
        self, sample_type: str, aperture: str, result: FittingResult
    ):
        """Update UI with fit results"""
        prefix = f"{sample_type}{aperture}"

        # Update amplitude
        if aperture == "CA":
            amp_widget = f"{prefix}_deltaPhi0Summary_doubleSpinBox"
        else:
            amp_widget = f"{prefix}_TSummary_doubleSpinBox"

        if hasattr(self, amp_widget):
            widget = getattr(self, amp_widget)
            widget.blockSignals(True)
            widget.setValue(result.params.amplitude)
            widget.blockSignals(False)

        # Update beamwaist
        bw_widget = f"{prefix}_beamwaistSummary_doubleSpinBox"
        if hasattr(self, bw_widget):
            widget = getattr(self, bw_widget)
            widget.blockSignals(True)
            widget.setValue(result.params.beamwaist * 1e6)
            widget.blockSignals(False)

        # Update zero level
        zl_widget = f"{prefix}_zeroLevel_doubleSpinBox"
        if hasattr(self, zl_widget):
            widget = getattr(self, zl_widget)
            widget.blockSignals(True)
            widget.setValue(result.params.zero_level)
            widget.blockSignals(False)

    def _on_fit_error(self, sample_aperture: str, error_msg: str):
        """Fit failed"""
        self.statusbar.showMessage(f"✗ Fit error: {error_msg}")
        QMessageBox.critical(
            self,
            "Fitting Error",
            f"Error while fitting {sample_aperture}:\n\n{error_msg}",
        )


def main():
    app = QApplication(sys.argv)

    ui_path = Path(__file__).parent / "window.ui"

    if not ui_path.exists():
        print(f"ERROR: window.ui not found at {ui_path}")
        sys.exit(1)

    window = ZScanMainWindow(str(ui_path))
    window.show()

    print("\n" + "=" * 70)
    print("Z-SCAN APPLICATION STARTED")
    print("=" * 70)
    print("Features:")
    print("  ✓ Parameters read from UI each fit")
    print("  ✓ Invalidation warnings only on actual param changes")
    print("  ✓ ROI selection with two vertical lines (3rd click resets)")
    print("  ✓ Sliders enabled per sample type on data load")
    print("  ✓ Absorption model visibility control")
    print("  ✓ Saturation model shows for SA/RSA only")
    print("  ✓ Fit summary displayed in UI")
    print("=" * 70 + "\n")

    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
