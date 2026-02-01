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
from typing import Any, Dict, Optional, Tuple

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtCore import QObject, QThread, pyqtSignal
from PyQt5.QtWidgets import QApplication, QFileDialog, QMainWindow, QMessageBox
from PyQt5.uic import loadUi

from zscan_config_manager_v2 import FittingConfig
from zscan_data_parser_v2 import (
    ProcessedZScanData,
    RawZScanData,
    ZScanFileParser,
    ZScanProcessor,
)
from zscan_physics_v2 import (
    ClosedAperturePhysics,
    FittingParams,
    FittingResult,
    OpenAperturePhysics,
)
from zscan_slider_controller_v2 import (
    ApertureType,
    SampleType,
    SliderConfig,
    SliderController,
)
from zscan_tab_manager_v2 import TabManager


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
                self.ax.lines[:] = [
                    line for line in self.ax.lines if line is not artist
                ]
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
            assert self.line1 is not None
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

    SLIDER_NAMES = {
        "silica": [
            "silicaCA_zeroLevel_slider",
            "silicaCA_DPhi0_slider",
            "silicaCA_centerPoint_slider",
            "silicaCA_Beamwaist_slider",
            "silicaCA_filterSize_slider",
        ],
        "solvent": [
            "solventCA_zeroLevel_slider",
            "solventCA_DPhi0_slider",
            "solventCA_centerPoint_slider",
            "solventCA_Beamwaist_slider",
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
            "sampleCA_Beamwaist_slider",
            "sampleCA_filterSize_slider",
            "sampleOA_zeroLevel_slider",
            "sampleOA_T_slider",
            "sampleOA_centerPoint_slider",
            "sampleOA_filterSize_slider",
        ],
    }

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
        """Enable/disable all sliders for a given sample type"""
        slider_list = UIParameterManager.SLIDER_NAMES.get(sample_type, [])

        for slider_name in slider_list:
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

        # ================================================================
        # Process data - THIS HANDLES REVERSAL INTERNALLY!
        # ================================================================
        processed = ZScanProcessor.normalize(raw)
        if processed is None:
            QMessageBox.critical(self, "Error", "Failed to normalize data")
            return
        
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

        # Connect sliders for this sample type
        sample_enum = SampleType(sample_type)
        self.slider_controller.connect_sliders(
            sample_enum, ApertureType.CA, self.on_slider_changed
        )
        self.slider_controller.connect_sliders(
            sample_enum, ApertureType.OA, self.on_slider_changed
        )

        # Plot initial raw data (no fit yet)
        for aperture in ["CA", "OA"]:
            self.plot_data(sample_type, aperture, processed, None)
        
        # ================================================================
        # Get initial curves for both CA and OA
        # ================================================================

        # READ ALL PARAMETERS FROM UI
        new_params = self.param_manager.get_all_params()

        # Infer initial parameters from data (use CA data for both)
        initial = ClosedAperturePhysics.infer_initial_params(
            processed.ca,
            processed.position_centered,
            new_params["wavelength_nm"],
        )

        # Process each aperture (CA first, then OA)
        for aperture in ["CA", "OA"]:
            result_new = None
            
            try:
                if aperture == "CA":
                    print(f"\n[load_file] Computing initial {aperture} curve...")
                    result_new = ClosedAperturePhysics.fit_ca_manual(
                        ca=processed.ca,
                        position_centered_mm=processed.position_centered,
                        wavelength_nm=processed.wavelength_nm,
                        params=initial,
                    )

                    if result_new is None:
                        print(f"✗ Initial {aperture} fit failed")
                        continue

                else:  # OA
                    # Get CA reference for beamwaist
                    ca_result = self.fit_results.get((sample_type, "CA"))
                    if ca_result is None:
                        print("⚠ Skipping initial OA: need CA reference first")
                        continue

                    print(f"\n[load_file] Computing initial {aperture} curve...")
                    result_new = OpenAperturePhysics.fit_oa_manual(
                        oa=processed.oa,
                        position_centered_mm=processed.position_centered,
                        wavelength_nm=processed.wavelength_nm,
                        beamwaist_m=ca_result.params.beamwaist,
                        beta=initial.amplitude,
                        zero_level=initial.zero_level,
                        centerpoint=initial.centerpoint,
                        d0_m=0.26,
                        ra_m=0.001,
                    )

                    if result_new is None:
                        print(f"✗ Initial {aperture} fit failed")
                        continue

            except Exception as e:
                print(f"✗ Error during initial {aperture} fit: {e}")
                import traceback
                traceback.print_exc()
                continue

            # ================================================================
            # Store and display the result
            # ================================================================
            if result_new is not None:
                # Store the result
                self.fit_results[(sample_type, aperture)] = result_new
                self.config.mark_fitted(sample_type, aperture)
                
                print(f"✓ Initial {aperture} curve computed:")
                print(f"  R²: {result_new.r_squared:.6f}")
                print(f"  Amplitude: {result_new.params.amplitude:+.6f}")
                print(f"  Beamwaist: {result_new.params.beamwaist*1e6:.2f} µm")
                print(f"  Zero level: {result_new.params.zero_level:.6f}")
                print(f"  Centerpoint: {result_new.params.centerpoint:.6f}")
                
                # Update result display spinboxes
                self._update_result_display(sample_type, aperture, result_new)
                
                # Update sliders to match initial fit
                self.slider_controller.set_slider_values(
                    sample_type=SampleType(sample_type),
                    aperture=ApertureType(aperture),
                    param_values={
                        "amplitude": result_new.params.amplitude,
                        "zero_level": result_new.params.zero_level,
                        "centerpoint": result_new.params.centerpoint,
                        "beamwaist": result_new.params.beamwaist,
                    },
                )
                
                # Plot with initial curve
                self.plot_data(sample_type, aperture, processed, result_new)
                
                print(f"✓ Initial {aperture} curve displayed\n")

        self.statusbar.showMessage(f"✓ Loaded {raw.sample_code}")

    def on_slider_changed(
        self,
        sample_type,
        aperture,
        param_name,
        physical_value,
        config: SliderConfig,
    ):
        """
        Handle slider value change - update fit curve in real-time
        """
        sample_val = sample_type.value
        aperture_val = aperture.value

        print(f"\n{'=' * 70}")
        print(f"SLIDER CHANGED: {param_name} → {physical_value}")
        print(f"{'=' * 70}")

        # DEBUGGING AMPLITUDE
        if param_name == "amplitude":
            print(f"\n{'=' * 70}")
            print("SLIDER CHANGED - AMPLITUDE DEBUG")
            print(f"{'=' * 70}")
            print(f"Physical value from slider: {physical_value}")
            print(f"Sign: {'POSITIVE' if physical_value > 0 else 'NEGATIVE'}")
            print(f"Absolute value: {abs(physical_value)}")

            # Check if fit_result has it
            fit_result = self.fit_results.get((sample_type, aperture))
            if fit_result:
                print(
                    f"Fit result amplitude BEFORE: {fit_result.params.amplitude}"
                )
                print(
                    f"Fit result amplitude sign: {'POSITIVE' if fit_result.params.amplitude > 0 else 'NEGATIVE'}"
                )
            print(f"{'=' * 70}\n")

        # ============================================================
        # STEP 1: Get current fit result and data
        # ============================================================
        fit_result = self.fit_results.get((sample_val, aperture_val))
        if fit_result is None:
            print(f"✗ No fit result for {sample_val} {aperture_val}")
            return

        data_map = {
            "silica": "silica_ca",
            "solvent": "solvent_ca",
            "sample": "sample_ca",
        }
        data = self.data.get(data_map[sample_val])
        if data is None:
            print(f"✗ No data for {sample_val}")
            return

        print("Current fit result:")
        print(f"  amplitude: {fit_result.params.amplitude:.6f}")
        print(f"  beamwaist: {fit_result.params.beamwaist * 1e6:.2f} µm")
        print(f"  zero_level: {fit_result.params.zero_level:.4f}")
        print(f"  centerpoint: {fit_result.params.centerpoint:.6f} m")

        # ============================================================
        # STEP 2: Update the parameter in fit_result.params
        # ============================================================
        if config.param_name == "zero_level":
            fit_result.params.zero_level = physical_value
        elif config.param_name == "DPhi0" or config.param_name == "T":
            print(
                f"Step 2: Update the amplitude parameter in fit_result with {physical_value}"
            )
            fit_result.params.amplitude = physical_value
        elif config.param_name == "centerpoint":
            fit_result.params.centerpoint = physical_value
        elif config.param_name == "beamwaist":
            fit_result.params.beamwaist = physical_value
        else:
            print(f"✗ Unknown parameter: {config.param_name}")
            return

        print(f"\nUpdated {config.param_name} to {physical_value}")
        print("New fit_result.params:")
        print(f"  amplitude: {fit_result.params.amplitude:.6f}")
        print(f"  beamwaist: {fit_result.params.beamwaist * 1e6:.2f} µm")
        print(f"  zero_level: {fit_result.params.zero_level:.4f}")
        print(f"  centerpoint: {fit_result.params.centerpoint:.6f} m")

        # ============================================================
        # STEP 3: Call fit function with SAME data, UPDATED params
        # ============================================================
        print("\nRecalculating fit with new parameter...")
        print(f"  Calling fit_{aperture_val}(")
        print(
            f"    data={len(data.ca if aperture_val == 'CA' else data.oa)} points"
        )
        print(
            f"    positions={data.position_centered[0]:.3f} to {data.position_centered[-1]:.3f} mm"
        )
        print(f"    wavelength={data.wavelength_nm} nm")
        print("    params=(updated)")
        print("  )")

        try:
            if aperture_val == "CA":
                result_new = ClosedAperturePhysics.fit_ca_manual(
                    ca=data.ca,  # SAME data
                    position_centered_mm=data.position_centered,  # SAME positions
                    wavelength_nm=data.wavelength_nm,  # SAME wavelength
                    params=fit_result.params,  # UPDATED params
                )

                if result_new is None:
                    print("✗ Fit failed")
                    return

                print("✓ Fit successful")
                print(f"  R²: {result_new.r_squared:.6f}")
                print(f"  χ²: {result_new.chi_squared:.6f}")

                print("ClosedAperturePhysics.fit_ca_manual() done.")
                print(
                    f"Fit results returned params.amplitude: {'POSITIVE' if result_new.params.amplitude > 0 else 'NEGATIVE'}"
                )

            else:  # OA
                # Get CA reference for beamwaist
                ca_result = self.fit_results.get((sample_val, "CA"))
                if ca_result is None:
                    print("✗ Need CA reference for beamwaist")
                    return

                result_new = OpenAperturePhysics.fit_oa_manual(
                    oa=data.oa,  # SAME data
                    position_centered_mm=data.position_centered,  # SAME positions
                    wavelength_nm=data.wavelength_nm,  # SAME wavelength
                    beamwaist_m=ca_result.params.beamwaist,  # From CA fit
                    beta=fit_result.params.amplitude,  # Updated amplitude
                    zero_level=fit_result.params.zero_level,  # Updated
                    centerpoint=fit_result.params.centerpoint,  # Updated
                    d0_m=0.26,  # Config
                    ra_m=0.001,  # Config
                )

                if result_new is None:
                    print("✗ Fit failed")
                    return

                print("✓ Fit successful")
                print(f"  R²: {result_new.r_squared:.6f}")
                print(f"  χ²: {result_new.chi_squared:.6f}")

        except Exception as e:
            print(f"✗ Error during fit: {e}")
            import traceback

            traceback.print_exc()
            return

        # ============================================================
        # STEP 4: Update fit_result with new curve
        # ============================================================
        fit_result.y_fit = result_new.y_fit  # Store the new curve
        fit_result.r_squared = result_new.r_squared
        fit_result.chi_squared = result_new.chi_squared
        if result_new.n2 is not None:
            fit_result.n2 = result_new.n2

        print(f"\n✓ Stored new y_fit curve ({len(fit_result.y_fit)} points)")

        # After updating fit_result.params.amplitude:
        print(f"Updated amplitude to: {fit_result.params.amplitude}")
        print(
            f"Sign preserved: {fit_result.params.amplitude == physical_value}"
        )

        # ============================================================
        # STEP 5: Update UI display values
        # ============================================================
        # Update spinbox if it exists
        if config.spinbox_name:
            spinbox = getattr(self, config.spinbox_name, None)
            if spinbox:
                spinbox.blockSignals(True)
                spinbox.setValue(physical_value)
                spinbox.blockSignals(False)
                print(f"✓ Updated spinbox: {config.spinbox_name}")

        # Update status bar
        self.statusbar.showMessage(
            f"{sample_val.upper()} {aperture_val} - {config.param_name}: {physical_value:.4f}"
        )

        self.capture_slider_curve_data("solvent", "CA")

        # Update canvas
        self.plot_data(sample_val, aperture_val, data, fit_result)

        # ============================================================
        # STEP 7: Update result display spinboxes
        # ============================================================
        try:
            prefix = f"{sample_val}{aperture_val}"

            # Amplitude spinbox
            if aperture_val == "CA":
                amp_widget_name = f"{prefix}_deltaPhi0Summary_doubleSpinBox"
            else:
                amp_widget_name = f"{prefix}_TSummary_doubleSpinBox"

            amp_widget = getattr(self, amp_widget_name, None)
            if amp_widget:
                amp_widget.blockSignals(True)
                amp_widget.setValue(fit_result.params.amplitude)
                amp_widget.blockSignals(False)

            # Beamwaist spinbox
            bw_widget_name = f"{prefix}_beamwaistSummary_doubleSpinBox"
            bw_widget = getattr(self, bw_widget_name, None)
            if bw_widget:
                bw_widget.blockSignals(True)
                bw_widget.setValue(fit_result.params.beamwaist * 1e6)
                bw_widget.blockSignals(False)

            # Zero level spinbox
            zl_widget_name = f"{prefix}_zeroLevel_doubleSpinBox"
            zl_widget = getattr(self, zl_widget_name, None)
            if zl_widget:
                zl_widget.blockSignals(True)
                zl_widget.setValue(fit_result.params.zero_level)
                zl_widget.blockSignals(False)

            print("✓ Updated result display spinboxes")

        except Exception as e:
            print(f"! Error updating spinboxes: {e}")

        print(f"{'=' * 70}")
        print("✓ SLIDER CHANGE COMPLETE")
        print(f"{'=' * 70}\n")

    def capture_slider_curve_data(self, sample_type="solvent", aperture="CA"):
        """
        After moving slider, call this to capture the new curve

        Call: self.capture_slider_curve_data('solvent', 'CA')
        """

        import numpy as np

        print("\n" + "=" * 80)
        print(f"SLIDER CURVE CAPTURE - {sample_type.upper()} {aperture}")
        print("=" * 80)

        data_map = {
            "silica": "silica_ca",
            "solvent": "solvent_ca",
            "sample": "sample_ca",
        }
        data = self.data.get(data_map[sample_type])
        fit_result = self.fit_results.get((sample_type, aperture))

        if not data or not fit_result:
            print("ERROR: No data or fit result")
            return

        print("\nCURRENT SLIDER STATE:")
        print(f"   Amplitude: {fit_result.params.amplitude:+.6f}")

        # Get the y_fit from fit_result
        if fit_result.y_fit is not None:
            y_fit = fit_result.y_fit

            print("\nCURVE SHAPE:")
            print(f"   Min: {np.min(y_fit):.6f}")
            print(f"   Max: {np.max(y_fit):.6f}")

            peak_idx = np.argmax(y_fit)
            valley_idx = np.argmin(y_fit)
            peak_pos = data.position_centered[peak_idx]
            valley_pos = data.position_centered[valley_idx]

            print(
                f"\n   Peak: {y_fit[peak_idx]:.6f} at position {peak_pos:.2f} mm"
            )
            print(
                f"   Valley: {y_fit[valley_idx]:.6f} at position {valley_pos:.2f} mm"
            )

            # Expected for negative amplitude
            if fit_result.params.amplitude < 0:
                print(
                    "\n   EXPECTED (amplitude < 0): Peak before focus, valley after"
                )
                print(
                    f"   ACTUAL: Peak at {peak_pos:.1f}mm, valley at {valley_pos:.1f}mm"
                )
                if peak_pos < 0 < valley_pos:
                    print("   ✓ CORRECT!")
                else:
                    print("   ✗ WRONG!")

            # Expected for positive amplitude
            else:
                print(
                    "\n   EXPECTED (amplitude > 0): Valley before focus, peak after"
                )
                print(
                    f"   ACTUAL: Valley at {valley_pos:.1f}mm, peak at {peak_pos:.1f}mm"
                )
                if valley_pos < 0 < peak_pos:
                    print("   ✓ CORRECT!")
                else:
                    print("   ✗ WRONG!")

        print("\n" + "=" * 80 + "\n")

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
            y = data.ca if aperture == "CA" else data.oa
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

    def _infer_amplitude_sign_from_data(
        self, ca: np.ndarray, position_centered_mm: np.ndarray
    ) -> int:
        """
        Infer amplitude sign from Z-position of peak and valley in centered coordinates.

        Centered coordinates mean: focus (Z=0) is at center of data range.
        Negative Z: before focus
        Positive Z: after focus

        Physics interpretation:
        - Peak before focus (negative Z), valley after focus (positive Z) → NEGATIVE
        - Valley before focus (negative Z), peak after focus (positive Z) → POSITIVE

        Args:
            ca: CA data
            position_centered_mm: Z-positions centered at focus (0 at focus)

        Returns:
            +1 for positive DPhi0, -1 for negative DPhi0
        """
        max_idx = np.argmax(ca)
        min_idx = np.argmin(ca)

        peak_z = position_centered_mm[max_idx]
        valley_z = position_centered_mm[min_idx]

        # Both valleys and peaks should be roughly symmetric around focus
        # Check which comes first along the Z-axis

        if peak_z < 0 and valley_z > 0:
            # Peak at negative Z (before focus), valley at positive Z (after focus)
            # → Self-defocusing (negative)
            return -1
        elif valley_z < 0 and peak_z > 0:
            # Valley at negative Z (before focus), peak at positive Z (after focus)
            # → Self-focusing (positive)
            return +1
        elif peak_z < 0 and valley_z < 0:
            # Both before focus - peak is closer to focus
            # → Self-defocusing (negative)
            return -1 if abs(peak_z) < abs(valley_z) else +1
        elif peak_z > 0 and valley_z > 0:
            # Both after focus - valley is closer to focus
            # → Self-focusing (positive)
            return +1 if abs(valley_z) < abs(peak_z) else -1
        else:
            # Shouldn't happen for well-behaved Z-scan data
            print("⚠ Warning: Unexpected peak/valley configuration")
            return +1  # Default to positive

    def on_fit_clicked(self, sample_type: str, aperture: str):
        """Handle fit button click. Call automatic fit."""
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
                data.ca,
                data.position_centered,
                new_params["wavelength_nm"],
            )

            # Create and start fitting thread
            fit_thread = FittingThread(
                ClosedAperturePhysics.fit_ca_automatic,
                f"{sample_type}_{aperture}",
                ca=data.ca,
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

            absorption_model_combo = f"{sample_type}OA_absorptionModel_comboBox"
            if hasattr(self, absorption_model_combo):
                selected_model = getattr(
                    self, absorption_model_combo
                ).currentText()
            else:
                selected_model = "2PA"

            fit_thread = FittingThread(
                OpenAperturePhysics.fit_oa_automatic,
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
                absorption_model=selected_model,
            )

        # Connect signals
        self.fit_thread = fit_thread
        fit_thread.emitter.progress.connect(self.statusbar.showMessage)
        fit_thread.emitter.finished.connect(self._on_fit_done)
        fit_thread.emitter.error.connect(self._on_fit_error)
        fit_thread.start()

    def _on_fit_done(self, sample_aperture: str, result: FittingResult):
        """Fit completed successfully - handle amplitude sign based on Z-position."""
        sample_type, aperture = sample_aperture.rsplit("_", 1)
        data_map = {
            "silica": "silica_ca",
            "solvent": "solvent_ca",
            "sample": "sample_ca",
        }

        if aperture == "CA":
            # ============================================================
            # NEW: Correct amplitude sign based on Z-position of peak/valley
            # ============================================================
            data = self.data.get(data_map[sample_type])
            if data is not None:
                # Use CENTERED position (focus at 0)
                inferred_sign = self._infer_amplitude_sign_from_data(
                    data.ca, data.position_centered
                )
                amplitude_before = result.params.amplitude

                # Apply sign correction
                if inferred_sign < 0 and result.params.amplitude > 0:
                    result.params.amplitude *= -1
                    print(
                        f"✓ Corrected amplitude sign (neg): "
                        f"{amplitude_before:+.6f} → {result.params.amplitude:+.6f}"
                    )
                elif inferred_sign > 0 and result.params.amplitude < 0:
                    result.params.amplitude *= -1
                    print(
                        f"✓ Corrected amplitude sign (pos): "
                        f"{amplitude_before:+.6f} → {result.params.amplitude:+.6f}"
                    )
                else:
                    print(
                        f"✓ Amplitude sign confirmed: "
                        f"{result.params.amplitude:+.6f} "
                        f"({'positive' if inferred_sign > 0 else 'negative'})"
                    )
            else:
                print("⚠ Could not verify amplitude sign (no data available)")

        # Store result and mark as fitted
        self.fit_results[(sample_type, aperture)] = result
        self.config.mark_fitted(sample_type, aperture)
        self.fitted_params_snapshot[(sample_type, aperture)] = (
            self.param_manager.get_all_params()
        )

        # AFTER storing fit_result, UPDATE SLIDERS AND SPINBOXES
        self._update_all_fit_displays(sample_type, aperture, result)

        # Then sync slider positions
        self._sync_sliders_to_fit(sample_type, aperture, result.params)

        # Update plot
        data = self.data[data_map[sample_type]]
        self.plot_data(sample_type, aperture, data, result)

        # Update results display
        self._update_result_display(sample_type, aperture, result)

        # Update ALL sliders after automatic fit
        param_values = {
            "amplitude": result.params.amplitude,
            "zero_level": result.params.zero_level,
            "centerpoint": result.params.centerpoint,
            "beamwaist": result.params.beamwaist,
        }
        print(param_values)
        self.slider_controller.set_slider_values(
            sample_type=SampleType(sample_type),
            aperture=ApertureType(aperture),
            param_values={
                "amplitude": result.params.amplitude,
                "zero_level": result.params.zero_level,
                "centerpoint": result.params.centerpoint,
                "beamwaist": result.params.beamwaist,
            },
        )

        # Update status
        self.statusbar.showMessage(
            f"✓ {sample_aperture} fitted | R²={result.r_squared:.4f}"
        )
        print(f"✓ Fit complete: {sample_aperture}")

    def _update_all_fit_displays(
        self, sample_type: str, aperture: str, result: FittingResult
    ):
        """Update ALL spinboxes and displays with fit result"""
        prefix = f"{sample_type}{aperture}"

        updates = {
            # Amplitude
            (
                "deltaPhi0Summary_doubleSpinBox"
                if aperture == "CA"
                else "TSummary_doubleSpinBox",
                result.params.amplitude,
            ),
            # Beamwaist
            (
                f"{prefix}_beamwaistSummary_doubleSpinBox",
                result.params.beamwaist * 1e6,
            ),
            # Zero level
            (f"{prefix}_zeroLevel_doubleSpinBox", result.params.zero_level),
        }

        for widget_name, value in updates:
            full_name = f"{prefix}_{widget_name}"
            if hasattr(self, full_name):
                widget = getattr(self, full_name)
                widget.blockSignals(True)
                widget.setValue(value)
                widget.blockSignals(False)

    def _sync_sliders_to_fit(
        self, sample_type: str, aperture: str, params: FittingParams
    ):
        """Sync all sliders to match fitted parameters"""
        self.slider_controller.set_slider_values(
            sample_type=SampleType(sample_type),
            aperture=ApertureType(aperture),
            param_values={
                "amplitude": params.amplitude,
                "zero_level": params.zero_level,
                "centerpoint": params.centerpoint,
                "beamwaist": params.beamwaist,
            },
        )
        print(f"✓ Sliders synced for {sample_type} {aperture}")

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
