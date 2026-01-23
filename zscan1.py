"""Z-Scan Measurement Analysis Application.

Main GUI application for automated Z-scan measurements and analysis.
Integrates modular components for physics, fitting, and data processing.
"""

import json
import logging
import os
import re
import sys
import time
import traceback
from datetime import datetime
from typing import List, Optional, Tuple

import matplotlib
import numpy as np
from numpy.typing import NDArray
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtCore import QObject, QThreadPool, QTimer
from PyQt5.QtGui import QColor, QPalette
from PyQt5.QtWidgets import QFileDialog, QMessageBox, QSlider
from scipy.signal import medfilt

# Import modular components
from lib.constants import (
    BEAMWAIST_UM_SCALE,
    POSITIONS_MM_SCALE,
    CUVETTE_PATH_LENGTH,
    INTEGRATION_STEPS,
    LASER_INTENSITY_GW_SCALE,
    MAX_DPHI0,
    N_COMPONENTS,
    OA_FITTING_PARAMS,
    RAYLEIGH_LENGTH_MM_SCALE,
    SILICA_BETA,
    SOLVENT_CA_SLIDER_CENTER_OFFSET,
    SOLVENT_OA_SLIDER_CENTER_OFFSET,
    SOLVENT_T_SLIDER_MAX,
)
from lib.cursors import BlittedCursor, SnappingCursor
from lib.figure import MplCanvas
from lib.fitting import Fitting
from lib.helpers import as_float, ensure_float
from lib.integration import Integration
from lib.mgmotor import MG17Motor
from lib.sample_fitting import SampleFitter
from lib.scientific_rounding import error_rounding
from lib.slider_mapper import CenterPointSliderMapper, ParameterSliderMapper
from lib.worker import Worker

matplotlib.rcParams.update({"font.size": 7})
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)


class Window(QtWidgets.QMainWindow):
    # INITIALIZATION
    def __init__(self):
        super(Window, self).__init__()
        uic.loadUi(os.path.join(os.path.dirname(__file__), "./window.ui"), self)

        # Additions to UI design
        # self.path = os.path.join("C:/z-scan/_wyniki/") # default main directory for z-scan data
        self.path = os.path.join(os.path.dirname(__file__), "data")
        self.mainDirectory_lineEdit.setText(self.path.replace("/", "\\"))
        self.dataDirectory_lineEdit.setText(self.path.replace("/", "\\"))

        self.solventOA_absorptionModel_label.setVisible(False)
        self.solventOA_absorptionModel_comboBox.setVisible(False)
        self.solventOA_fixROI_checkBox.setVisible(False)

        self.solventOA_saturationModel_label.setVisible(False)
        self.solventOA_saturationModel_comboBox.setVisible(False)

        self.timing_and_threading()
        self.states()
        self.additional_variables()
        self.additional_objects()
        self.value_change_triggers()
        self.slider_triggers()
        self.clicker_triggers()
        self.timer_triggers()

        self.initialize_workflows()
        self.initialize_slider_mappers()

        # SHOW THE APP WINDOW
        self.show()

    def initialize_workflows(self) -> None:
        """Initialize all fitting workflows.

        Workflows are initialized with a window_provider callback.
        This keeps workflows pure and testable (headless) while allowing UI integration.
        """
        from workflows.silica_workflow import SilicaCAWorkflow
        from workflows.solvent_workflow import (
            SolventCAWorkflow,
            SolventOAWorkflow,
        )
        from workflows.sample_workflow import SampleCAWorkflow, SampleOAWorkflow

        # window_provider is a callable that returns the window object
        # This allows workflows to get the window only when needed, but work without it
        window_provider = lambda: self

        self.silica_ca_workflow = SilicaCAWorkflow(
            self.state, window_provider=window_provider
        )

        self.solvent_ca_workflow = SolventCAWorkflow(
            self.state, window_provider=window_provider
        )
        self.solvent_oa_workflow = SolventOAWorkflow(
            self.state, window_provider=window_provider
        )

        self.sample_ca_workflow = SampleCAWorkflow(
            self.state, window_provider=window_provider
        )
        self.sample_oa_workflow = SampleOAWorkflow(
            self.state, window_provider=window_provider
        )

    def initialize_slider_mappers(self) -> None:
        # Initialize slider mappers
        self.t_mapper = ParameterSliderMapper(
            OA_FITTING_PARAMS["T"]["min"],
            OA_FITTING_PARAMS["T"]["max"],
            100,  # slider max
        )
        self.center_mapper = CenterPointSliderMapper(
            slider_max=self.silicaCA_centerPoint_slider.maximum(),
            center_range=(
                self.zscanRange_doubleSpinBox.value() / POSITIONS_MM_SCALE / 2
            ),
        )

    def timing_and_threading(self):
        self.timer = QTimer()
        self.threadpool = QThreadPool()
        print(
            f"Multithreading with maximum {self.threadpool.maxThreadCount()} threads"
        )

    def states(self):
        # State management now delegated to StateManager
        from lib.state_manager import StateManager

        self.state = StateManager()
        # Access states via: self.state.running, self.state.silica_autofit_done, etc.

    def additional_variables(self) -> None:
        """
        Initialize all auxiliary variables used across the application.
        Grouped by purpose for clarity and maintainability.
        """

        # ----------------------------------------------------------------------
        # Data management
        # ----------------------------------------------------------------------
        from lib.data_manager import DataManager

        self.data_mgr: DataManager = DataManager()

        self.silica_data_set: list[NDArray[np.float64]] = []
        self.solvent_data_set: list[NDArray[np.float64]] = []
        self.sample_data_set: list[NDArray[np.float64]] = []

        # ----------------------------------------------------------------------
        # Motor & scan state
        # ----------------------------------------------------------------------
        self.motor_list: List[str] = []
        self.offset: int = 0

        self.previous_end_pos: float = self.endPos_doubleSpinBox.value()
        self.previous_start_pos: float = self.startPos_doubleSpinBox.value()
        self.previous_steps_scan: int = self.stepsScan_spinBox.value()

        # ----------------------------------------------------------------------
        # Header / file correctness
        # ----------------------------------------------------------------------
        self.header_correct: bool = False

        # ----------------------------------------------------------------------
        # Cursor positions for manual fitting
        # ----------------------------------------------------------------------
        self.silicaCA_cursorPositions: List[Tuple[float, float]] = []
        self.solventCA_cursorPositions: List[Tuple[float, float]] = []
        self.solventOA_cursorPositions: List[Tuple[float, float]] = []

        # ----------------------------------------------------------------------
        # Fitting line drawn flags
        # ----------------------------------------------------------------------
        self.silicaCA_fittingLine_drawn: bool = False
        self.solventCA_fittingLine_drawn: bool = False
        self.solventOA_fittingLine_drawn: bool = False
        self.sampleCA_fittingLine_drawn: bool = False
        self.sampleOA_fittingLine_drawn: bool = False

        # ----------------------------------------------------------------------
        # Beam properties
        # ----------------------------------------------------------------------

        self.laserI0: float = 0.0
        self.laserI0Error: float | None = None
        self.laserI0Precision: int | None = None

        self.numericalAperture: float = 0.0
        self.numericalApertureError: float | None = None
        self.numericalAperturePrecision: int | None = None

        # ----------------------------------------------------------------------
        # Silica fitting parameters
        # ----------------------------------------------------------------------
        self.silicaCA_zeroLevel: float = 0.0
        self.silicaCA_zeroLevelError: float | None = None
        self.silicaCA_zeroLevelPrecision: int | None = None

        self.silicaCA_centerPoint: float = 0.0
        self.silicaCA_centerPointError: float | None = None
        self.silicaCA_centerPointPrecision: int | None = None

        self.silicaCA_DPhi0: float = 0.0
        self.silicaCA_DPhi0Error: float | None = None
        self.silicaCA_DPhi0Precision: int | None = None

        self.silicaCA_beamwaist: float = 0.0
        self.silicaCA_beamwaistError: float | None = None
        self.silicaCA_beamwaistPrecision: int | None = None

        self.silica_rayleighLength: float = 0.0
        self.silica_rayleighLengthError: float | None = None
        self.silica_rayleighLengthPrecision: int | None = None

        # ----------------------------------------------------------------------
        # Solvent fitting parameters
        # ----------------------------------------------------------------------
        self.solventCA_zeroLevel: float = 0.0
        self.solventCA_zeroLevelError: float | None = None
        self.solventCA_zeroLevelPrecision: int | None = None

        self.solventCA_centerPoint: float = 0.0
        self.solventCA_centerPointError: float | None = None
        self.solventCA_centerPointPrecision: int | None = None

        self.solventCA_DPhi0: float = 0.0
        self.solventCA_DPhi0Error: float | None = None
        self.solventCA_DPhi0Precision: int | None = None

        self.solventCA_beamwaist: float = 0.0
        self.solventCA_beamwaistError: float | None = None
        self.solventCA_beamwaistPrecision: int | None = None

        self.solvent_n2: float = 0.0
        self.solvent_n2Error: float | None = None
        self.solvent_n2Precision: float | None = None

        self.solvent_rayleighLength: float = 0.0
        self.solvent_rayleighLengthError: float | None = None
        self.solvent_rayleighLengthPrecision: int | None = None

        # ----------------------------------------------------------------------
        # Sample fitting parameters
        # ----------------------------------------------------------------------
        self.sampleCA_DPhi0: float = 0.0
        self.sampleCA_beamwaist: float = 50e-6
        self.sampleCA_zeroLevel: float = 1.0
        self.sampleCA_centerPoint: int = 0

        self.sampleOA_DPhi0: float = 0.0
        self.sampleOA_beamwaist: float = 50e-6
        self.sampleOA_zeroLevel: float = 1.0
        self.sampleOA_centerPoint: int = 0

        # ----------------------------------------------------------------------
        # Autofit tracking flags
        # ----------------------------------------------------------------------
        self.silica_autofit_done: bool = False
        self.solventCA_autofit_done: bool = False
        self.solventOA_autofit_done: bool = False
        self.sampleCA_autofit_done: bool = False
        self.sampleOA_autofit_done: bool = False

    def additional_objects(self):
        # UI manager for initialization
        from lib.ui_manager import UIManager

        self.ui_mgr = UIManager(self)
        self.ui_mgr.setup_initial_visibility()

        # Motor control
        self.ocx = MG17Motor()
        ocx_layout = self.mg17motor_control_vlayout
        ocx_layout.addWidget(self.ocx)

        # Charts
        self.initialize_measurement_charts()
        self.initialize_fitting_charts()

        # Solvents list
        self.load_solvents()

    def load_solvents(self, caller=""):
        # Populate Solvent combobox with data from file
        if caller == "":
            try:
                json_file = open(
                    os.path.join(os.path.dirname(__file__), "solvents.json")
                )

            except FileNotFoundError:
                self.showdialog(
                    "Warning",
                    "solvents.json not found in the default location. Select the file.",
                )
                path = os.path.abspath(
                    os.path.dirname(__file__)
                )  # this is where solvents.json is expected to be
                file = QFileDialog.getOpenFileName(
                    self, "Open File", path, filter="JSON file (*.json)"
                )
                if file[0] != "":  # if dialog was not cancelled
                    json_file = open(file[0])
                else:
                    json_file = ""

            finally:
                if json_file != "":
                    self.solvents = json.load(json_file)
                    json_file.close()
                    self.solventName_comboBox.addItems(
                        [key for key in self.solvents.keys()]
                    )
                    self.solvent_autocomplete()

        elif caller == "LoadSolvents":
            path = os.path.abspath(
                os.path.dirname(__file__)
            )  # this is where solvents.json is expected to be
            file = QFileDialog.getOpenFileName(
                self, "Open File", path, filter="JSON file (*.json)"
            )
            if file[0] != "":  # if dialog was not cancelled
                with open(file[0]) as json_file:
                    self.solvents = json.load(json_file)
                    json_file.close()
                    self.solventName_comboBox.addItems(
                        [key for key in self.solvents.keys()]
                    )
                    self.solvent_autocomplete()
            else:
                return

    def value_change_triggers(self):
        # Measurement Tab
        self.endPos_doubleSpinBox.editingFinished.connect(
            lambda: self.measurement_plot_rescale("end")
        )
        self.startPos_doubleSpinBox.editingFinished.connect(
            lambda: self.measurement_plot_rescale("start")
        )
        self.stepsScan_spinBox.valueChanged.connect(
            self.measurement_plot_rescale
        )

        # Data saving Tab
        self.concentration_dataSavingTab_doubleSpinBox.editingFinished.connect(
            lambda: self.concentration_dataSavingTab_doubleSpinBox.setText(
                self.concentration_dataSavingTab_doubleSpinBox.text() + " %"
            )
            if "%" not in self.concentration_dataSavingTab_doubleSpinBox.text()
            else self.concentration_dataSavingTab_doubleSpinBox.text()
        )
        self.wavelength_dataSavingTab_doubleSpinBox.editingFinished.connect(
            lambda: self.wavelength_dataSavingTab_doubleSpinBox.setText(
                self.wavelength_dataSavingTab_doubleSpinBox.text() + " nm"
            )
            if "nm" not in self.wavelength_dataSavingTab_doubleSpinBox.text()
            else self.wavelength_dataSavingTab_doubleSpinBox.text()
        )

        # Data fitting Tab
        self.solventName_comboBox.currentIndexChanged.connect(
            self.solvent_autocomplete
        )
        self.zscanRange_doubleSpinBox.editingFinished.connect(
            self.set_new_positions
        )

    def slider_triggers(self):
        """Connect slider value changes to their handlers.

        CRITICAL: Slider changes must trigger BOTH plot updates AND fit summary updates.
        """
        # Silica CA sliders
        self.slider_fit_manually_connect(
            self.silicaCA_RayleighLength_slider, "Connect"
        )
        self.slider_fit_manually_connect(
            self.silicaCA_centerPoint_slider, "Connect"
        )
        self.slider_fit_manually_connect(
            self.silicaCA_zeroLevel_slider, "Connect"
        )
        self.slider_fit_manually_connect(self.silicaCA_DPhi0_slider, "Connect")

        # Add fit summary update trigger for Silica CA
        self.silicaCA_RayleighLength_slider.valueChanged.connect(
            lambda: self.update_fit_summary_from_sliders("Silica", "CA")
        )
        self.silicaCA_centerPoint_slider.valueChanged.connect(
            lambda: self.update_fit_summary_from_sliders("Silica", "CA")
        )
        self.silicaCA_zeroLevel_slider.valueChanged.connect(
            lambda: self.update_fit_summary_from_sliders("Silica", "CA")
        )
        self.silicaCA_DPhi0_slider.valueChanged.connect(
            lambda: self.update_fit_summary_from_sliders("Silica", "CA")
        )

        self.silicaCA_filterSize_slider.valueChanged.connect(
            lambda: self.reduce_noise_in_data(
                self.silica_data_set, ftype="Silica", stype="CA"
            )
        )

        # Solvent CA sliders
        self.solventCA_RayleighLength_slider.valueChanged.connect(
            lambda: self.fit_manually(ftype="Solvent", stype="CA")
        )
        self.solventCA_centerPoint_slider.valueChanged.connect(
            lambda: self.fit_manually(ftype="Solvent", stype="CA")
        )
        self.solventCA_zeroLevel_slider.valueChanged.connect(
            lambda: self.fit_manually(ftype="Solvent", stype="CA")
        )
        self.solventCA_DPhi0_slider.valueChanged.connect(
            lambda: self.fit_manually(ftype="Solvent", stype="CA")
        )
        # Add fit summary update for Solvent CA
        self.solventCA_RayleighLength_slider.valueChanged.connect(
            lambda: self.update_fit_summary_from_sliders("Solvent", "CA")
        )
        self.solventCA_centerPoint_slider.valueChanged.connect(
            lambda: self.update_fit_summary_from_sliders("Solvent", "CA")
        )
        self.solventCA_zeroLevel_slider.valueChanged.connect(
            lambda: self.update_fit_summary_from_sliders("Solvent", "CA")
        )
        self.solventCA_DPhi0_slider.valueChanged.connect(
            lambda: self.update_fit_summary_from_sliders("Solvent", "CA")
        )

        self.solventCA_filterSize_slider.valueChanged.connect(
            lambda: self.reduce_noise_in_data(
                self.solvent_data_set, ftype="Solvent", stype="CA"
            )
        )

        # Solvent OA sliders (already have handlers, add summary update)
        self.solventOA_centerPoint_slider.valueChanged.connect(
            self.solventOA_centerPoint_slider_moved
        )
        self.solventOA_centerPoint_slider.valueChanged.connect(
            lambda: self.update_fit_summary_from_sliders("Solvent", "OA")
        )

        self.solventOA_zeroLevel_slider.valueChanged.connect(
            lambda: self.fit_manually(ftype="Solvent", stype="OA")
        )
        self.solventOA_zeroLevel_slider.valueChanged.connect(
            lambda: self.update_fit_summary_from_sliders("Solvent", "OA")
        )

        self.solventOA_T_slider.valueChanged.connect(
            self.solventOA_T_slider_moved
        )
        self.solventOA_T_slider.valueChanged.connect(
            lambda: self.update_fit_summary_from_sliders("Solvent", "OA")
        )

        try:
            self.solventOA_T_doubleSpinBox.valueChanged.connect(
                self.solventOA_T_spin_changed
            )
            self.solventOA_T_doubleSpinBox.valueChanged.connect(
                lambda: self.update_fit_summary_from_sliders("Solvent", "OA")
            )
        except Exception:
            pass

        try:
            self.solventOA_centerPoint_slider.valueChanged.disconnect()
        except Exception:
            pass
        self.solventOA_centerPoint_slider.valueChanged.connect(
            self.solventOA_centerPoint_slider_moved
        )
        self.solventOA_centerPoint_slider.valueChanged.connect(
            lambda: self.update_fit_summary_from_sliders("Solvent", "OA")
        )

        try:
            self.solventOA_centerPoint_doubleSpinBox.valueChanged.connect(
                self.solventOA_centerPoint_spin_changed
            )
            self.solventOA_centerPoint_doubleSpinBox.valueChanged.connect(
                lambda: self.update_fit_summary_from_sliders("Solvent", "OA")
            )
        except Exception:
            pass

        self.solventOA_filterSize_slider.valueChanged.connect(
            lambda: self.reduce_noise_in_data(
                self.solvent_data_set, ftype="Solvent", stype="OA"
            )
        )

        # Sample CA sliders
        self.sampleCA_centerPoint_slider.valueChanged.connect(
            lambda: self.fit_manually(ftype="Sample", stype="CA")
        )
        self.sampleCA_centerPoint_slider.valueChanged.connect(
            lambda: self.update_fit_summary_from_sliders("Sample", "CA")
        )

        self.sampleCA_zeroLevel_slider.valueChanged.connect(
            lambda: self.fit_manually(ftype="Sample", stype="CA")
        )
        self.sampleCA_zeroLevel_slider.valueChanged.connect(
            lambda: self.update_fit_summary_from_sliders("Sample", "CA")
        )

        self.sampleCA_deltaTpv_slider.valueChanged.connect(
            lambda: self.fit_manually(ftype="Sample", stype="CA")
        )
        self.sampleCA_deltaTpv_slider.valueChanged.connect(
            lambda: self.update_fit_summary_from_sliders("Sample", "CA")
        )

        self.sampleCA_filterSize_slider.valueChanged.connect(
            lambda: self.reduce_noise_in_data(
                self.sample_data_set, ftype="Sample", stype="CA"
            )
        )

        # Sample OA sliders
        self.sampleOA_centerPoint_slider.valueChanged.connect(
            lambda: self.fit_manually(ftype="Sample", stype="OA")
        )
        self.sampleOA_centerPoint_slider.valueChanged.connect(
            lambda: self.update_fit_summary_from_sliders("Sample", "OA")
        )

        self.sampleOA_zeroLevel_slider.valueChanged.connect(
            lambda: self.fit_manually(ftype="Sample", stype="OA")
        )
        self.sampleOA_zeroLevel_slider.valueChanged.connect(
            lambda: self.update_fit_summary_from_sliders("Sample", "OA")
        )

        self.sampleOA_filterSize_slider.valueChanged.connect(
            lambda: self.reduce_noise_in_data(
                self.sample_data_set, ftype="Sample", stype="OA"
            )
        )

    def update_fit_summary_from_sliders(self, ftype: str, stype: str) -> None:
        """Update fit summary display when sliders are moved.

        This is called on every slider change to keep fit summary in sync
        with the current curve shown on screen.
        """
        try:
            # Read current slider positions and extract parameters
            self.get_curve_interpretation(
                ftype, stype, "from_geometry", on_data_load=False
            )
            # Update display only (no errors since this is manual fitting)
            self.set_fit_summary(ftype, stype, caller="manual")
        except Exception as e:
            logger.debug(
                f"Could not update fit summary for {ftype} {stype}: {e}"
            )

    def solventOA_T_slider_moved(self):
        """Handler: map slider -> physical T, update spinbox, then trigger manual OA fit."""
        try:
            phys_t = self.t_mapper.slider_to_physical(
                self.solventOA_T_slider.value()
            )
            self.solventOA_T = phys_t

            # Update spinbox
            try:
                self.solventOA_T_doubleSpinBox.blockSignals(True)
                self.solventOA_T_doubleSpinBox.setValue(phys_t)
            finally:
                try:
                    self.solventOA_T_doubleSpinBox.blockSignals(False)
                except Exception:
                    pass

            # Trigger manual fit
            self.fit_manually(ftype="Solvent", stype="OA")

        except Exception:
            self.fit_manually(ftype="Solvent", stype="OA")

    def solventOA_T_spin_changed(self):
        """Handler: map spinbox -> slider position."""
        try:
            t_min = OA_FITTING_PARAMS["T"]["min"]
            t_max = OA_FITTING_PARAMS["T"]["max"]
            val = self.solventOA_T_doubleSpinBox.value()
            if t_max == t_min:
                mapped = 0
            else:
                mapped = (
                    (val - t_min)
                    / (t_max - t_min)
                    * self.solventOA_T_slider.maximum()
                )
            # Update slider without emitting to avoid double-fit
            try:
                self.solventOA_T_slider.blockSignals(True)
                self.solventOA_T_slider.setValue(int(round(mapped)))
            finally:
                try:
                    self.solventOA_T_slider.blockSignals(False)
                except Exception:
                    pass
            self.solventOA_T = val
        except Exception:
            pass

    def solventOA_centerPoint_slider_moved(self):
        """Handler: map centerpoint slider -> physical center value, update spinbox and trigger manual fit."""
        try:
            slider_val = self.solventOA_centerPoint_slider.value()
            # slider maps 0..100 -> centerpoint -50..50 (see values in lib/constants.py)
            center_val = slider_val - SOLVENT_OA_SLIDER_CENTER_OFFSET
            self.solventOA_centerPoint = center_val
            try:
                self.solventOA_centerPoint_doubleSpinBox.blockSignals(True)
                self.solventOA_centerPoint_doubleSpinBox.setValue(center_val)
            finally:
                try:
                    self.solventOA_centerPoint_doubleSpinBox.blockSignals(False)
                except Exception:
                    pass
            # Trigger manual fit for OA (use activated_by None)
            self.fit_manually(ftype="Solvent", stype="OA")
        except Exception:
            self.fit_manually(ftype="Solvent", stype="OA")

    def solventOA_centerPoint_spin_changed(self):
        """Handler: map centerpoint spinbox -> slider position."""
        try:
            val = float(self.solventOA_centerPoint_doubleSpinBox.value())
            mapped = int(round(val + SOLVENT_OA_SLIDER_CENTER_OFFSET))
            try:
                self.solventOA_centerPoint_slider.blockSignals(True)
                self.solventOA_centerPoint_slider.setValue(mapped)
            finally:
                try:
                    self.solventOA_centerPoint_slider.blockSignals(False)
                except Exception:
                    pass
            self.solventOA_centerPoint = val
        except Exception:
            pass

            # Sample
        self.sampleCA_centerPoint_slider.valueChanged.connect(
            lambda: self.fit_manually(ftype="Sample", stype="CA")
        )
        self.sampleCA_zeroLevel_slider.valueChanged.connect(
            lambda: self.fit_manually(ftype="Sample", stype="CA")
        )
        self.sampleCA_deltaTpv_slider.valueChanged.connect(
            lambda: self.fit_manually(ftype="Sample", stype="CA")
        )
        self.sampleCA_filterSize_slider.valueChanged.connect(
            lambda: self.reduce_noise_in_data(
                self.sample_data_set, ftype="Sample", stype="CA"
            )
        )

        self.sampleOA_centerPoint_slider.valueChanged.connect(
            lambda: self.fit_manually(ftype="Sample", stype="OA")
        )
        self.sampleOA_zeroLevel_slider.valueChanged.connect(
            lambda: self.fit_manually(ftype="Sample", stype="OA")
        )
        self.sampleOA_filterSize_slider.valueChanged.connect(
            lambda: self.reduce_noise_in_data(
                self.sample_data_set, ftype="Sample", stype="OA"
            )
        )

    def clicker_triggers(self):
        # Menu triggers
        self.actionLoadSolvents.triggered.connect(
            lambda: self.load_solvents(caller="LoadSolvents")
        )
        self.actionLight.triggered.connect(self.changeSkinLight)
        self.actionDark.triggered.connect(self.changeSkinDark)

        # Top bar
        self.exitRed_pushButton.clicked.connect(self.stop_experiment)

        # Measurement Tab
        self.clear_pushButton.clicked.connect(self.measurement_clear)
        self.focusAt_comboBox.currentTextChanged.connect(
            lambda: self.measurement_plot_rescale(
                self.focusAt_comboBox.currentText()
            )
        )
        self.initialize_pushButton.clicked.connect(self.initialize)
        self.run_pushButton.clicked.connect(self.set_to_start)

        # Save data Tab
        self.chooseDirectory_pushButton.clicked.connect(
            lambda: self.choose_dir(caller="DataSaving")
        )
        self.saveData_pushButton.clicked.connect(self.data_save)
        self.sendToFit_pushButton.clicked.connect(
            lambda: self.data_loader(
                caller="Current Measurement",
                ftype=self.cuvetteType_comboBox.currentText(),
            )
        )
        self.update_pushButton.clicked.connect(self.update_data_and_filenames)

        # Data fitting Tab
        # Files
        self.customDirectory_pushButton.clicked.connect(
            lambda: self.choose_dir(caller="DataFitting")
        )
        self.customSilicaFile_pushButton.clicked.connect(
            lambda: self.data_loader(caller="Load From File", ftype="Silica")
        )
        self.customSolventFile_pushButton.clicked.connect(
            lambda: self.data_loader(caller="Load From File", ftype="Solvent")
        )
        self.customSampleFile_pushButton.clicked.connect(
            lambda: self.data_loader(caller="Load From File", ftype="Sample")
        )
        # General parameters
        self.customApertureDiameter_checkBox.stateChanged.connect(
            lambda: self.enable_custom("ApertureDiameter")
        )
        self.customApertureToFocusDistance_checkBox.stateChanged.connect(
            lambda: self.enable_custom("ApertureDistance")
        )
        self.customSilicaThickness_checkBox.stateChanged.connect(
            lambda: self.enable_custom("SilicaThickness")
        )
        self.customWavelength_checkBox.stateChanged.connect(
            lambda: self.enable_custom("Wavelength")
        )
        self.customZscanRange_checkBox.stateChanged.connect(
            lambda: self.enable_custom("ZscanRange")
        )
        # Sample properties
        self.customConcentration_checkBox.stateChanged.connect(
            lambda: self.enable_custom("Concentration")
        )

        # Silica CA tab
        self.silicaCA_fit_pushButton.clicked.connect(
            lambda: self.fit_automatically(ftype="Silica", stype="CA")
        )
        self.silicaCA_fixROI_checkBox.stateChanged.connect(
            lambda: self.enable_cursors(ftype="Silica", stype="CA")
        )

        # Solvent CA tab
        self.solventCA_fit_pushButton.clicked.connect(
            lambda: self.fit_automatically(ftype="Solvent", stype="CA")
        )
        self.solventCA_fixROI_checkBox.stateChanged.connect(
            lambda: self.enable_cursors(ftype="Solvent", stype="CA")
        )
        self.solventCA_customBeamwaist_checkBox.stateChanged.connect(
            lambda: self.enable_custom("SolventBeamwaist")
        )

        # Solvent OA tab
        self.solventOA_fit_pushButton.clicked.connect(
            lambda: self.fit_automatically(ftype="Solvent", stype="OA")
        )
        self.solventOA_fixROI_checkBox.stateChanged.connect(
            lambda: self.enable_cursors(ftype="Solvent", stype="OA")
        )
        self.solventOA_isAbsorption_checkBox.stateChanged.connect(
            lambda: self.toggle_absorption_model(ftype="Solvent")
        )
        self.solventOA_absorptionModel_comboBox.currentIndexChanged.connect(
            lambda: self.toggle_saturation_model(ftype="Solvent")
        )
        self.solventOA_customCenterPoint_checkBox.stateChanged.connect(
            lambda: self.enable_custom("SolventCenterPoint")
        )

        # Sample CA tab
        self.sampleCA_fit_pushButton.clicked.connect(
            lambda: self.fit_automatically(ftype="Sample", stype="CA")
        )
        self.sampleCA_fixROI_checkBox.stateChanged.connect(
            lambda: self.enable_cursors(ftype="Sample", stype="CA")
        )
        self.sampleCA_customBeamwaist_checkBox.stateChanged.connect(
            lambda: self.enable_custom("SampleBeamwaist")
        )

        # Sample OA tab
        self.sampleOA_fit_pushButton.clicked.connect(
            lambda: self.fit_automatically(ftype="Sample", stype="OA")
        )
        self.sampleOA_fixROI_checkBox.stateChanged.connect(
            lambda: self.enable_cursors(ftype="Sample", stype="OA")
        )
        self.sampleOA_isAbsorption_checkBox.stateChanged.connect(
            lambda: self.toggle_absorption_model(ftype="Sample")
        )
        self.sampleOA_absorptionModel_comboBox.currentIndexChanged.connect(
            lambda: self.toggle_saturation_model(ftype="Sample")
        )
        self.sampleOA_customCenterPoint_checkBox.stateChanged.connect(
            lambda: self.enable_custom("SampleCenterPoint")
        )

    def timer_triggers(self):  # (started with Initalize button click)
        self.timer.timeout.connect(self.motion_detection)

        self.start_timer()

    def initialize(self, *args, **kwargs):
        self.state.initializing = True

        # Initialize detectors
        device_name = "/Dev1"
        self.detector_core_name = device_name + "/ai"
        self.number_of_channels_used = 3

        print("Detectors initialized")

        # Initialize data dictionaries and apply empty data to lines
        for (
            type,
            chart,
        ) in (
            self.charts.items()
        ):  # e.g.: take the tuple ("relative", "self.rel_chart")
            for chan_no in range(self.number_of_channels_used):
                self.data[type].update(
                    {chan_no: []}
                )  # then fill "relative" dictionary in "self.data" dictionary with pairs of channel number and list of values.
                (line,) = chart.axes.plot(
                    self.data["positions"], self.data[type][chan_no], marker="."
                )  # add empty line for each channel on the "relative" chart
                self.measurement_lines[type].update(
                    {chan_no: line}
                )  # update the "lines" dictionary with line for each channel

        # Initialize translation stage motor
        self.motor_detection_and_homing()

        self.state.initialized = True
        self.state.initializing = False

    # INITIALIZE CHARTS
    def initialize_measurement_charts(self):
        """Initializes charts in the 'Measurement' Tab
        and sets their scales using 'measurement_plot_rescale()' method."""
        # "Measurement charts"
        self.rel_chart = MplCanvas(self)
        layout_rel = self.relative_layout
        layout_rel.addWidget(self.rel_chart)

        self.abs_chart = MplCanvas(self)
        self.abs_chart.axes.set_ylabel("Amplitude (V)")
        layout_abs = self.absolute_layout
        layout_abs.addWidget(self.abs_chart)

        self.rms_text = self.abs_chart.axes.text(
            0.87,
            0.9,
            f"RMS noise = {self.data_mgr.rms_value * 100:.3f}%",
            transform=self.abs_chart.axes.transAxes,
            bbox=dict(boxstyle="round", facecolor="white", alpha=1),
        )

        self.charts = {"relative": self.rel_chart, "absolute": self.abs_chart}
        # initialize empty lines and data dictionaries
        self.measurement_lines = {"relative": {}, "absolute": {}}

        # Use DataManager for data storage
        self.data = self.data_mgr.data

        self.measurement_plot_rescale()

        print("Canvas loaded")

    def initialize_fitting_charts(self):
        # Silica chart
        self.silicaCA_figure = MplCanvas(self)
        self.silicaCA_figure.axes.plot(
            [], [], marker="o", ms=5, linestyle="", color="tab:blue"
        )
        self.silicaCA_figure.axes.set_title("Silica - closed aperture")
        self.silicaCA_figure.axes.set_ylabel("Norm. trasmittance")
        self.silicaCA_figure.axes.set_position([0.135, 0.105, 0.825, 0.825])
        layout_ca_silica = self.silicaCA_layout
        layout_ca_silica.addWidget(self.silicaCA_figure)

        self.silicaOA_figure = MplCanvas(self)
        self.silicaOA_figure.axes.plot(
            [], [], marker="o", ms=5, linestyle="", color="tab:blue"
        )
        self.silicaOA_figure.axes.set_title("Silica - open aperture")
        self.silicaOA_figure.axes.set_ylabel("Norm. trasmittance")
        self.silicaOA_figure.axes.set_position([0.135, 0.105, 0.825, 0.825])
        layout_oa_silica = self.silicaOA_layout
        layout_oa_silica.addWidget(self.silicaOA_figure)

        # Solvent charts
        self.solventCA_figure = MplCanvas(self)
        self.solventCA_figure.axes.plot(
            [], [], marker="o", ms=5, linestyle="", color="tab:blue"
        )
        self.solventCA_figure.axes.set_title("Solvent - closed aperture")
        self.solventCA_figure.axes.set_ylabel("Norm. trasmittance")
        self.solventCA_figure.axes.set_position([0.135, 0.105, 0.825, 0.825])
        layout_ca_solvent = self.solventCA_layout
        layout_ca_solvent.addWidget(self.solventCA_figure)

        self.solventOA_figure = MplCanvas(self)
        self.solventOA_figure.axes.plot(
            [], [], marker="o", ms=5, linestyle="", color="tab:blue"
        )
        self.solventOA_figure.axes.set_title("Solvent - open aperture")
        self.solventOA_figure.axes.set_ylabel("Norm. trasmittance")
        self.solventOA_figure.axes.set_position([0.135, 0.105, 0.825, 0.825])
        layout_oa_solvent = self.solventOA_layout
        layout_oa_solvent.addWidget(self.solventOA_figure)

        # Sample charts
        self.sampleCA_figure = MplCanvas(self)
        self.sampleCA_figure.axes.plot(
            [], [], marker="o", ms=5, linestyle="", color="tab:blue"
        )
        self.sampleCA_figure.axes.set_title("Sample - closed aperture")
        self.sampleCA_figure.axes.set_ylabel("Norm. trasmittance")
        self.sampleCA_figure.axes.set_position(
            [0.135, 0.105, 0.825, 0.825]
        )  # [left, bottom, width, height]
        layout_ca_sample = self.sampleCA_layout
        layout_ca_sample.addWidget(self.sampleCA_figure)

        self.sampleOA_figure = MplCanvas(self)
        self.sampleOA_figure.axes.plot(
            [], [], marker="o", ms=5, linestyle="", color="tab:blue"
        )
        self.sampleOA_figure.axes.set_title("Sample - open aperture")
        self.sampleOA_figure.axes.set_ylabel("Norm. trasmittance")
        self.sampleOA_figure.axes.set_position([0.135, 0.105, 0.825, 0.825])
        layout_oa_sample = self.sampleOA_layout
        layout_oa_sample.addWidget(self.sampleOA_figure)

        self.fitting_charts = {
            "Silica": {"CA": self.silicaCA_figure, "OA": self.silicaOA_figure},
            "Solvent": {
                "CA": self.solventCA_figure,
                "OA": self.solventOA_figure,
            },
            "Sample": {"CA": self.sampleCA_figure, "OA": self.sampleOA_figure},
        }

    # MOTOR NAVIGATION
    def motion_detection(self):
        if self.state.initializing:
            self.clearLED_pushButton.setEnabled(False)
            self.initLED_pushButton.setEnabled(True)
            self.initialize_pushButton.setEnabled(False)
            self.runLED_pushButton.setEnabled(False)
        elif self.state.running:
            self.clearLED_pushButton.setEnabled(False)
            self.initLED_pushButton.setEnabled(False)
            self.runLED_pushButton.setEnabled(True)
        elif self.state.clearing:
            self.clearLED_pushButton.setEnabled(True)
            self.initLED_pushButton.setEnabled(False)
            self.runLED_pushButton.setEnabled(False)
        else:
            self.clearLED_pushButton.setEnabled(False)
            self.initLED_pushButton.setEnabled(False)
            self.runLED_pushButton.setEnabled(False)

        if hasattr(self, "motor"):
            if self.motor.is_in_motion:
                self.clear_pushButton.setEnabled(False)
                self.run_pushButton.setEnabled(False)
                self.waitLED_pushButton.setEnabled(True)

                self.stepsScan_spinBox.setEnabled(False)
                self.samplesStep_spinBox.setEnabled(False)

            else:
                if self.state.initializing:
                    self.clear_pushButton.setEnabled(False)
                    self.run_pushButton.setEnabled(False)
                    self.waitLED_pushButton.setEnabled(True)

                    self.startPos_doubleSpinBox.setEnabled(False)
                    self.endPos_doubleSpinBox.setEnabled(False)
                    self.stepsScan_spinBox.setEnabled(False)
                    self.samplesStep_spinBox.setEnabled(False)

                elif self.state.running:
                    self.clear_pushButton.setEnabled(False)
                    self.run_pushButton.setEnabled(False)
                    self.waitLED_pushButton.setEnabled(True)

                    self.startPos_doubleSpinBox.setEnabled(False)
                    self.endPos_doubleSpinBox.setEnabled(False)
                    self.stepsScan_spinBox.setEnabled(False)
                    self.samplesStep_spinBox.setEnabled(False)

                elif self.state.clearing:
                    self.clear_pushButton.setEnabled(False)
                    self.run_pushButton.setEnabled(False)
                    self.waitLED_pushButton.setEnabled(True)

                    self.startPos_doubleSpinBox.setEnabled(False)
                    self.endPos_doubleSpinBox.setEnabled(False)
                    self.stepsScan_spinBox.setEnabled(False)
                    self.samplesStep_spinBox.setEnabled(False)

                else:
                    self.clear_pushButton.setEnabled(True)
                    self.run_pushButton.setEnabled(True)
                    self.waitLED_pushButton.setEnabled(False)

                    self.startPos_doubleSpinBox.setEnabled(True)
                    self.endPos_doubleSpinBox.setEnabled(True)
                    self.stepsScan_spinBox.setEnabled(True)
                    self.samplesStep_spinBox.setEnabled(True)

        if not self.state.data_acquisition_complete:
            self.update_pushButton.setEnabled(False)
            self.sendToFit_pushButton.setEnabled(False)
        else:
            self.update_pushButton.setEnabled(True)

    def motor_detection_and_homing(self, *args, **kwargs):
        motor_id = 40180184

        import thorlabs_apt as apt

        try:
            self.motor = apt.Motor(motor_id)
            self.ocx.configure(motor_id)
            print(f"Motor {motor_id} connected")

        except Exception:
            print("Motion controller not found!")
            print(
                "Try closing other programs that may be accessing the device."
            )
            print("Try re-plugging the USB device and then this program.")

        self.mpositioner = MotorPositioner()
        self.motor.backlash_distance = 0
        self.thread_it(self.mpositioner.movehome)

    def set_custom_pos(self):
        self.thread_it(self.mpositioner.movetocustompos)

    def set_to_start(self):
        self.state.running = True
        self.thread_it(self.mpositioner.movetostart)

    # GUI EVENTS TIMING
    def start_timer(self):
        self.timer.start(
            100
        )  # How often various actions or states are triggered or checked for

    # def stop_timer(self):
    #    self.timer.stop()

    def stop_experiment(self):
        if hasattr(self, "motor"):
            self.motor.stop_profiled()
            self.state.experiment_stopped = True

            self.state.initializing = False
            self.state.running = False
            self.state.clearing = False

            self.measurement_clear()

    # DATA ACQUISITION AND DISPLAY
    def measurement_clear(
        self,
    ):  # clears all in the first two tabs (Measurement and Data Saving)
        self.state.clearing = True
        self.state.data_acquisition_complete = False

        self.data["positions"] = []  # reset positions to empty list

        self.data_mgr.rms_value = 0.0
        self.rms_text.set_text(
            f"RMS noise = {self.data_mgr.rms_value * 100:.3f}%"
        )

        for (
            type,
            chart,
        ) in (
            self.charts.items()
        ):  # e.g.: take the tuple ("relative", "self.rel_chart")
            for chan_no in range(self.number_of_channels_used):
                self.data[type][
                    chan_no
                ] = []  # then fill "relative" dictionary in "self.data" dictionary empty list of values per channel

                # UPDATE LINES INSTEAD OF DELETING AND REINSTANTIATING
                self.measurement_lines[type][chan_no].set_xdata(
                    self.data["positions"]
                )  # and set data of lines in "lines" dictionary to empty lists of positions and values
                self.measurement_lines[type][chan_no].set_ydata(
                    self.data[type][chan_no]
                )

            chart.axes.relim()
            chart.axes.autoscale_view()
            chart.draw_idle()

        self.rawLogData_textBrowser.clear()
        self.fullLogData_textBrowser.clear()

        self.saveData_pushButton.setEnabled(False)

        self.state.clearing = False

    def measurement_plot_rescale(self, who_called=""):
        """Rescales plots in the 'Measurement' Tab based on
        values given in the 'Measurement control' panel\n
        If user changes the focus of chart at specific line, rescaling fits the charts
        to display the line in its min-max y-range."""

        match who_called:
            case "start":
                if (
                    self.startPos_doubleSpinBox.value()
                    >= self.endPos_doubleSpinBox.value()
                ):
                    self.startPos_doubleSpinBox.setValue(
                        self.previous_start_pos
                    )
            case "end":
                if (
                    self.endPos_doubleSpinBox.value()
                    <= self.startPos_doubleSpinBox.value()
                ):
                    self.endPos_doubleSpinBox.setValue(self.previous_end_pos)

            case _:  # if the who_called value is any of the values in the list of "Focus At" dropdown OR "ANYTHING ELSE"
                if not self.state.initialized:  # prevent the problem of accessing data for rescale when there is no data created yet.
                    return

                called_by = self.focusAt_comboBox.currentText()
                for type, chart in self.charts.items():
                    match called_by:
                        case "All":
                            chart.axes.relim()
                            chart.axes.autoscale()
                        case "Closed":
                            try:
                                chart.axes.set_ylim(
                                    min(self.data[type][0]),
                                    max(self.data[type][0]),
                                )
                                chart.axes.relim()
                            except ValueError:
                                # print(f'\nValueError: min() arg is an empty sequence.\nThis message was initiated by user clicking on {called_by}.\nNothing has to be done.')
                                return
                        case "Reference":
                            try:
                                chart.axes.set_ylim(
                                    min(self.data[type][1]),
                                    max(self.data[type][1]),
                                )
                                chart.axes.relim()
                            except ValueError:
                                # print(f'\nValueError: min() arg is an empty sequence.\nThis message was initiated by user clicking on {called_by}.\nNothing has to be done.')
                                return
                        case "Open":
                            try:
                                chart.axes.set_ylim(
                                    min(self.data[type][2]),
                                    max(self.data[type][2]),
                                )
                                chart.axes.relim()
                            except ValueError:
                                # print(f'\nValueError: min() arg is an empty sequence.\nThis message was initiated by user clicking on {called_by}.\nNothing has to be done.')
                                return
                        case _:  # HERE GOES "ANYTHING ELSE"
                            return

                    chart.draw_idle()

        self.focalPoint_doubleSpinBox.setValue(
            (
                self.endPos_doubleSpinBox.value()
                + self.startPos_doubleSpinBox.value()
            )
            / 2
        )
        self.zscanRange_measurementTab_doubleSpinBox.setValue(
            np.abs(
                self.endPos_doubleSpinBox.value()
                - self.startPos_doubleSpinBox.value()
            )
        )

        self.offset = (
            self.endPos_doubleSpinBox.value()
            - self.startPos_doubleSpinBox.value()
        ) / self.stepsScan_spinBox.value()
        for chart in self.charts.values():
            chart.axes.set_xlim(
                self.startPos_doubleSpinBox.value() - self.offset,
                self.endPos_doubleSpinBox.value() + self.offset,
            )

            chart.axes.relim()
            chart.draw_idle()

    def create_raw_log_line(self, step):
        if window.where_to_start == "end":
            new_step = np.abs(step - 200)
            line = f"{new_step:4d}{' '}"
        else:
            line = f"{step:4d}{' '}"

        for chan_no in range(self.number_of_channels_used):
            line += f"{2 * ' '}{self.data['absolute'][chan_no][step]:12.8f}{9 * ' '}"
            # line += f"{2*' '}{abs_chart_data[chan_no][step]:12.8f}{9*' '}"

            if chan_no == self.number_of_channels_used - 1:
                line += f"{2 * ' '}{0:12.8f}"

        if self.where_to_start == "end":
            old_lines = self.rawLogData_textBrowser.toPlainText()
            self.rawLogData_textBrowser.clear()
            self.rawLogData_textBrowser.append(line)
            self.rawLogData_textBrowser.append(old_lines)
        else:
            self.rawLogData_textBrowser.append(line)

    # SELECT DIRECTORY
    def choose_dir(self, caller=""):
        if caller == "DataSaving":
            path = self.mainDirectory_lineEdit.text()
        elif caller == "DataFitting":
            path = self.dataDirectory_lineEdit.text()

        p = QFileDialog.getExistingDirectory(self, "Select a directory", path)
        if (
            p != ""
        ):  # This keeps old path in directory QLineEdit lines, if dialog is closed with Cancel
            if caller == "DataSaving":
                self.mainDirectory_lineEdit.setText(p.replace("/", "\\"))
            elif caller == "DataFitting":
                self.dataDirectory_lineEdit.setText(p.replace("/", "\\"))

    # DATA SAVING
    def data_reverse(self):
        if self.where_to_start == "end":  # and self.data_reversed == False:
            for chan_no in range(self.number_of_channels_used):
                self.data["absolute"][chan_no] = np.flip(
                    self.data["absolute"][0], axis=0
                )
            # self.data_reversed = True

    def update_data_and_filenames(self):
        self.fullLogData_textBrowser.clear()

        # DATA PREVIEW
        # Full description header
        # if self.experimentDescription_plainTextEdit.toPlainText() != "": # this is for full log to look nicer
        #    self.experimentDescription_plainTextEdit.appendPlainText("")

        # full log header
        header = (
            "Z-scan Measurement\n"  # line 0
            # f"Sample type: {self.cuvetteType_comboBox.currentText()}\n"
            f"Code: {self.codeOfSample_lineEdit.text()}\n"  # line 1
            f"Silica thickness: {self.silicaThickness_dataSavingTab_doubleSpinBox.text()}\n"  # line 2
            f"Concentration: {self.concentration_dataSavingTab_doubleSpinBox.text()}\n"  # line 3
            f"Wavelength: {self.wavelength_dataSavingTab_doubleSpinBox.text()}\n"  # line 4
            f"{self.experimentDescription_plainTextEdit.toPlainText()}\n"  # line 5
            "--------------------------------------------------------------------------------------------\n\n"  # line 6
            f"Starting pos: {self.startPos_doubleSpinBox.value()}\n"  # line 7
            f"Ending pos: {self.endPos_doubleSpinBox.value()}\n"  # line 8
            "CH1:   Closed aperture\n"  # line 9
            "CH2:   Reference\n"  # line 10
            "CH3:   Open aperture\n"  # line 11
            "CH4:   Empty channel\n\n"  # line 12
            "--------------------------------------------------------------------------------------------\n\n"  # line 13
            "SNo.  [V] Voltage Max        [V] Voltage Max        [V] Voltage Max        [V] Voltage Max\n\n"  # line 14
            "--------------------------------------------------------------------------------------------\n"
        )  # line 15

        raw_log = self.rawLogData_textBrowser.toPlainText()
        self.fullLogData_textBrowser.append(header)
        self.fullLogData_textBrowser.append(raw_log)

        # log data
        if self.state.data_acquisition_complete:
            raw_log_data = np.genfromtxt(raw_log.split("\n"))
            self.data["absolute"] = {
                0: list(raw_log_data[:, 1]),
                1: list(raw_log_data[:, 2]),
                2: list(raw_log_data[:, 3]),
            }
            self.data_set = (
                list(raw_log_data[:, 0]),
                self.data["absolute"][0],
                self.data["absolute"][1],
                self.data["absolute"][2],
            )  # , data[:,3] not using the last column with zeros

            self.saveData_pushButton.setEnabled(True)
        else:
            self.saveData_pushButton.setEnabled(False)

        # FILENAMES
        now = datetime.now()
        self.cur_date = now.strftime("%Y_%m_%d")
        self.cur_time = now.strftime("%H_%M")
        sample_type = self.codeOfSample_lineEdit.text()

        concentration = self.concentration_dataSavingTab_doubleSpinBox.value()
        conc_hyphen = str(f"{concentration:.2f}").replace(".", "-")

        wavelength = self.wavelength_dataSavingTab_doubleSpinBox.value()
        wavel_hyphen = str(f"{wavelength:.1f}").replace(".", "-")

        self.rawLogFilename_lineEdit.setText(
            f"{self.cur_date}__{self.cur_time}__{sample_type}_{conc_hyphen}_{wavel_hyphen}.txt"
        )
        self.fullLogFilename_lineEdit.setText(
            f"{self.cur_date}__{self.cur_time}__{sample_type}_{conc_hyphen}_{wavel_hyphen}_2.txt"
        )

        self.files = (
            self.rawLogFilename_lineEdit.text(),
            self.fullLogFilename_lineEdit.text(),
        )

    def data_save(self):
        self.accurate_path = os.path.join(
            self.mainDirectory_lineEdit.text(), self.cur_date
        )

        try:
            if not os.path.exists(self.accurate_path):
                os.mkdir(self.accurate_path)

            for file in self.files:
                with open(os.path.join(self.accurate_path, file), "w") as f:
                    if self.files.index(file) == 0:
                        f.write(self.rawLogData_textBrowser.toPlainText())
                    else:
                        f.write(self.fullLogData_textBrowser.toPlainText())

            self.saveData_pushButton.setEnabled(False)
            self.sendToFit_pushButton.setEnabled(True)
            self.showdialog("Info", "Files successfully written!")

        except PermissionError:
            self.showdialog(
                "Error",
                "Permission denied!\nCannot write the file in this directory.",
            )

    # DATA FITTING
    def data_loader(self, caller: str, ftype: str) -> None:
        """Load data either from file or from current measurement ('caller' argument).\n
        'ftype' is passed by proper button for loading from file, or by proper option in "Cuvette type" combobox in "Data Saving" tab"""
        logger.debug(
            f"Window.data_loader running with params: {caller=}, {ftype=}"
        )

        if caller == "Current Measurement":
            # Fill in names in QLineEdits
            self.dataDirectory_lineEdit.setText(self.accurate_path + "\\")
            self.mainTabs.setCurrentIndex(2)

            match ftype:
                case "Silica":
                    self.silicaFilename_lineEdit.setText(
                        self.rawLogFilename_lineEdit.text()
                    )
                    self.fittingTabs.setCurrentIndex(0)
                    self.silicaAperture_tabWidget.setCurrentIndex(0)
                case "Solvent":
                    self.solventFilename_lineEdit.setText(
                        self.rawLogFilename_lineEdit.text()
                    )
                    self.fittingTabs.setCurrentIndex(1)
                    self.solventAperture_tabWidget.setCurrentIndex(0)
                case "Sample":
                    self.sampleFilename_lineEdit.setText(
                        self.rawLogFilename_lineEdit.text()
                    )
                    self.fittingTabs.setCurrentIndex(2)
                    self.sampleAperture_tabWidget.setCurrentIndex(0)

            # Get data
            self.data_set = (
                range(self.stepsScan_spinBox.value() + 1),
                self.data["absolute"][0],
                self.data["absolute"][1],
                self.data["absolute"][2],
            )  # , data[:,3] not using the last column with zeros

            # Read parameters
            self.read_header_params(caller="Current Measurement", ftype=ftype)
            self.data_display(self.data_set, ftype)

            self.switch_fitting_to_on_state(ftype)
            if ftype == "Silica":
                self.silica_autofit_done = (
                    False  # this is to start afresh with fitting
                )

            self.fit_manually(ftype, stype="CA")
            self.fit_manually(ftype, stype="OA")

        elif caller == "Load From File":
            # Load data
            def load_data():
                """Fills proper frame in GUI with file information, sets current tab to lead the user,\n
                loads data on screen, toggles to active state the fitting controls and calls for the initial fit\n
                for both CA and OA traces.
                """
                logger.debug("Window.data_loader.load_data: loading")
                # Fill in names in QLineEdits
                match ftype:
                    case "Silica":
                        self.silicaFilename_lineEdit.setText(fname)
                        self.fittingTabs.setCurrentIndex(0)
                        self.silicaAperture_tabWidget.setCurrentIndex(0)
                    case "Solvent":
                        self.solventFilename_lineEdit.setText(fname)
                        self.fittingTabs.setCurrentIndex(1)
                        self.solventAperture_tabWidget.setCurrentIndex(0)
                    case "Sample":
                        self.sampleFilename_lineEdit.setText(fname)
                        self.fittingTabs.setCurrentIndex(2)
                        self.sampleAperture_tabWidget.setCurrentIndex(0)

                # Get data
                data = np.genfromtxt(file_path, skip_header=last_header_line)
                self.data_set = (
                    data[:, 0],
                    data[:, 1],
                    data[:, 2],
                    data[:, 3],
                )  # , data[:,4] not using the last column with zeros

                # Read parameters
                logger.debug(
                    "Window.data_loader.load_data: reading header parameters"
                )
                self.read_header_params(caller="Load From File", ftype=ftype)
                logger.debug("Window.data_loader.load_data: plotting")
                self.data_display(self.data_set, ftype)

                # Initialize sliders and fit summary from geometry
                logger.debug(
                    "Window.data_loader.load_data: initializing sliders from geometry"
                )
                self._initialize_sliders_from_geometry(ftype)

                self.switch_fitting_to_on_state(ftype)
                if ftype == "Silica":
                    self.silica_autofit_done = (
                        False  # this is to start afresh with fitting
                    )

                self.state._file_loading = True
                self.fit_manually(ftype, stype="CA")
                self.fit_manually(ftype, stype="OA")
                self.state._file_loading = False

            file_path, _ = QFileDialog.getOpenFileName(
                self,
                "Select full description file",
                os.path.normpath(self.dataDirectory_lineEdit.text()),
            )
            # User selected a file
            if file_path:  # This keeps old filename in given file type QLineEdit lines, if dialog is closed with Cancel
                logger.debug("Window.data_loader: file_path correct")

                fname = file_path.split("/")[-1]
                self.dataDirectory_lineEdit.setText(
                    "\\".join(file_path.split("/")[:-1])
                )

                # Header check and manipulation
                with open(file_path, "r") as file:
                    # Correct header must include these and a beacon at the end in the form of "SNo" substring
                    required_matches = [
                        "Concentration",
                        "Wavelength",
                        "Starting pos",
                        "Ending pos",
                        "SNo",
                    ]
                    optional_matches = ["Silica thickness"]
                    header_matches = len(required_matches) * [False]
                    self.header = []
                    last_header_line = 0

                    for line_no, l in enumerate(file):
                        for match in required_matches:
                            if re.match(
                                r"\b%s\b" % match, l
                            ):  # lookup whole words
                                header_matches[
                                    required_matches.index(match)
                                ] = True
                                self.header.append(l.strip())

                                if (
                                    match == "SNo"
                                ):  # This is the header end beacon
                                    last_header_line = line_no + 3

                        opt_matched = 0
                        for opt_match in optional_matches:
                            if re.match(
                                r"\b%s\b" % opt_match, l
                            ):  # lookup whole words
                                header_matches.append(True)
                                self.header.append(l.strip())
                                opt_matched += 1

                    if len(self.header) != 0:
                        if (
                            last_header_line != 0
                        ):  # if the header end beacon was found
                            self.header_correct = all(
                                header_matches[: -(1 + opt_matched)]
                            )  # check if required matches are satisfied

                            if not self.header_correct:
                                self.showdialog(
                                    "Warning",
                                    (
                                        "Header corrupted!\n\n"
                                        "Some parameters might have been loaded improperly.\n\nPut measurement parameters manually."
                                    ),
                                )

                                self.customSilicaThickness_checkBox.setChecked(
                                    True
                                )
                                self.customWavelength_checkBox.setChecked(True)
                                self.customZscanRange_checkBox.setChecked(True)

                            load_data()

                        else:
                            self.header_correct = False
                            self.showdialog(
                                "Warning",
                                (
                                    "Header corrupted!\n\n"
                                    "Load a file without the header."
                                ),
                            )

                    else:
                        self.header_correct = False
                        self.showdialog(
                            "Warning",
                            (
                                "Header not found!\n\n"
                                "Set measurement parameters by yourself or load a file with the header."
                            ),
                        )

                        self.customSilicaThickness_checkBox.setChecked(True)
                        self.customWavelength_checkBox.setChecked(True)
                        self.customZscanRange_checkBox.setChecked(True)

                        load_data()

    def _initialize_sliders_from_geometry(self, ftype: str) -> None:
        """Extract geometry from loaded data and set sliders to initial positions.

        This is called only once when data is loaded - it doesn't affect other sample types.
        """
        try:
            logger.debug(
                "Window._initialize_sliders_from_geometry: trying to get_general_parameters"
            )
            self.get_general_parameters()

            match ftype:
                case "Silica":
                    # Extract geometry from silica data
                    dphi0, bw, zr = self.params_from_geometry("Silica", "CA")
                    self.silicaCA_DPhi0 = dphi0
                    self.silicaCA_beamwaist = bw
                    self.silica_rayleighLength = zr

                    logger.debug(f"""Window._initialize_sliders_from_geometry: retrieved from geometry:\n
                                 {self.silicaCA_DPhi0=}\n{self.silicaCA_beamwaist=}\n{self.silica_rayleighLength=}""")

                    # Set slider positions from geometry
                    self.silicaCA_zeroLevel_slider.blockSignals(True)
                    self.silicaCA_zeroLevel_slider.setValue(
                        int(round(1.0 * 100))
                    )  # Start at 1.0
                    self.silicaCA_zeroLevel_slider.blockSignals(False)

                    self.silicaCA_DPhi0_slider.blockSignals(True)
                    self.silicaCA_DPhi0_slider.setValue(
                        int(
                            round(
                                dphi0
                                / MAX_DPHI0
                                * self.silicaCA_DPhi0_slider.maximum()
                            )
                        )
                    )
                    self.silicaCA_DPhi0_slider.blockSignals(False)

                    self.silicaCA_RayleighLength_slider.blockSignals(True)
                    self.silicaCA_RayleighLength_slider.setValue(
                        int(
                            round(
                                zr
                                / (self.z_range / 2)
                                * self.silicaCA_RayleighLength_slider.maximum()
                            )
                        )
                    )
                    self.silicaCA_RayleighLength_slider.blockSignals(False)

                    self.silicaCA_centerPoint_slider.blockSignals(True)
                    self.silicaCA_centerPoint_slider.setValue(
                        int(round(self.silica_nop / 2))
                    )
                    self.silicaCA_centerPoint_slider.blockSignals(False)

                    # Update fit summary with geometry values
                    self.set_fit_summary("Silica", "CA", caller="manual")

                case "Solvent":
                    # Extract geometry from solvent data
                    dphi0, bw, zr = self.params_from_geometry("Solvent", "CA")
                    self.solventCA_DPhi0 = dphi0
                    self.solventCA_beamwaist = bw
                    self.solvent_rayleighLength = zr

                    # Set slider positions
                    self.solventCA_zeroLevel_slider.blockSignals(True)
                    self.solventCA_zeroLevel_slider.setValue(
                        int(round(1.0 * 100))
                    )
                    self.solventCA_zeroLevel_slider.blockSignals(False)

                    self.solventCA_DPhi0_slider.blockSignals(True)
                    self.solventCA_DPhi0_slider.setValue(
                        int(
                            round(
                                dphi0
                                / MAX_DPHI0
                                * self.solventCA_DPhi0_slider.maximum()
                            )
                        )
                    )
                    self.solventCA_DPhi0_slider.blockSignals(False)

                    self.solventCA_RayleighLength_slider.blockSignals(True)
                    self.solventCA_RayleighLength_slider.setValue(
                        int(
                            round(
                                zr
                                / (self.z_range / 2)
                                * self.solventCA_RayleighLength_slider.maximum()
                            )
                        )
                    )
                    self.solventCA_RayleighLength_slider.blockSignals(False)

                    self.solventCA_centerPoint_slider.blockSignals(True)
                    self.solventCA_centerPoint_slider.setValue(
                        int(round(self.solvent_nop / 2))
                    )
                    self.solventCA_centerPoint_slider.blockSignals(False)

                    # Update fit summary
                    self.set_fit_summary("Solvent", "CA", caller="manual")

                case "Sample":
                    # Extract geometry from sample data
                    dphi0, bw, zr = self.params_from_geometry("Sample", "CA")
                    self.sampleCA_DPhi0 = dphi0
                    self.sampleCA_beamwaist = bw

                    # Set slider positions (similar to solvent/silica)
                    self.sampleCA_zeroLevel_slider.blockSignals(True)
                    self.sampleCA_zeroLevel_slider.setValue(
                        int(round(1.0 * 100))
                    )
                    self.sampleCA_zeroLevel_slider.blockSignals(False)

                    self.sampleCA_deltaTpv_slider.blockSignals(True)
                    self.sampleCA_deltaTpv_slider.setValue(
                        int(
                            round(
                                dphi0
                                / MAX_DPHI0
                                * self.sampleCA_deltaTpv_slider.maximum()
                            )
                        )
                    )
                    self.sampleCA_deltaTpv_slider.blockSignals(False)

                    self.sampleCA_centerPoint_slider.blockSignals(True)
                    self.sampleCA_centerPoint_slider.setValue(
                        int(round(self.sample_nop / 2))
                    )
                    self.sampleCA_centerPoint_slider.blockSignals(False)

        except Exception as e:
            logger.warning(
                f"Could not initialize sliders from geometry for {ftype}: {e}"
            )

    def enable_custom(self, o: str):
        """Toggles readOnly parameter on `o` element from GUI. Uses match-case structure with `o` parameter to match\n
        to speed up processing at each call.

        Args:
            o (str): case for parameter to toggle
        """
        match o:
            case "ApertureDiameter":
                if self.customApertureDiameter_checkBox.isChecked():
                    self.apertureDiameter_doubleSpinBox.setReadOnly(False)
                else:
                    self.apertureDiameter_doubleSpinBox.setReadOnly(True)
            case "ApertureDistance":
                if self.customApertureToFocusDistance_checkBox.isChecked():
                    self.apertureToFocusDistance_doubleSpinBox.setReadOnly(
                        False
                    )
                else:
                    self.apertureToFocusDistance_doubleSpinBox.setReadOnly(True)
            case "SilicaThickness":
                if self.customSilicaThickness_checkBox.isChecked():
                    self.silicaThickness_dataFittingTab_doubleSpinBox.setReadOnly(
                        False
                    )
                else:
                    self.silicaThickness_dataFittingTab_doubleSpinBox.setReadOnly(
                        True
                    )
            case "Wavelength":
                if self.customWavelength_checkBox.isChecked():
                    self.wavelength_dataFittingTab_doubleSpinBox.setReadOnly(
                        False
                    )
                else:
                    self.wavelength_dataFittingTab_doubleSpinBox.setReadOnly(
                        True
                    )
            case "ZscanRange":
                if self.customZscanRange_checkBox.isChecked():
                    self.zscanRange_doubleSpinBox.setReadOnly(False)
                else:
                    self.zscanRange_doubleSpinBox.setReadOnly(True)
            case "Concentration":
                if self.customConcentration_checkBox.isChecked():
                    self.concentration_dataFittingTab_doubleSpinBox.setReadOnly(
                        False
                    )
                else:
                    self.concentration_dataFittingTab_doubleSpinBox.setReadOnly(
                        True
                    )
            case "SolventBeamwaist":
                if not self.solventCA_customBeamwaist_checkBox.isChecked():
                    self.solventCA_RayleighLength_slider.setEnabled(False)
                    if hasattr(window, "silicaCA_beamwaist"):
                        self.solventCA_RayleighLength_slider.valueChanged.disconnect()
                        self.solventCA_RayleighLength_slider.setValue(
                            int(
                                round(
                                    self.silicaCA_beamwaist * BEAMWAIST_UM_SCALE
                                )
                            )
                        )
                        self.solventCA_RayleighLength_slider.valueChanged.connect(
                            lambda: self.fit_manually(
                                ftype="Solvent", stype="CA"
                            )
                        )
                        self.solventCA_beamwaistSummary_doubleSpinBox.setValue(
                            self.silicaCA_beamwaist * BEAMWAIST_UM_SCALE
                        )

                else:
                    self.solventCA_RayleighLength_slider.setEnabled(True)
            case "SolventCenterPoint":  # OA case
                if not self.solventOA_customCenterPoint_checkBox.isChecked():
                    # When custom center is disabled, lock OA center to CA center
                    self.solventOA_centerPoint_slider.setEnabled(False)
                    try:
                        # Use CA center value where available
                        ca_center = getattr(self, "solventCA_centerPoint", 0)
                        self.solventOA_centerPoint = ca_center
                        # Update spinbox and slider to reflect CA center
                        try:
                            self.solventOA_centerPoint_doubleSpinBox.blockSignals(
                                True
                            )
                            self.solventOA_centerPoint_doubleSpinBox.setValue(
                                ca_center
                            )
                        finally:
                            try:
                                self.solventOA_centerPoint_doubleSpinBox.blockSignals(
                                    False
                                )
                            except Exception:
                                pass

                        try:
                            mapped = int(
                                round(
                                    ca_center + SOLVENT_OA_SLIDER_CENTER_OFFSET
                                )
                            )
                            self.solventOA_centerPoint_slider.blockSignals(True)
                            self.solventOA_centerPoint_slider.setValue(mapped)
                        finally:
                            try:
                                self.solventOA_centerPoint_slider.blockSignals(
                                    False
                                )
                            except Exception:
                                pass
                    except Exception:
                        pass
                else:
                    self.solventOA_centerPoint_slider.setEnabled(True)

    def slider_fit_manually_connect(self, current_slider: QSlider, mode):
        """This function is intended to connect and disconnect other sliders on demand, to prevent them from updating the fitting line:
        1) "Disconnect" means current slider should be kept and all others should disconnect from their method
        2) "Connect" means reconnect all sliders to their method"""

        silicaCA_sliders = [
            self.silicaCA_RayleighLength_slider,
            self.silicaCA_centerPoint_slider,
            self.silicaCA_zeroLevel_slider,
            self.silicaCA_DPhi0_slider,
        ]
        # silicaOA_sliders = [self.silicaOA_deltaZpv_slider, self.silicaOA_centerPoint_slider, self.silicaOA_zeroLevel_slider, self.silicaOA_deltaTpv_slider]
        # solventCA_sliders = [self.solventCA_RayleighLength_slider, self.solventCA_centerPoint_slider, self.solventCA_zeroLevel_slider, self.solventCA_DPhi0_slider]
        # solventOA_sliders = [self.solventOA_deltaZpv_slider, self.solventOA_centerPoint_slider, self.solventOA_zeroLevel_slider, self.solventOA_deltaTpv_slider]
        # sampleCA_sliders = [self.sampleCA_deltaZpv_slider, self.sampleCA_centerPoint_slider, self.sampleCA_zeroLevel_slider, self.sampleCA_deltaTpv_slider]
        # sampleOA_sliders = [self.sampleOA_deltaZpv_slider, self.sampleOA_centerPoint_slider, self.sampleOA_zeroLevel_slider, self.sampleOA_deltaTpv_slider]

        available_sliders = [
            silicaCA_sliders
        ]  # , silicaOA_sliders, solventCA_sliders, solventOA_sliders, sampleCA_sliders, sampleOA_sliders]
        self.disconnected_sliders = "All"

        match mode:
            case "Disconnect":
                for sliders in available_sliders:
                    try:
                        sliders.remove(current_slider)
                        for s in sliders:
                            s.disconnect()
                            self.disconnected_sliders = "silicaCA"

                    except ValueError:
                        pass  # If not found

            case "Connect":
                match self.disconnected_sliders:
                    case "silicaCA":
                        for s in available_sliders[0]:
                            s.valueChanged.connect(
                                lambda: self.fit_manually(
                                    ftype="Silica",
                                    stype="CA",
                                    activated_by=self.disconnected_sliders,
                                )
                            )
                    case "silicaOA":
                        for s in available_sliders[1]:
                            s.valueChanged.connect(
                                lambda: self.fit_manually(
                                    ftype="Silica",
                                    stype="OA",
                                    activated_by=self.disconnected_sliders,
                                )
                            )
                    case "solventCA":
                        for s in available_sliders[2]:
                            s.valueChanged.connect(
                                lambda: self.fit_manually(
                                    ftype="Solvent",
                                    stype="CA",
                                    activated_by=self.disconnected_sliders,
                                )
                            )
                    case "solventOA":
                        for s in available_sliders[3]:
                            s.valueChanged.connect(
                                lambda: self.fit_manually(
                                    ftype="Solvent",
                                    stype="OA",
                                    activated_by=self.disconnected_sliders,
                                )
                            )
                    case "sampleCA":
                        for s in available_sliders[4]:
                            s.valueChanged.connect(
                                lambda: self.fit_manually(
                                    ftype="Sample",
                                    stype="CA",
                                    activated_by=self.disconnected_sliders,
                                )
                            )
                    case "sampleOA":
                        for s in available_sliders[5]:
                            s.valueChanged.connect(
                                lambda: self.fit_manually(
                                    ftype="Sample",
                                    stype="OA",
                                    activated_by=self.disconnected_sliders,
                                )
                            )
                    case "All":
                        for s in available_sliders[0]:
                            s.valueChanged.connect(
                                lambda: self.fit_manually(
                                    ftype="Silica", stype="CA"
                                )
                            )
                        # for s in available_sliders[1]:
                        #     s.valueChanged.connect(lambda: self.fit_manually(ftype="Silica", stype="OA"))
                        # for s in available_sliders[2]:
                        #     s.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="CA"))
                        # for s in available_sliders[3]:
                        #     s.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="OA"))
                        # for s in available_sliders[4]:
                        #     s.valueChanged.connect(lambda: self.fit_manually(ftype="Sample", stype="CA"))
                        # for s in available_sliders[5]:
                        #     s.valueChanged.connect(lambda: self.fit_manually(ftype="Sample", stype="OA"))

            case None:
                pass

    def toggle_absorption_model(self, ftype):
        match ftype:
            case "Solvent":
                if not self.solventOA_isAbsorption_checkBox.isChecked():
                    self.solventOA_absorptionModel_label.setVisible(True)
                    self.solventOA_absorptionModel_comboBox.setVisible(True)
                    self.solventOA_fixROI_checkBox.setVisible(True)
                else:
                    self.solventOA_absorptionModel_label.setVisible(False)
                    self.solventOA_absorptionModel_comboBox.setVisible(False)
                    self.solventOA_fixROI_checkBox.setVisible(False)
                    # When assuming no absorption, disable T controls and set OA plot to flat 1
                    try:
                        self.solventOA_T_slider.setEnabled(False)
                    except Exception:
                        pass
                    try:
                        self.solventOA_T_doubleSpinBox.setEnabled(False)
                    except Exception:
                        pass
                    try:
                        line = self.solventOA_figure.axes.get_lines()[0]
                        line.set_ydata(np.ones_like(line.get_ydata()))
                        self.solventOA_figure.axes.relim()
                        self.solventOA_figure.axes.autoscale_view()
                        self.solventOA_figure.draw_idle()
                    except Exception:
                        pass

            case "Sample":
                if not self.sampleOA_isAbsorption_checkBox.isChecked():
                    self.sampleOA_absorptionModel_label.setVisible(True)
                    self.sampleOA_absorptionModel_comboBox.setVisible(True)
                    self.sampleOA_fixROI_checkBox.setVisible(True)
                else:
                    self.sampleOA_absorptionModel_label.setVisible(False)
                    self.sampleOA_absorptionModel_comboBox.setVisible(False)
                    self.sampleOA_fixROI_checkBox.setVisible(False)

    def toggle_saturation_model(self, ftype):
        match ftype:
            case "Solvent":
                if not self.solventOA_isAbsorption_checkBox.isChecked():
                    models = ["SA", "2PA+SA"]  # , "RSA"]
                    if (
                        self.solventOA_absorptionModel_comboBox.currentText()
                        in models
                    ):
                        self.solventOA_saturationModel_label.setVisible(True)
                        self.solventOA_saturationModel_comboBox.setVisible(True)
                    else:
                        self.solventOA_saturationModel_label.setVisible(False)
                        self.solventOA_saturationModel_comboBox.setVisible(
                            False
                        )
            case "Sample":
                if not self.sampleOA_isAbsorption_checkBox.isChecked():
                    models = ["SA", "2PA+SA"]  # , "RSA"]
                    if (
                        self.sampleOA_absorptionModel_comboBox.currentText()
                        in models
                    ):
                        self.sampleOA_saturationModel_label.setVisible(True)
                        self.sampleOA_saturationModel_comboBox.setVisible(True)
                    else:
                        self.sampleOA_saturationModel_label.setVisible(False)
                        self.sampleOA_saturationModel_comboBox.setVisible(False)

    def switch_fitting_to_on_state(self, ftype: str):
        match ftype:
            case "Silica":
                self.silicaCA_RayleighLength_slider.setEnabled(True)
                self.silicaCA_centerPoint_slider.setEnabled(True)
                self.silicaCA_zeroLevel_slider.setEnabled(True)
                self.silicaCA_DPhi0_slider.setEnabled(True)
                self.silicaCA_fit_pushButton.setEnabled(True)
                self.silicaCA_filterSize_slider.setEnabled(True)

            case "Solvent":
                if self.solventCA_customBeamwaist_checkBox.isChecked():
                    self.solventCA_RayleighLength_slider.setEnabled(True)
                else:
                    self.solventCA_RayleighLength_slider.setEnabled(False)
                self.solventCA_centerPoint_slider.setEnabled(True)
                self.solventCA_zeroLevel_slider.setEnabled(True)
                self.solventCA_DPhi0_slider.setEnabled(True)
                self.solventCA_fit_pushButton.setEnabled(True)
                self.solventCA_filterSize_slider.setEnabled(True)
                self.solventCA_customBeamwaist_checkBox.setEnabled(True)

                # Initialize OA center slider/spinbox to 0 (center of data) when Solvent data loads
                try:
                    self.solventOA_centerPoint = 0
                    self.solventOA_centerPoint_doubleSpinBox.blockSignals(True)
                    self.solventOA_centerPoint_doubleSpinBox.setValue(0)
                    self.solventOA_centerPoint_doubleSpinBox.blockSignals(False)
                    self.solventOA_centerPoint_slider.blockSignals(True)
                    self.solventOA_centerPoint_slider.setValue(
                        SOLVENT_OA_SLIDER_CENTER_OFFSET
                    )  # slider SOLVENT_OA_SLIDER_CENTER_OFFSET -> center 0
                    self.solventOA_centerPoint_slider.blockSignals(False)
                except Exception:
                    pass

                self.solventOA_centerPoint_slider.setEnabled(True)
                self.solventOA_zeroLevel_slider.setEnabled(True)
                self.solventOA_T_slider.setEnabled(True)
                self.solventOA_centerPoint_doubleSpinBox.setEnabled(True)
                self.solventOA_zeroLevel_doubleSpinBox.setEnabled(True)
                self.solventOA_T_doubleSpinBox.setEnabled(True)
                self.solventOA_fit_pushButton.setEnabled(True)
                self.solventOA_filterSize_slider.setEnabled(True)
                self.solventOA_isAbsorption_checkBox.setEnabled(True)
                self.solventOA_customCenterPoint_checkBox.setEnabled(True)

            case "Sample":
                pass

    def set_fit_summary(self, ftype: str, stype: str, caller=""):
        """Sets proper number of digits in summary display and shows parameters (with errors) after (automatic) fitting.

        Args:
            ftype (str): Sample type ("Silica","Solvent","Sample")
            stype (str): Z-scan mode ("CA", "OA")
            caller (str, optional): Which fit type has called this function. Defaults to "".
        """

        # DISPLAY VALUES (ALREADY ROUNDED)
        match ftype:
            case "Silica":
                # DPhi0
                if self.silicaCA_DPhi0 is not None:
                    self.silicaCA_deltaPhi0Summary_doubleSpinBox.setValue(
                        self.silicaCA_DPhi0
                    )

                # Laser intensity [GW/cm²]
                if self.laserI0 is not None:
                    self.silicaCA_laserIntensitySummary_doubleSpinBox.setValue(
                        self.laserI0 * LASER_INTENSITY_GW_SCALE
                    )

                # Beamwaist [µm]
                if self.silicaCA_beamwaist is not None:
                    self.silicaCA_beamwaistSummary_doubleSpinBox.setValue(
                        self.silicaCA_beamwaist * BEAMWAIST_UM_SCALE
                    )

                # Rayleigh range [mm]
                if self.silicaCA_beamwaist is not None:
                    self.silicaCA_rayleighRangeSummary_doubleSpinBox.setValue(
                        np.pi
                        * self.silicaCA_beamwaist**2
                        / self.lda
                        * RAYLEIGH_LENGTH_MM_SCALE
                    )

                # Numerical aperture
                if self.numericalAperture is not None:
                    self.numericalAperture_doubleSpinBox.setValue(
                        self.numericalAperture
                    )

            case "Solvent":
                if stype == "CA":
                    self.solventCA_deltaPhi0Summary_doubleSpinBox.setValue(
                        self.solventCA_DPhi0
                    )
                    if self.solvent_n2:
                        self.solventCA_n2Summary_doubleSpinBox.setValue(
                            self.solvent_n2 * 1e13 * 1e9
                        )  # display n2 in multiples of 10^-9 cm^2/GW
                        self.solventOA_n2Summary_doubleSpinBox.setValue(
                            self.solvent_n2 * 1e13 * 1e9
                        )  # display n2 in multiples of 10^-9 cm^2/GW

                    if (
                        self.solventCA_customBeamwaist_checkBox.isChecked()
                        and self.solventCA_beamwaist
                    ):
                        self.solventCA_beamwaistSummary_doubleSpinBox.setValue(
                            self.solventCA_beamwaist * BEAMWAIST_UM_SCALE
                        )  # [um] radius in focal point
                    elif self.silicaCA_beamwaist:
                        self.solventCA_beamwaistSummary_doubleSpinBox.setValue(
                            self.silicaCA_beamwaist * BEAMWAIST_UM_SCALE
                        )

                    if self.solvent_rayleighLength:
                        self.solventCA_rayleighRangeSummary_doubleSpinBox.setValue(
                            self.solvent_rayleighLength
                            * RAYLEIGH_LENGTH_MM_SCALE
                        )  # [mm]

                elif stype == "OA":
                    pass
                    # self.solventOA_TSummary_doubleSpinBox.setValue(self.solventOA_T)
                    # TEMPORARILY DISABLED
                    # self.solventOA_betaSummary_doubleSpinBox.setValue(self.solvent_beta) # CUVETTE_PATH_LENGTH is the solvent/sample thickness assuming no one-photon absorption

            case "Sample":
                pass

        if caller == "manual":
            # The errors are unknown until automatic fitting, so set to #.##
            match ftype:
                case "Silica":
                    if not self.silica_autofit_done:
                        self.silicaCA_deltaPhi0ErrorSummary_label.setText(
                            "#.##"
                        )
                        self.silicaCA_laserIntensityErrorSummary_label.setText(
                            "#.##"
                        )
                        self.silicaCA_beamwaistErrorSummary_label.setText(
                            "#.##"
                        )
                        self.silicaCA_rayleighRangeErrorSummary_label.setText(
                            "#.##"
                        )
                        self.numericalApertureErrorSummary_label.setText("#.##")
                case "Solvent":
                    if not self.solventCA_autofit_done:
                        self.solventCA_deltaPhi0ErrorSummary_label.setText(
                            "#.##"
                        )
                        self.solventCA_n2ErrorSummary_label.setText("#.##")
                        self.solventOA_n2ErrorSummary_label.setText("#.##")
                        self.solventCA_beamwaistErrorSummary_label.setText(
                            "#.##"
                        )
                        self.solventCA_rayleighRangeErrorSummary_label.setText(
                            "#.##"
                        )
                    # if self.solventOA_autofit_done == False:
                    #     self.solventOA_TErrorSummary_label.setText("#.##")
                    #     self.solventOA_betaErrorSummary_label.setText('#.##')
                    #     self.solventOA_gammaErrorSummary_label.setText('#.##')
                case "Sample":
                    if not self.sampleCA_autofit_done:
                        self.sampleCA_deltaPhi0ErrorSummary_label.setText(
                            "#.##"
                        )
                        self.sampleCA_n2ErrorSummary_label.setText("#.##")
                        self.sampleOA_n2ErrorSummary_label.setText("#.##")
                        self.sampleCA_beamwaistErrorSummary_label.setText(
                            "#.##"
                        )
                        self.sampleCA_rayleighRangeErrorSummary_label.setText(
                            "#.##"
                        )
                    if not self.sampleOA_autofit_done:
                        self.sampleOA_TErrorSummary_label.setText("#.##")
                        self.sampleOA_betaErrorSummary_label.setText("#.##")
                        self.sampleOA_gammaErrorSummary_label.setText("#.##")

        # ROUND THE NUMBERS AND SET PRECISION OF THE DISPLAYED VALUES
        elif caller == "auto":
            try:
                match ftype:
                    case "Silica":
                        # --- decimals / precisions ---
                        dphi_prec = self.silicaCA_DPhi0Precision or 0
                        laser_prec = self.laserI0Precision or 0
                        bw_prec = self.silicaCA_beamwaistPrecision or 0
                        zr_prec = self.silica_rayleighLengthPrecision or 0
                        na_prec = self.numericalAperturePrecision or 0

                        self.silicaCA_deltaPhi0Summary_doubleSpinBox.setDecimals(
                            dphi_prec
                        )
                        self.silicaCA_laserIntensitySummary_doubleSpinBox.setDecimals(
                            laser_prec
                        )
                        self.silicaCA_beamwaistSummary_doubleSpinBox.setDecimals(
                            max(bw_prec - 6, 0)
                        )
                        self.silicaCA_rayleighRangeSummary_doubleSpinBox.setDecimals(
                            max(zr_prec - 3, 0)
                        )
                        self.numericalAperture_doubleSpinBox.setDecimals(
                            na_prec
                        )

                        # --- error labels ---
                        if (
                            self.silicaCA_DPhi0Error is not None
                            and dphi_prec > 0
                        ):
                            self.silicaCA_deltaPhi0ErrorSummary_label.setText(
                                f"{self.silicaCA_DPhi0Error:.{dphi_prec}f}"
                            )

                        if self.laserI0Error is not None and laser_prec > 0:
                            self.silicaCA_laserIntensityErrorSummary_label.setText(
                                f"{self.laserI0Error * LASER_INTENSITY_GW_SCALE:.{laser_prec}f}"
                            )

                        if (
                            self.silicaCA_beamwaistError is not None
                            and bw_prec > 0
                        ):
                            self.silicaCA_beamwaistErrorSummary_label.setText(
                                f"{self.silicaCA_beamwaistError * BEAMWAIST_UM_SCALE:.{max(bw_prec - 6, 0)}f}"
                            )

                        if (
                            self.silica_rayleighLengthError is not None
                            and zr_prec > 0
                        ):
                            self.silicaCA_rayleighRangeErrorSummary_label.setText(
                                f"{self.silica_rayleighLengthError * RAYLEIGH_LENGTH_MM_SCALE:.{max(zr_prec - 3, 0)}f}"
                            )

                        if (
                            self.numericalApertureError is not None
                            and na_prec > 0
                        ):
                            self.numericalApertureErrorSummary_label.setText(
                                f"{self.numericalApertureError:.{na_prec}f}"
                            )

                    case "Solvent":
                        (
                            self.solventCA_DPhi0,
                            self.solventCA_DPhi0Error,
                            self.solventCA_DPhi0Precision,
                        ) = error_rounding(
                            self.solventCA_minimizerResult.params[
                                "DPhi0"
                            ].value,
                            self.solventCA_minimizerResult.params[
                                "DPhi0"
                            ].stderr,
                        )
                        (
                            self.solventCA_beamwaist,
                            self.solventCA_beamwaistError,
                            self.solventCA_beamwaistPrecision,
                        ) = error_rounding(
                            self.solventCA_minimizerResult.params[
                                "Beamwaist"
                            ].value,
                            self.solventCA_minimizerResult.params[
                                "Beamwaist"
                            ].stderr,
                        )
                        (
                            self.solventCA_centerPoint,
                            self.solventCA_centerPointError,
                            self.solventCA_centerPointPrecision,
                        ) = error_rounding(
                            self.solventCA_minimizerResult.params[
                                "Center"
                            ].value,
                            self.solventCA_minimizerResult.params[
                                "Center"
                            ].stderr,
                        )
                        (
                            self.solventCA_zeroLevel,
                            self.solventCA_zeroLevelError,
                            self.solventCA_zeroLevelPrecision,
                        ) = error_rounding(
                            self.solventCA_minimizerResult.params["Zero"].value,
                            self.solventCA_minimizerResult.params[
                                "Zero"
                            ].stderr,
                        )

                        self.calculate_derived_parameters_errors(ftype)
                        (
                            self.solvent_n2,
                            self.solvent_n2Error,
                            self.solvent_n2Precision,
                        ) = error_rounding(
                            self.solvent_n2, self.solvent_n2_error
                        )
                        (
                            self.solvent_rayleighLength,
                            self.solvent_zR_error,
                            self.solvent_zR_precision,
                        ) = error_rounding(
                            self.solvent_n2, self.solvent_n2_error
                        )

                        # set display precision
                        if self.solventCA_DPhi0Precision:
                            self.solventCA_deltaPhi0Summary_doubleSpinBox.setDecimals(
                                self.solventCA_DPhi0Precision
                            )

                        if self.solvent_n2Precision:
                            self.solventCA_n2Summary_doubleSpinBox.setDecimals(
                                self.solvent_n2Precision
                            )
                            self.solventOA_n2Summary_doubleSpinBox.setDecimals(
                                self.solvent_n2Precision
                            )

                        if self.solventCA_beamwaistPrecision:
                            self.solventCA_beamwaistSummary_doubleSpinBox.setDecimals(
                                max(self.solventCA_beamwaistPrecision - 6, 0)
                            )

                        if self.solvent_zR_precision:
                            self.solventCA_rayleighRangeSummary_doubleSpinBox.setDecimals(
                                max(self.solvent_zR_precision - 3, 0)
                            )

                        # display error values
                        self.solventCA_deltaPhi0ErrorSummary_label.setText(
                            f"{self.solventCA_DPhi0Error:.{self.solventCA_DPhi0Precision}f}"
                        )
                        if (
                            self.solvent_n2_error is not None
                            and self.solvent_n2Precision is not None
                            and self.solvent_n2Precision > 0
                        ):
                            scaled_n2_err = self.solvent_n2_error * 1e13 * 1e12
                            self.solventCA_n2ErrorSummary_label.setText(
                                f"{scaled_n2_err:.{self.solvent_n2Precision}f}"
                            )
                            self.solventOA_n2ErrorSummary_label.setText(
                                f"{scaled_n2_err:.{self.solvent_n2Precision}f}"
                            )

                        if (
                            self.solventCA_beamwaistError is not None
                            and self.solventCA_beamwaistPrecision is not None
                            and self.solventCA_beamwaistPrecision > 0
                        ):
                            bw_prec = max(
                                self.solventCA_beamwaistPrecision - 6, 0
                            )
                            self.solventCA_beamwaistErrorSummary_label.setText(
                                f"{self.solventCA_beamwaistError * BEAMWAIST_UM_SCALE:.{bw_prec}f}"
                            )

                        if (
                            self.solvent_zR_error is not None
                            and self.solvent_zR_precision is not None
                            and self.solvent_zR_precision > 0
                        ):
                            zr_prec = max(self.solvent_zR_precision - 3, 0)
                            self.solventCA_rayleighRangeErrorSummary_label.setText(
                                f"{self.solvent_zR_error * RAYLEIGH_LENGTH_MM_SCALE:.{zr_prec}f}"
                            )

                    case "Sample":
                        pass

            except Exception as e:
                logging.error(traceback.format_exc())
                self.showdialog(
                    "Error",
                    "The fit converged, but another error occurred. Try using different initial parameters.",
                )
        else:
            pass

    def get_general_parameters(self) -> None:
        """Gets the values from 'General Parameters' GUI frame and assigns them to variables with basic SI units (mm -> m):\n
        `self.l_silica`, `self.lda`, `self.z_range`, `self.ra`, `self.d0`\n
        and then calculate parameters that depend only on them and not on fitting process:\n
        `self.silica_n2`
        """
        logger.debug("Window.get_general_parameters")
        from lib.validation import ParameterValidator

        l_silica = self.silicaThickness_dataFittingTab_doubleSpinBox.value()
        lda = self.wavelength_dataFittingTab_doubleSpinBox.value()
        z_range = self.zscanRange_doubleSpinBox.value()
        aperture = self.apertureDiameter_doubleSpinBox.value()

        # Validate
        valid, error = ParameterValidator.validate_general_parameters(
            lda, l_silica, aperture, z_range
        )
        if not valid:
            self.showdialog("Warning", f"Parameter validation: {error}")
            return

        # Convert UI values to basic SI units (meters)
        self.l_silica = l_silica * 1e-3  # [m] silica thickness
        self.lda = lda * 1e-9  # [m] wavelength
        self.z_range = z_range * 1e-3  # [m] z-scan range
        self.ra = aperture / 2 * 1e-3  # [m] CA aperture radius
        self.d0 = (
            self.apertureToFocusDistance_doubleSpinBox.value() * 1e-3
        )  # [m] distance from focal point to aperture plane
        self.silica_n2 = (
            2.8203e-20 - 3e-27 / (self.lda) + 2e-33 / (self.lda) ** 2
        )  # [m2/W] Bandar A. Babgi formula (based on David Milam tables for n2)
        logger.debug(f"""Window.get_general_parameters retrieved
        {self.l_silica=}\n{self.lda=}\n{self.z_range=}\n{self.ra=}\n{self.d0=}\n{self.silica_n2=}""")

    def get_curve_interpretation(
        self, ftype: str, stype: str, from_what: str, on_data_load: bool = False
    ) -> None:
        """Interprets the curve based on its geometry (manual fitting)
        or based on automatic fitting parameters (automatic fitting),
        and performs scientific error estimation.

        NOTE: Includes bounds checking to prevent unrealistic values from
        breaking the GUI display.
        """

        # Helper: safe rounding wrapper
        def safe_round(
            value: float, error: Optional[float]
        ) -> tuple[float, Optional[float], Optional[int]]:
            try:
                v, e, p = error_rounding(value, error)
                return float(v), e, p
            except Exception:
                return float(value), None, None

        def is_valid_fit(p):
            return all(
                p[name].value is not None and p[name].stderr is not None
                for name in ["Zero", "Center", "DPhi0", "Beamwaist"]
            )

        # Helper: clamp values to reasonable bounds to prevent display issues
        def clamp_to_bounds(
            value: float, name: str, ftype: str, stype: str
        ) -> float:
            """Clamp unrealistic values to keep GUI functional."""
            if name == "DPhi0":
                # DPhi0 should be reasonable - clamp to [-10, 10]
                clamped = max(min(value, 10.0), -10.0)
                if clamped != value:
                    logger.warning(
                        f"{ftype} {stype}: DPhi0={value} out of bounds, clamped to {clamped}"
                    )
                return clamped

            elif name == "beamwaist":
                # Beamwaist should be 1-200 µm, i.e., 1e-6 to 200e-6 m
                clamped = max(min(value, 200e-6), 1e-6)
                if clamped != value:
                    logger.warning(
                        f"{ftype} {stype}: beamwaist={value:.2e} out of bounds, clamped to {clamped:.2e}"
                    )
                return clamped

            elif name == "zero_level":
                # Zero level should be 0.5-2.0
                clamped = max(min(value, 2.0), 0.5)
                if clamped != value:
                    logger.warning(
                        f"{ftype} {stype}: zero_level={value} out of bounds, clamped to {clamped}"
                    )
                return clamped

            elif name == "centerpoint":
                # Center point should be within data range, roughly -50 to 50
                clamped = max(min(value, 50.0), -50.0)
                if clamped != value:
                    logger.warning(
                        f"{ftype} {stype}: centerpoint={value} out of bounds, clamped to {clamped}"
                    )
                return clamped

            return value

        logger.debug(f"""Window.get_curve_interpretation: running with params:\n" \
        {ftype=}, {stype=}, {from_what=}, {on_data_load=}""")
        match ftype:
            # ============================================================
            #                           SILICA
            # ============================================================
            case "Silica":
                match from_what:
                    # -------------------------
                    # Manual fitting
                    # -------------------------
                    case "from_geometry":
                        logger.debug(
                            "Window.get_curve_interpretation: reading from silicaCA_zeroLevel slider"
                        )
                        self.silicaCA_zeroLevel = (
                            self.silicaCA_zeroLevel_slider.value() / 100
                        )
                        logger.debug(
                            "Window.get_curve_interpretation: reading from silicaCA_centerPoint slider"
                        )
                        self.silicaCA_centerPoint = np.round(
                            self.silicaCA_centerPoint_slider.value()
                            - self.silica_nop / 2
                        )

                        if self.state._file_loading:
                            logger.debug(
                                f"Window.get_curve_interpretation: {self.state._file_loading=}, reading params_from_geometry"
                            )
                            (
                                self.silicaCA_DPhi0,
                                self.silicaCA_beamwaist,
                                self.silica_rayleighLength,
                            ) = self.params_from_geometry(ftype, "CA")
                        else:
                            logger.debug(
                                f"Window.get_curve_interpretation: {self.state._file_loading=}, reading from sliders"
                            )
                            self.silicaCA_DPhi0 = (
                                self.silicaCA_DPhi0_slider.value()
                                * MAX_DPHI0
                                / self.silicaCA_DPhi0_slider.maximum()
                            )
                            self.silica_rayleighLength = (
                                self.silicaCA_RayleighLength_slider.value()
                                * self.z_range
                                / 2
                                / self.silicaCA_RayleighLength_slider.maximum()
                            )
                            self.silicaCA_beamwaist = np.sqrt(
                                self.silica_rayleighLength * self.lda / np.pi
                            )
                        logger.debug(f"""Window.get_curve_interpretation: values read:\n
                        {self.silicaCA_DPhi0=}\n{self.silica_rayleighLength=}\n{self.silicaCA_beamwaist=}""")

                    # -------------------------
                    # Automatic fitting
                    # -------------------------
                    case "from_autofit":
                        p = self.silicaCA_minimizerResult.params

                        # Central values with bounds checking
                        self.silicaCA_zeroLevel = clamp_to_bounds(
                            p["Zero"].value, "zero_level", ftype, stype
                        )
                        self.silicaCA_centerPoint = clamp_to_bounds(
                            p["Center"].value, "centerpoint", ftype, stype
                        )
                        self.silicaCA_DPhi0 = clamp_to_bounds(
                            p["DPhi0"].value, "DPhi0", ftype, stype
                        )
                        self.silicaCA_beamwaist = clamp_to_bounds(
                            p["Beamwaist"].value, "beamwaist", ftype, stype
                        )
                        self.silica_rayleighLength = float(
                            np.pi * self.silicaCA_beamwaist**2 / self.lda
                        )

                # Common for both geometry/autofit
                self.laserI0 = (
                    self.silicaCA_DPhi0
                    * self.lda
                    / (2 * np.pi * self.l_silica * self.silica_n2)
                )
                logger.debug(f"{self.laserI0=}")

                try:
                    self.numericalAperture = (
                        self.silicaCA_beamwaist / self.silica_rayleighLength
                    )
                except ZeroDivisionError:
                    self.numericalAperture = 0
                logger.debug(f"{self.numericalAperture=}")

                # -------------------------
                # Error estimation
                # -------------------------
                match from_what:
                    case "from_geometry":
                        logger.debug(
                            f"Window.get_curve_interpretation: {from_what=} : omitting error estimation"
                        )
                        pass  # No errors for manual fitting

                    case "from_autofit":
                        p = self.silicaCA_minimizerResult.params

                        fit_valid = is_valid_fit(p)
                        if not fit_valid:
                            logger.warning(
                                f"{ftype} {stype}: Fit parameters missing standard errors"
                            )
                            return

                        # Zero level
                        (
                            self.silicaCA_zeroLevel,
                            self.silicaCA_zeroLevelError,
                            self.silicaCA_zeroLevelPrecision,
                        ) = safe_round(p["Zero"].value, p["Zero"].stderr)

                        # Center point
                        (
                            self.silicaCA_centerPoint,
                            self.silicaCA_centerPointError,
                            self.silicaCA_centerPointPrecision,
                        ) = safe_round(p["Center"].value, p["Center"].stderr)

                        # DPhi0
                        (
                            self.silicaCA_DPhi0,
                            self.silicaCA_DPhi0Error,
                            self.silicaCA_DPhi0Precision,
                        ) = safe_round(p["DPhi0"].value, p["DPhi0"].stderr)

                        # Beamwaist
                        (
                            self.silicaCA_beamwaist,
                            self.silicaCA_beamwaistError,
                            self.silicaCA_beamwaistPrecision,
                        ) = safe_round(
                            p["Beamwaist"].value, p["Beamwaist"].stderr
                        )

                        # Rayleigh length error
                        if self.silicaCA_beamwaistError is None:
                            self.silica_rayleighLengthError = None
                        else:
                            self.silica_rayleighLengthError = (
                                np.pi
                                / 2
                                / self.lda
                                * p["Beamwaist"].value
                                * self.silicaCA_beamwaistError
                            )

                        (
                            self.silica_rayleighLength,
                            self.silica_rayleighLengthError,
                            self.silica_rayleighLengthPrecision,
                        ) = safe_round(
                            self.silica_rayleighLength,
                            self.silica_rayleighLengthError,
                        )

                        # Laser intensity error
                        if self.silicaCA_DPhi0Error is None:
                            self.laserI0Error = None
                        else:
                            self.laserI0Error = (
                                self.silicaCA_DPhi0Error
                                * self.lda
                                / (2 * np.pi * self.l_silica * self.silica_n2)
                            )

                        # Convert to GW/cm² for rounding
                        I0_scaled = self.laserI0 * LASER_INTENSITY_GW_SCALE
                        I0err_scaled = (
                            None
                            if self.laserI0Error is None
                            else self.laserI0Error * LASER_INTENSITY_GW_SCALE
                        )

                        (
                            self.laserI0,
                            self.laserI0Error,
                            self.laserI0Precision,
                        ) = safe_round(I0_scaled, I0err_scaled)

                        # Convert back to W/m²
                        self.laserI0 *= 1e13
                        if self.laserI0Error is not None:
                            self.laserI0Error *= 1e13

                        # Numerical aperture error
                        if (
                            self.silicaCA_beamwaistError is None
                            or self.silica_rayleighLengthError is None
                        ):
                            self.numericalApertureError = None
                        else:
                            self.numericalApertureError = (
                                self.silicaCA_beamwaistError
                                / self.silica_rayleighLength
                                + self.silicaCA_beamwaist
                                * self.silica_rayleighLengthError
                                / self.silica_rayleighLength**2
                            )

                        (
                            self.numericalAperture,
                            self.numericalApertureError,
                            self.numericalAperturePrecision,
                        ) = safe_round(
                            self.numericalAperture, self.numericalApertureError
                        )

            # ============================================================
            #                           SOLVENT
            # ============================================================
            case "Solvent":
                match from_what:
                    case "from_geometry":
                        self.solventCA_zeroLevel = (
                            self.solventCA_zeroLevel_slider.value() / 100
                        )
                        self.solventCA_centerPoint = np.round(
                            self.solventCA_centerPoint_slider.value()
                            - self.solvent_nop / 2
                        )

                        if on_data_load:
                            (
                                self.solventCA_DPhi0,
                                self.solventCA_beamwaist,
                                self.solvent_rayleighLength,
                            ) = self.params_from_geometry(ftype, "CA")
                        else:
                            self.solventCA_DPhi0 = (
                                self.solventCA_DPhi0_slider.value()
                                * MAX_DPHI0
                                / self.solventCA_DPhi0_slider.maximum()
                            )
                            self.solvent_rayleighLength = (
                                self.solventCA_RayleighLength_slider.value()
                                * self.z_range
                                / 2
                                / self.solventCA_RayleighLength_slider.maximum()
                            )
                            self.solventCA_beamwaist = np.sqrt(
                                self.solvent_rayleighLength * self.lda / np.pi
                            )

                        # If OA sliders are present and user activated custom center, read it
                        try:
                            if (
                                hasattr(self, "solventOA_centerPoint_slider")
                                and self.solventOA_customCenterPoint_checkBox.isChecked()
                            ):
                                self.solventOA_centerPoint = (
                                    self.solventOA_centerPoint_slider.value()
                                    - SOLVENT_OA_SLIDER_CENTER_OFFSET
                                )
                        except Exception:
                            pass

                        if not self.solventCA_customBeamwaist_checkBox.isChecked():
                            self.solventCA_beamwaist = self.silicaCA_beamwaist
                            self.solvent_rayleighLength = (
                                self.silica_rayleighLength
                            )

                    case "from_autofit":
                        p = self.solventCA_minimizerResult.params

                        fit_valid = is_valid_fit(p)
                        if not fit_valid:
                            logger.warning(
                                f"{ftype} {stype}: Fit parameters missing standard errors"
                            )
                            return

                        # Central values with bounds checking
                        self.solventCA_zeroLevel = clamp_to_bounds(
                            p["Zero"].value, "zero_level", ftype, stype
                        )
                        self.solventCA_centerPoint = clamp_to_bounds(
                            p["Center"].value, "centerpoint", ftype, stype
                        )
                        self.solventCA_DPhi0 = clamp_to_bounds(
                            p["DPhi0"].value, "DPhi0", ftype, stype
                        )
                        self.solventCA_beamwaist = clamp_to_bounds(
                            p["Beamwaist"].value, "beamwaist", ftype, stype
                        )
                        self.solvent_rayleighLength = float(
                            np.pi * self.solventCA_beamwaist**2 / self.lda
                        )

                        # For OA: sync centerpoint if custom is disabled
                        try:
                            if (
                                hasattr(
                                    self, "solventOA_customCenterPoint_checkBox"
                                )
                                and not self.solventOA_customCenterPoint_checkBox.isChecked()
                            ):
                                self.solventOA_centerPoint = (
                                    self.solventCA_centerPoint
                                )
                        except Exception:
                            pass

                # Solvent nonlinear index
                self.solvent_n2 = (
                    self.solventCA_DPhi0
                    / self.silicaCA_DPhi0
                    * self.silica_n2
                    * self.l_silica
                    / 0.001
                )

                # Error estimation
                match from_what:
                    case "from_geometry":
                        pass

                    case "from_autofit":
                        p = self.solventCA_minimizerResult.params

                        (
                            self.solventCA_zeroLevel,
                            self.solventCA_zeroLevelError,
                            self.solventCA_zeroLevelPrecision,
                        ) = safe_round(p["Zero"].value, p["Zero"].stderr)

                        (
                            self.solventCA_centerPoint,
                            self.solventCA_centerPointError,
                            self.solventCA_centerPointPrecision,
                        ) = safe_round(p["Center"].value, p["Center"].stderr)

                        (
                            self.solventCA_DPhi0,
                            self.solventCA_DPhi0Error,
                            self.solventCA_DPhi0Precision,
                        ) = safe_round(p["DPhi0"].value, p["DPhi0"].stderr)

                        (
                            self.solventCA_beamwaist,
                            self.solventCA_beamwaistError,
                            self.solventCA_beamwaistPrecision,
                        ) = safe_round(
                            p["Beamwaist"].value, p["Beamwaist"].stderr
                        )

                        # Rayleigh length error
                        if self.solventCA_beamwaistError is None:
                            self.solvent_rayleighLengthError = None
                        else:
                            self.solvent_rayleighLengthError = (
                                np.pi
                                / 2
                                / self.lda
                                * p["Beamwaist"].value
                                * self.solventCA_beamwaistError
                            )

                        (
                            self.solvent_rayleighLength,
                            self.solvent_rayleighLengthError,
                            self.solvent_rayleighLengthPrecision,
                        ) = safe_round(
                            self.solvent_rayleighLength,
                            self.solvent_rayleighLengthError,
                        )

            # ============================================================
            #                           SAMPLE
            # ============================================================
            case "Sample":
                if self.sampleCA_fittingLine_drawn:
                    self.sampleCA_zeroLevel = (
                        self.sampleCA_zeroLevel_slider.value() / 100
                    )
                    self.sampleCA_centerPoint = (
                        self.sampleCA_centerPoint_slider.value()
                        - SOLVENT_CA_SLIDER_CENTER_OFFSET
                    )

                if self.sampleOA_fittingLine_drawn:
                    self.sampleOA_zeroLevel = (
                        self.sampleOA_zeroLevel_slider.value() / 100
                    )
                    self.sampleOA_centerPoint = (
                        self.sampleOA_centerPoint_slider.value()
                        - SOLVENT_OA_SLIDER_CENTER_OFFSET
                    )

    def params_from_geometry(self, ftype, stype):
        """Using equations given by van Stryland and Sheik-Bahae in
        http://www.phys.unm.edu/msbahae/publications/z-scan.pdf
        return the fitting curve physical interpretation."""
        match ftype:
            case "Silica":
                if stype == "CA":
                    dataset = self.silica_data_set[0], self.silica_data_set[1]
                else:
                    dataset = self.silica_data_set[0], self.silica_data_set[3]
            case "Solvent":
                if stype == "CA":
                    dataset = self.solvent_data_set[0], self.solvent_data_set[1]
                else:
                    dataset = self.solvent_data_set[0], self.solvent_data_set[3]
            case "Sample":
                if stype == "CA":
                    dataset = self.sample_data_set[0], self.sample_data_set[1]
                else:
                    dataset = self.sample_data_set[0], self.sample_data_set[3]

        fit_x, fit_y = dataset[0], dataset[1]

        if stype == "CA":
            try:
                ymax_pos = np.where(fit_y == np.max(fit_y))[0][0]
                ymin_pos = np.where(fit_y == np.min(fit_y))[0][0]

                if ymax_pos > ymin_pos:
                    deltaTpv_sign = 1
                elif ymax_pos < ymin_pos:
                    deltaTpv_sign = -1
                else:
                    deltaTpv_sign = 0

                deltaTpv = deltaTpv_sign * abs(np.max(fit_y) - np.min(fit_y))
                deltaPhi0 = deltaTpv / 0.406
                deltaZpv = abs(fit_x[ymax_pos] - fit_x[ymin_pos])  # [m]
                rayleighLength = deltaZpv / 1.7  # [m] Rayleigh length
                beamwaist = np.sqrt(
                    rayleighLength * self.lda / np.pi
                )  # [m] beam radius in focal point

                return deltaPhi0, beamwaist, rayleighLength

            except TypeError:
                print(
                    "Something is wrong while interpreting the closed aperture curve."
                )

        elif stype == "OA":
            pass

    def set_sliders_positions(self, ftype, stype):
        """Called by `self.fit_automatically`.\n
        This method recalculates value to reflect given slider position.
        """
        match ftype:
            case "Silica":
                # Silica CA sliders
                # Because change in slider value triggers the 'fit_manually' method, it has to be disabled.
                self.silicaCA_RayleighLength_slider.valueChanged.disconnect()
                self.silicaCA_centerPoint_slider.valueChanged.disconnect()
                self.silicaCA_zeroLevel_slider.valueChanged.disconnect()
                self.silicaCA_DPhi0_slider.valueChanged.disconnect()

                self.silicaCA_RayleighLength_slider.setValue(
                    int(
                        round(
                            self.silica_rayleighLength
                            * self.silicaCA_RayleighLength_slider.maximum()
                            / (self.z_range / 2)
                        )
                    )
                )
                self.silicaCA_centerPoint_slider.setValue(
                    int(
                        round(
                            self.silicaCA_centerPoint
                            + self.silicaCA_centerPoint_slider.maximum() / 2
                        )
                    )
                )
                self.silicaCA_zeroLevel_slider.setValue(
                    int(round(self.silicaCA_zeroLevel * 100))
                )
                self.silicaCA_DPhi0_slider.setValue(
                    int(
                        round(
                            self.silicaCA_DPhi0
                            / np.pi
                            * self.silicaCA_DPhi0_slider.maximum()
                        )
                    )
                )

                # And now when all is updated by the 'fit_automatically', reconnect the sliders to their slots
                self.slider_triggers()

            case "Solvent":
                if stype == "CA":
                    # Solvent CA sliders
                    # Because change in slider value triggers the 'fit_manually' method, it has to be disabled.
                    self.solventCA_RayleighLength_slider.valueChanged.disconnect()
                    self.solventCA_centerPoint_slider.valueChanged.disconnect()
                    self.solventCA_zeroLevel_slider.valueChanged.disconnect()
                    self.solventCA_DPhi0_slider.valueChanged.disconnect()

                    if not self.solventCA_customBeamwaist_checkBox.isChecked():
                        self.solventCA_RayleighLength_slider.setValue(
                            self.silicaCA_RayleighLength_slider.value()
                        )
                    else:
                        self.solventCA_RayleighLength_slider.setValue(
                            int(
                                round(
                                    self.solventCA_beamwaist
                                    * BEAMWAIST_UM_SCALE
                                )
                            )
                        )
                    self.solventCA_centerPoint_slider.setValue(
                        int(
                            round(
                                self.solventCA_centerPoint
                                + SOLVENT_CA_SLIDER_CENTER_OFFSET
                            )
                        )
                    )
                    self.solventCA_zeroLevel_slider.setValue(
                        int(round(self.solventCA_zeroLevel * 100))
                    )
                    self.solventCA_DPhi0_slider.setValue(
                        int(
                            round(
                                self.solventCA_DPhi0
                                * 0.406
                                * self.solventCA_DPhi0_slider.maximum()
                                / 5
                            )
                        )
                    )

                    # And now when all is updated by the 'fit_automatically', reconnect the sliders to their slots
                    self.solventCA_RayleighLength_slider.valueChanged.connect(
                        lambda: self.fit_manually(ftype="Solvent", stype="CA")
                    )
                    self.solventCA_centerPoint_slider.valueChanged.connect(
                        lambda: self.fit_manually(ftype="Solvent", stype="CA")
                    )
                    self.solventCA_zeroLevel_slider.valueChanged.connect(
                        lambda: self.fit_manually(ftype="Solvent", stype="CA")
                    )
                    self.solventCA_DPhi0_slider.valueChanged.connect(
                        lambda: self.fit_manually(ftype="Solvent", stype="CA")
                    )

                elif stype == "OA":
                    # Solvent OA sliders
                    # Because change in slider value triggers the 'fit_manually' method, it has to be disabled.
                    self.solventOA_centerPoint_slider.valueChanged.disconnect()
                    self.solventOA_zeroLevel_slider.valueChanged.disconnect()
                    self.solventOA_T_slider.valueChanged.disconnect()

                    self.solventOA_centerPoint_slider.setValue(
                        int(
                            round(
                                self.solventOA_centerPoint
                                + SOLVENT_OA_SLIDER_CENTER_OFFSET
                            )
                        )
                    )
                    self.solventOA_zeroLevel_slider.setValue(
                        int(round(self.solventOA_zeroLevel * 100))
                    )
                    # Map physical T (which may have arbitrary min/max) into slider range
                    try:
                        t_min = OA_FITTING_PARAMS["T"]["min"]
                        t_max = OA_FITTING_PARAMS["T"]["max"]
                        if t_max == t_min:
                            mapped = 0
                        else:
                            mapped = (
                                (self.solventOA_T - t_min)
                                / (t_max - t_min)
                                * self.solventOA_T_slider.maximum()
                            )
                        self.solventOA_T_slider.setValue(int(round(mapped)))
                    except Exception:
                        # Fallback to legacy behaviour if config missing
                        self.solventOA_T_slider.setValue(
                            int(
                                round(
                                    self.solventOA_T
                                    * self.solventOA_T_slider.maximum()
                                    / SOLVENT_T_SLIDER_MAX
                                )
                            )
                        )

                    # And now when all is updated by the 'fit_automatically', reconnect the sliders to their slots
                    self.solventOA_centerPoint_slider.valueChanged.connect(
                        self.solventOA_centerPoint_slider_moved
                    )
                    self.solventOA_zeroLevel_slider.valueChanged.connect(
                        lambda: self.fit_manually(ftype="Solvent", stype="OA")
                    )
                    # connect T slider to dedicated handler that updates spinbox then triggers manual fit
                    self.solventOA_T_slider.valueChanged.connect(
                        self.solventOA_T_slider_moved
                    )

            case "Sample":
                pass

    def calculate_derived_parameters(self, ftype):
        '''Calculates `laserI0` for `ftype`="Silica", `n2` and `rayleigh length` for `ftype`="Solvent/Sample"'''

        match ftype:
            case "Silica":
                pass  # posprzątane
            case "Solvent":
                self.solvent_rayleighLength = (
                    np.pi / self.lda * self.solventCA_beamwaist**2
                )  # [m] Rayleigh length

    def calculate_derived_parameters_errors(self, ftype):
        match ftype:
            case "Silica":
                pass  # posprzątane
            case "Solvent":
                if hasattr(self, "solventCA_DPhi0Error"):
                    self.solvent_zR_error = (
                        np.pi
                        / 2
                        / self.lda
                        * self.solventCA_minimizerResult.params[
                            "Beamwaist"
                        ].value
                        * self.solventCA_beamwaistError
                    )  # [m] Rayleigh length
                    self.solvent_n2_error = (
                        self.solventCA_DPhi0Error
                        / self.silicaCA_minimizerResult.params["DPhi0"].value
                        * self.silica_n2
                        * self.l_silica
                        + self.solventCA_minimizerResult.params["DPhi0"].value
                        * self.silicaCA_DPhi0Error
                        / self.silicaCA_minimizerResult.params["DPhi0"].value
                        ** 2
                        * self.silica_n2
                        * self.l_silica
                    )  # [m2/W]

                else:
                    self.solvent_zR_error = 0
                    self.solvent_n2_error = 0
                    self.showdialog(
                        "Error",
                        "Errors were not estimated. Try different region of interest.",
                    )

                self.solvent_beta_error = None

            case "Sample":
                pass

    def data_display(self, data_set, ftype: str):
        """Called by data_loader(). Displays:\n
        1) CA column divided by OA column to remove influence of OA on CA as 'CA'
        2) OA column as 'OA'\n
        to separate data processing for n2 and beta nonlinear coefficients."""
        self.get_general_parameters()

        def basic_data_manipulation(
            data: tuple[list[float], list[float], list[float], list[float]],
        ) -> Tuple[
            NDArray[np.float64],
            NDArray[np.float64],
            NDArray[np.float64],
            NDArray[np.float64],
        ]:
            """Reference division, separation of OA signal from CA fluctuation signal, normalization

            Args:
                data (any): four-column data containing datapoints, CA data, Reference data and OA data

            Returns:
                Tuple: Data ready for separated determination of n2 and beta nonlinear coefficients
            """
            positions0, ca_data0, ref_data0, oa_data0 = data

            # Divide data by reference
            ca_list = [ca / ref for ca, ref in zip(ca_data0, ref_data0)]
            oa_list = [oa / ref for oa, ref in zip(oa_data0, ref_data0)]

            # centralize positions
            nop = len(positions0)
            positions: NDArray[np.float64] = np.array(
                [
                    (self.z_range * zz / nop - self.z_range / 2)
                    for zz in range(nop)
                ],
                dtype=float,
            )

            def normalize(arr: list[float]) -> NDArray[np.float64]:
                zero_lvl = (
                    np.mean(arr[0:10]) + np.mean(arr[len(arr) - 10 :])
                ) / 2
                out: NDArray[np.float64] = (
                    np.asarray(arr, dtype=float) / zero_lvl
                )
                return out

            ca_data: NDArray[np.float64] = normalize(ca_list)
            oa_data: NDArray[np.float64] = normalize(oa_list)
            ref_data: NDArray[np.float64] = np.asarray(ref_data0, dtype=float)

            return positions, ca_data, ref_data, oa_data

        self.positions, self.ca_data, self.ref_data, self.oa_data = (
            basic_data_manipulation(data_set)
        )

        # Create reference-corrected 'data_set' variable for each ftype for further reference
        match ftype:
            case "Silica":
                self.silica_data_set = [
                    self.positions,
                    self.ca_data,
                    self.ref_data,
                    self.oa_data,
                ]
                self.silica_nop = len(self.positions)
            case "Solvent":
                self.solvent_data_set = [
                    self.positions,
                    self.ca_data,
                    self.ref_data,
                    self.oa_data,
                ]
                self.solvent_nop = len(self.positions)
            case "Sample":
                self.sample_data_set = [
                    self.positions,
                    self.ca_data,
                    self.ref_data,
                    self.oa_data,
                ]
                self.sample_nop = len(self.positions)

        # Display data on proper figures
        # Use reference-corrected data for fitting
        match ftype:
            case "Silica":
                self.update_datafitting_plotlimits(
                    self.silica_data_set, ftype=ftype
                )

            case "Solvent":
                self.update_datafitting_plotlimits(
                    self.solvent_data_set, ftype=ftype
                )

            case "Sample":
                self.update_datafitting_plotlimits(
                    self.sample_data_set, ftype=ftype
                )

    def set_new_positions(self):
        """Takes current values in 'General Parameters' from GUI and updates data_set[0] (positions) to meet new Z-scan range.\n
        Updates fitting plots limits to visualize the change.
        """
        self.get_general_parameters()

        if hasattr(self, "silica_data_set"):
            self.silica_data_set[0] = np.array(
                [
                    (self.z_range * zz / self.silica_nop - self.z_range / 2)
                    * 1000
                    for zz in range(self.silica_nop)
                ]
            )  # [mm] update positions with newly-read Z-scan range value
            self.update_datafitting_plotlimits(
                self.silica_data_set, ftype="Silica"
            )

        if hasattr(self, "solvent_data_set"):
            self.solvent_data_set[0] = np.array(
                [
                    (self.z_range * zz / self.solvent_nop - self.z_range / 2)
                    * 1000
                    for zz in range(self.solvent_nop)
                ]
            )  # [mm] update positions with newly-read Z-scan range value
            self.update_datafitting_plotlimits(
                self.solvent_data_set, ftype="Solvent"
            )

        if hasattr(self, "sample_data_set"):
            self.sample_data_set[0] = np.array(
                [
                    (self.z_range * zz / self.sample_nop - self.z_range / 2)
                    * 1000
                    for zz in range(self.sample_nop)
                ]
            )  # [mm] update positions with newly-read Z-scan range value
            self.update_datafitting_plotlimits(
                self.sample_data_set, ftype="Sample"
            )

    def update_datafitting_plotlimits(self, data_set: list, ftype: str) -> None:
        """Updates limits of the data fitting plots for given `ftype` (Silica, Solvent, Sample)

        Args:
            data_set (list): four-column reference-corrected data
            ftype (str): parameter holding information of sample type (Silica, Solvent, Sample)
        """
        self.get_general_parameters()  # Ensures working with currently typed-in General Parameters values from GUI
        padding_vertical = 0.01

        def set_limits(axes, line, direction, padding):
            if direction == "vertical":
                axes.set_ylim(
                    top=np.max(line.get_ydata()) * (1 + padding),
                    bottom=np.min(line.get_ydata()) * (1 - padding),
                )

        match ftype:
            case "Silica":
                # Closed aperture
                line = self.silicaCA_figure.axes.get_lines()[0]
                line.set_xdata(data_set[0] * POSITIONS_MM_SCALE)
                # line.set_ydata(ca_data)
                line.set_ydata(
                    [ca / oa for ca, oa in zip(data_set[1], data_set[3])]
                )

                self.silicaCA_figure.axes.set_xlim(
                    left=-self.z_range / 2 * POSITIONS_MM_SCALE,
                    right=self.z_range / 2 * POSITIONS_MM_SCALE,
                )  # displayed in mm
                set_limits(
                    self.silicaCA_figure.axes,
                    line,
                    "vertical",
                    padding_vertical,
                )

                # Open aperture
                line = self.silicaOA_figure.axes.get_lines()[0]
                line.set_xdata(data_set[0] * POSITIONS_MM_SCALE)
                line.set_ydata(data_set[3])

                self.silicaOA_figure.axes.set_xlim(
                    left=-self.z_range / 2 * POSITIONS_MM_SCALE,
                    right=self.z_range / 2 * POSITIONS_MM_SCALE,
                )  # displayed in mm
                set_limits(
                    self.silicaOA_figure.axes,
                    line,
                    "vertical",
                    padding_vertical,
                )
                # self.silicaOA_figure.axes.set_ylim(top=np.max(line.get_ydata())*(1+margin_vertical),bottom=np.min(line.get_ydata())*(1-margin_vertical))

                # Update
                self.silicaCA_figure.axes.relim()
                self.silicaOA_figure.axes.relim()

                self.silicaCA_figure.axes.autoscale_view()
                self.silicaOA_figure.axes.autoscale_view()

                self.silicaCA_figure.draw_idle()
                self.silicaOA_figure.draw_idle()

            case "Solvent":
                # Closed aperture
                line = self.solventCA_figure.axes.get_lines()[0]
                line.set_xdata(data_set[0] * POSITIONS_MM_SCALE)
                # line.set_ydata(ca_data)
                line.set_ydata(
                    [ca / oa for ca, oa in zip(data_set[1], data_set[3])]
                )

                self.solventCA_figure.axes.set_xlim(
                    left=-self.z_range / 2 * POSITIONS_MM_SCALE,
                    right=self.z_range / 2 * POSITIONS_MM_SCALE,
                )  # displayed in mm
                set_limits(
                    self.solventCA_figure.axes,
                    line,
                    "vertical",
                    padding_vertical,
                )
                # self.solventCA_figure.axes.set_ylim(top=np.max(line.get_ydata())*(1+margin_vertical),bottom=np.min(line.get_ydata())*(1-margin_vertical))

                # Open aperture
                line = self.solventOA_figure.axes.get_lines()[0]
                line.set_xdata(data_set[0] * POSITIONS_MM_SCALE)
                line.set_ydata(data_set[3])

                self.solventOA_figure.axes.set_xlim(
                    left=-self.z_range / 2 * POSITIONS_MM_SCALE,
                    right=self.z_range / 2 * POSITIONS_MM_SCALE,
                )  # displayed in mm
                set_limits(
                    self.solventOA_figure.axes,
                    line,
                    "vertical",
                    padding_vertical,
                )
                # self.solventOA_figure.axes.set_ylim(top=np.max(line.get_ydata())*(1+margin_vertical),bottom=np.min(line.get_ydata())*(1-margin_vertical))

                # Update
                self.solventCA_figure.axes.relim()
                self.solventOA_figure.axes.relim()

                self.solventCA_figure.axes.autoscale_view()
                self.solventOA_figure.axes.autoscale_view()

                self.solventCA_figure.draw_idle()
                self.solventOA_figure.draw_idle()

            case "Sample":
                # Closed aperture
                line = self.sampleCA_figure.axes.get_lines()[0]
                line.set_xdata(data_set[0] * POSITIONS_MM_SCALE)
                # line.set_ydata(ca_data)
                line.set_ydata(
                    [ca / oa for ca, oa in zip(data_set[1], data_set[3])]
                )

                self.sampleCA_figure.axes.set_xlim(
                    left=-self.z_range / 2 * POSITIONS_MM_SCALE,
                    right=self.z_range / 2 * POSITIONS_MM_SCALE,
                )  # displayed in mm
                set_limits(
                    self.sampleCA_figure.axes,
                    line,
                    "vertical",
                    padding_vertical,
                )
                # self.sampleCA_figure.axes.set_ylim(top=np.max(line.get_ydata())*1.04,bottom=np.min(line.get_ydata())*0.96)

                # Open aperture
                line = self.sampleOA_figure.axes.get_lines()[0]
                line.set_xdata(data_set[0] * POSITIONS_MM_SCALE)
                line.set_ydata(data_set[3])

                self.sampleOA_figure.axes.set_xlim(
                    left=-self.z_range / 2 * POSITIONS_MM_SCALE,
                    right=self.z_range / 2 * POSITIONS_MM_SCALE,
                )  # displayed in mm
                set_limits(
                    self.sampleOA_figure.axes,
                    line,
                    "vertical",
                    padding_vertical,
                )
                # self.sampleOA_figure.axes.set_ylim(top=np.max(line.get_ydata())*1.04,bottom=np.min(line.get_ydata())*0.96)

                # Update
                self.sampleCA_figure.axes.relim()
                self.sampleOA_figure.axes.relim()

                self.sampleCA_figure.axes.autoscale_view()
                self.sampleOA_figure.axes.autoscale_view()

                self.sampleCA_figure.draw_idle()
                self.sampleOA_figure.draw_idle()

    def reduce_noise_in_data(self, data_set, ftype, stype) -> None:
        match ftype:
            case "Silica":
                if stype == "CA":
                    filter_size = self.silicaCA_filterSize_slider.value()
                elif stype == "OA":
                    filter_size = self.silicaOA_filterSize_slider.value()
            case "Solvent":
                if stype == "CA":
                    filter_size = self.solventCA_filterSize_slider.value()
                elif stype == "OA":
                    filter_size = self.solventOA_filterSize_slider.value()
            case "Sample":
                if stype == "CA":
                    filter_size = self.sampleCA_filterSize_slider.value()
                elif stype == "OA":
                    filter_size = self.sampleOA_filterSize_slider.value()

        if (filter_size % 2 == 1 or filter_size == 0) and data_set is not None:
            try:
                match stype:
                    case "CA":
                        if filter_size > 0:
                            filtered_y = medfilt(data_set[1], filter_size)
                        else:
                            filtered_y = data_set[1]
                    case "OA":
                        if filter_size > 0:
                            filtered_y = medfilt(data_set[3], filter_size)
                        else:
                            filtered_y = data_set[3]

                match ftype:
                    case "Silica":
                        if stype == "CA":
                            line = self.silicaCA_figure.axes.get_lines()[0]
                            line.set_ydata(filtered_y)

                            self.silicaCA_figure.axes.relim()
                            self.silicaCA_figure.axes.autoscale_view()
                            self.silicaCA_figure.draw_idle()
                        elif stype == "OA":
                            pass

                    case "Solvent":
                        if stype == "CA":
                            line = self.solventCA_figure.axes.get_lines()[0]
                            line.set_ydata(filtered_y)

                            self.solventCA_figure.axes.relim()
                            self.solventCA_figure.axes.autoscale_view()
                            self.solventCA_figure.draw_idle()
                        elif stype == "OA":
                            line = self.solventOA_figure.axes.get_lines()[0]
                            line.set_ydata(filtered_y)

                            self.solventOA_figure.axes.relim()
                            self.solventOA_figure.axes.autoscale_view()
                            self.solventOA_figure.draw_idle()

                    case "Sample":
                        pass

            except ValueError:
                self.showdialog(
                    "Error",
                    "Possibly too few data points!\nMinimum required is 16 datapoints.",
                )

    def enable_cursors(self, ftype: str, stype: str) -> None:
        match ftype:
            case "Silica":
                if stype == "CA":
                    if self.silicaCA_fixROI_checkBox.isChecked():
                        self.silicaCA_cursor_positioner = BlittedCursor(
                            self.silicaCA_figure.axes,
                            color="magenta",
                            linewidth=2,
                        )
                        self.on_mouse_move = self.silicaCA_figure.mpl_connect(
                            "motion_notify_event",
                            self.silicaCA_cursor_positioner.on_mouse_move,
                        )
                        self.on_mouse_click = self.silicaCA_figure.mpl_connect(
                            "button_press_event",
                            lambda event: self.collect_cursor_clicks(
                                event, ftype
                            ),
                        )
                    else:
                        try:
                            self.silicaCA_cursor_positioner.vertical_line.remove()
                            self.silicaCA_cursor_positioner.horizontal_line.remove()
                        except (AttributeError, ValueError):
                            print("Specified cross-hair doesn't exist.")
                        finally:
                            self.silicaCA_figure.mpl_disconnect(
                                self.on_mouse_move
                            )
                            self.silicaCA_figure.mpl_disconnect(
                                self.on_mouse_click
                            )
                            self.silicaCA_cursorPositions.clear()
                        try:
                            self.silicaCA_verline1.remove()
                            self.silicaCA_verline2.remove()
                        except (AttributeError, ValueError):
                            print("Specified cross-hair doesn't exist.")
                        finally:
                            self.silicaCA_figure.draw_idle()

                elif stype == "OA":
                    pass
            case "Solvent":
                if stype == "CA":
                    if self.solventCA_fixROI_checkBox.isChecked():
                        self.solventCA_cursor_positioner = BlittedCursor(
                            self.solventCA_figure.axes,
                            color="magenta",
                            linewidth=2,
                        )
                        self.on_mouse_move = self.solventCA_figure.mpl_connect(
                            "motion_notify_event",
                            self.solventCA_cursor_positioner.on_mouse_move,
                        )
                        self.on_mouse_click = self.solventCA_figure.mpl_connect(
                            "button_press_event",
                            lambda event: self.collect_cursor_clicks(
                                event, ftype, stype
                            ),
                        )
                    else:
                        try:
                            self.solventCA_cursor_positioner.vertical_line.remove()
                            self.solventCA_cursor_positioner.horizontal_line.remove()
                        except (AttributeError, ValueError):
                            print("Specified cross-hair doesn't exist.")
                        finally:
                            self.solventCA_figure.mpl_disconnect(
                                self.on_mouse_move
                            )
                            self.solventCA_figure.mpl_disconnect(
                                self.on_mouse_click
                            )
                            self.solventCA_cursorPositions.clear()
                        try:
                            self.solventCA_verline1.remove()
                            self.solventCA_verline2.remove()
                        except (AttributeError, ValueError):
                            print("Specified cross-hair doesn't exist.")
                        finally:
                            self.solventCA_figure.draw_idle()

                elif stype == "OA":
                    if self.solventOA_fixROI_checkBox.isChecked():
                        self.solventOA_cursor_positioner = BlittedCursor(
                            self.solventOA_figure.axes,
                            color="magenta",
                            linewidth=2,
                        )
                        self.on_mouse_move = self.solventOA_figure.mpl_connect(
                            "motion_notify_event",
                            self.solventOA_cursor_positioner.on_mouse_move,
                        )
                        self.on_mouse_click = self.solventOA_figure.mpl_connect(
                            "button_press_event",
                            lambda event: self.collect_cursor_clicks(
                                event, ftype, stype
                            ),
                        )
                    else:
                        try:
                            self.solventOA_cursor_positioner.vertical_line.remove()
                            self.solventOA_cursor_positioner.horizontal_line.remove()
                        except (AttributeError, ValueError):
                            print("Specified cross-hair doesn't exist.")
                        finally:
                            self.solventOA_figure.mpl_disconnect(
                                self.on_mouse_move
                            )
                            self.solventOA_figure.mpl_disconnect(
                                self.on_mouse_click
                            )
                            self.solventOA_cursorPositions.clear()
                        try:
                            self.solventOA_verline1.remove()
                            self.solventOA_verline2.remove()
                        except (AttributeError, ValueError):
                            print("Specified cross-hair doesn't exist.")
                        finally:
                            self.solventOA_figure.draw_idle()

            case "Sample":
                if stype == "CA":
                    self.sampleCA_cursor_positioner = BlittedCursor(
                        self.sampleCA_figure.axes, color="magenta", linewidth=2
                    )
                elif stype == "OA":
                    self.sampleOA_cursor_positioner = BlittedCursor(
                        self.sampleOA_figure.axes, color="magenta", linewidth=2
                    )

    def collect_cursor_clicks(self, event, ftype: str, stype="") -> None:
        x, y = event.xdata, event.ydata
        match ftype:
            case "Silica":
                if len(self.silicaCA_cursorPositions) < 2:
                    self.silicaCA_cursorPositions.append((x, y))
                    if len(self.silicaCA_cursorPositions) == 1:
                        self.silicaCA_verline1 = (
                            self.silicaCA_figure.axes.axvline(
                                x, color="orange", linewidth=2
                            )
                        )
                    else:
                        self.silicaCA_verline2 = (
                            self.silicaCA_figure.axes.axvline(
                                x, color="orange", linewidth=2
                            )
                        )
                else:
                    self.silicaCA_cursorPositions.clear()
                    self.silicaCA_cursorPositions.append((x, y))
                    self.silicaCA_verline1.set_xdata(x)
                    self.silicaCA_verline2.set_xdata(None)

                self.silicaCA_figure.draw_idle()

            case "Solvent":
                if stype == "CA":
                    if len(self.solventCA_cursorPositions) < 2:
                        self.solventCA_cursorPositions.append((x, y))
                        if len(self.solventCA_cursorPositions) == 1:
                            self.solventCA_verline1 = (
                                self.solventCA_figure.axes.axvline(
                                    x, color="orange", linewidth=2
                                )
                            )
                        else:
                            self.solventCA_verline2 = (
                                self.solventCA_figure.axes.axvline(
                                    x, color="orange", linewidth=2
                                )
                            )
                    else:
                        self.solventCA_cursorPositions.clear()
                        self.solventCA_cursorPositions.append((x, y))
                        self.solventCA_verline1.set_xdata(x)
                        self.solventCA_verline2.set_xdata(None)

                    self.solventCA_figure.draw_idle()

                elif stype == "OA":
                    if len(self.solventOA_cursorPositions) < 2:
                        self.solventOA_cursorPositions.append((x, y))
                        if len(self.solventOA_cursorPositions) == 1:
                            self.solventOA_verline1 = (
                                self.solventOA_figure.axes.axvline(
                                    x, color="orange", linewidth=2
                                )
                            )
                        else:
                            self.solventOA_verline2 = (
                                self.solventOA_figure.axes.axvline(
                                    x, color="orange", linewidth=2
                                )
                            )
                    else:
                        self.solventOA_cursorPositions.clear()
                        self.solventOA_cursorPositions.append((x, y))
                        self.solventOA_verline1.set_xdata(x)
                        self.solventOA_verline2.set_xdata(None)

                    self.solventOA_figure.draw_idle()

            case "Sample":
                pass

    def draw_fitting_line(self, ftype: str, stype: str) -> None:
        match ftype:
            case "Silica":
                if not self.silicaCA_fittingLine_drawn:
                    (self.silica_fitting_line_ca,) = (
                        self.silicaCA_figure.axes.plot(
                            self.silica_data_set[0], self.result, "r"
                        )
                    )
                    self.silicaCA_figure.draw_idle()

                    self.silicaCA_fittingLine_drawn = True
                else:
                    if hasattr(self, "silica_data_set"):
                        self.silica_fitting_line_ca.set_xdata(
                            self.silica_data_set[0]
                        )
                        self.silica_fitting_line_ca.set_ydata(self.result)

                        self.silicaCA_figure.axes.set_xlim(
                            left=-self.z_range / 2 * 1000,
                            right=self.z_range / 2 * 1000,
                        )  # displayed in mm
                        self.silicaCA_figure.axes.relim()
                        self.silicaCA_figure.axes.autoscale_view()
                        self.silicaCA_figure.draw_idle()

            case "Solvent":
                if stype == "CA":
                    if not self.solventCA_fittingLine_drawn:
                        (self.solvent_fitting_line_ca,) = (
                            self.solventCA_figure.axes.plot(
                                self.solvent_data_set[0], self.result, "r"
                            )
                        )
                        self.solventCA_figure.draw_idle()

                        self.solventCA_fittingLine_drawn = True
                    else:
                        if hasattr(self, "solvent_data_set"):
                            self.solvent_fitting_line_ca.set_xdata(
                                self.solvent_data_set[0]
                            )
                            self.solvent_fitting_line_ca.set_ydata(self.result)

                            self.solventCA_figure.axes.set_xlim(
                                left=-self.z_range / 2 * 1000,
                                right=self.z_range / 2 * 1000,
                            )  # displayed in mm
                            self.solventCA_figure.axes.relim()
                            self.solventCA_figure.axes.autoscale_view()
                            self.solventCA_figure.draw_idle()

                elif stype == "OA":
                    if not self.solventOA_fittingLine_drawn:
                        (self.solvent_fitting_line_oa,) = (
                            self.solventOA_figure.axes.plot(
                                self.solvent_data_set[0], self.result, "r"
                            )
                        )
                        self.solventOA_figure.draw_idle()

                        self.solventOA_fittingLine_drawn = True
                    else:
                        if hasattr(self, "solvent_data_set"):
                            self.solvent_fitting_line_oa.set_xdata(
                                self.solvent_data_set[0]
                            )
                            self.solvent_fitting_line_oa.set_ydata(self.result)

                            self.solventOA_figure.axes.set_xlim(
                                left=-self.z_range / 2 * 1000,
                                right=self.z_range / 2 * 1000,
                            )  # displayed in mm
                            self.solventOA_figure.axes.relim()
                            self.solventOA_figure.axes.autoscale_view()
                            self.solventOA_figure.draw_idle()

            case "Sample":
                if stype == "CA":
                    if not self.sampleCA_fittingLine_drawn:
                        (self.sample_fitting_line_ca,) = (
                            self.sampleCA_figure.axes.plot(
                                self.sample_data_set[0], self.result, "r"
                            )
                        )
                        self.sampleCA_figure.draw_idle()

                        self.sampleCA_fittingLine_drawn = True
                    else:
                        if hasattr(self, "sample_data_set"):
                            self.sample_fitting_line_ca.set_xdata(
                                self.sample_data_set[0]
                            )
                            self.sample_fitting_line_ca.set_ydata(self.result)

                            self.sampleCA_figure.axes.set_xlim(
                                left=-self.z_range / 2 * 1000,
                                right=self.z_range / 2 * 1000,
                            )  # displayed in mm
                            self.sampleCA_figure.axes.relim()
                            self.sampleCA_figure.axes.autoscale_view()
                            self.sampleCA_figure.draw_idle()

                elif stype == "OA":
                    if not self.sampleOA_fittingLine_drawn:
                        (self.sample_fitting_line_oa,) = (
                            self.sampleOA_figure.axes.plot(
                                self.sample_data_set[0], self.result, "r"
                            )
                        )
                        self.sampleOA_figure.draw_idle()

                        self.sampleOA_fittingLine_drawn = True
                    else:
                        if hasattr(self, "sample_data_set"):
                            self.sample_fitting_line_oa.set_xdata(
                                self.sample_data_set[0]
                            )
                            self.sample_fitting_line_oa.set_ydata(self.result)

                            self.sampleOA_figure.axes.set_xlim(
                                left=-self.z_range / 2 * 1000,
                                right=self.z_range / 2 * 1000,
                            )  # displayed in mm
                            self.sampleOA_figure.axes.relim()
                            self.sampleOA_figure.axes.autoscale_view()
                            self.sampleOA_figure.draw_idle()

    def fit_automatically(self, ftype: str, stype: str) -> None:
        """Fit using optimization. Delegates to workflows.

        CRITICAL: Extract initial parameters from data geometry BEFORE optimization.
        This gives the optimizer excellent initial guesses.
        """

        # Precondition: require silica to be fitted first
        if ftype != "Silica" and not getattr(
            self, "silica_autofit_done", False
        ):
            self.showdialog("Info", "Fit silica first.")
            return

        # Precondition: require solvent to be fitted first (for sample)
        if ftype == "Sample" and not getattr(
            self, "solventCA_autofit_done", False
        ):
            self.showdialog("Info", "Fit solvent first.")
            return

        self.get_general_parameters()

        # CRITICAL: Extract initial parameters from data geometry
        # on_data_load=True forces params_from_geometry() to be called
        self.get_curve_interpretation(
            ftype, stype, "from_geometry", on_data_load=True
        )

        try:
            match (ftype, stype):
                case ("Silica", "CA"):
                    result = self._fit_automatically_silica_ca()
                case ("Solvent", "CA"):
                    result = self._fit_automatically_solvent_ca()
                case ("Solvent", "OA"):
                    result = self._fit_automatically_solvent_oa()
                case ("Sample", "CA"):
                    result = self._fit_automatically_sample_ca()
                case ("Sample", "OA"):
                    result = self._fit_automatically_sample_oa()
                case _:
                    self.showdialog(
                        "Error", f"Unsupported fitting type: {ftype} {stype}"
                    )
                    return

            if result and result.converged:
                self.result = result.curve
                self.draw_fitting_line(ftype, stype)
                # Extract results from fitting
                self.get_curve_interpretation(ftype, stype, "from_autofit")
                self.set_sliders_positions(ftype, stype)
                self.set_fit_summary(ftype, stype, caller="auto")
                # Update autofit flags
                match ftype:
                    case "Silica":
                        self.silica_autofit_done = True
                    case "Solvent":
                        self.solventCA_autofit_done = True
                    case "Sample":
                        match stype:
                            case "CA":
                                self.sampleCA_autofit_done = True
                            case "OA":
                                self.sampleOA_autofit_done = True
            elif result:
                self.showdialog("Warning", result.message)
            else:
                self.showdialog(
                    "Error", f"Fitting returned None for {ftype} {stype}"
                )

        except Exception as e:
            import logging

            logging.error(
                f"fit_automatically failed for {ftype} {stype}: {e}",
                exc_info=True,
            )
            self.showdialog("Error", f"Automatic fitting failed: {str(e)}")

    def _fit_automatically_silica_ca(self):
        """Silica CA automatic fitting via workflow."""
        from workflows.silica_workflow import SilicaCAFittingParams

        line = self.silicaCA_figure.axes.get_lines()[0]

        params = SilicaCAFittingParams(
            dphi0=self.silicaCA_DPhi0,
            beamwaist=self.silicaCA_beamwaist,
            zero_level=self.silicaCA_zeroLevel,
            center_point=self.silicaCA_centerPoint,
            z_positions=self.silica_data_set[0],
            ca_data=line.get_ydata(),
            n2=self.silica_n2,
            d0=self.d0,
            ra=self.ra,
            wavelength=self.lda,
            z_range=self.z_range,
        )

        result = self.silica_ca_workflow.fit_automatic(params)

        if result.converged:
            # Store minimizer result for error extraction in get_curve_interpretation
            self.silicaCA_minimizerResult = result.minimizer_result

        return result

    def _fit_automatically_solvent_ca(self):
        """Solvent CA automatic fitting via workflow."""
        from workflows.solvent_workflow import SolventCAFittingParams

        line = self.solventCA_figure.axes.get_lines()[0]

        params = SolventCAFittingParams(
            dphi0=self.solventCA_DPhi0,
            beamwaist=self.solventCA_beamwaist,
            zero_level=self.solventCA_zeroLevel,
            center_point=self.solventCA_centerPoint,
            z_positions=self.solvent_data_set[0],
            ca_data=line.get_ydata(),
            n2=getattr(self, "solvent_n2", 0),
            d0=self.d0,
            ra=self.ra,
            wavelength=self.lda,
            z_range=self.z_range,
        )

        result = self.solvent_ca_workflow.fit_automatic(params)

        if result.converged:
            self.solventCA_minimizerResult = result.minimizer_result
            self.solventCA_DPhi0 = result.params["dphi0"]
            self.solventCA_beamwaist = result.params["beamwaist"]
            self.solventCA_centerPoint = result.params["center"]

            # Sync OA center if custom is disabled
            try:
                if (
                    hasattr(self, "solventOA_customCenterPoint_checkBox")
                    and not self.solventOA_customCenterPoint_checkBox.isChecked()
                ):
                    self.solventOA_centerPoint = self.solventCA_centerPoint
                    self.solventOA_centerPoint_doubleSpinBox.blockSignals(True)
                    self.solventOA_centerPoint_doubleSpinBox.setValue(
                        self.solventOA_centerPoint
                    )
                    self.solventOA_centerPoint_doubleSpinBox.blockSignals(False)

                    mapped = int(
                        round(
                            self.solventOA_centerPoint
                            + SOLVENT_OA_SLIDER_CENTER_OFFSET
                        )
                    )
                    self.solventOA_centerPoint_slider.blockSignals(True)
                    self.solventOA_centerPoint_slider.setValue(mapped)
                    self.solventOA_centerPoint_slider.blockSignals(False)
            except Exception:
                pass

        return result

    def _fit_automatically_solvent_oa(self):
        """Solvent OA automatic fitting via workflow."""
        from workflows.solvent_workflow import SolventOAFittingParams

        # Check if should skip
        if self.solventOA_isAbsorption_checkBox.isChecked():
            line = self.solventOA_figure.axes.get_lines()[0]
            line.set_ydata(np.ones_like(line.get_xdata()))
            self.solventOA_figure.draw_idle()
            self.solventOA_fittingLine_drawn = False
            return None

        line = self.solventOA_figure.axes.get_lines()[0]

        # Determine center point to use
        try:
            if (
                hasattr(self, "solventOA_customCenterPoint_checkBox")
                and self.solventOA_customCenterPoint_checkBox.isChecked()
            ):
                center = float(self.solventOA_centerPoint_doubleSpinBox.value())
            else:
                center = float(self.solventCA_centerPoint)
        except Exception:
            center = self.solventCA_centerPoint

        self.solventOA_centerPoint = center

        params = SolventOAFittingParams(
            t_value=getattr(self, "solventOA_T", 0.5),
            beamwaist=self.solventCA_beamwaist,
            zero_level=self.solventOA_zeroLevel,
            center_point=center,
            z_positions=self.solvent_data_set[0],
            oa_data=line.get_ydata(),
            n2=getattr(self, "solvent_n2", 0),
            d0=self.d0,
            ra=self.ra,
            wavelength=self.lda,
            z_range=self.z_range,
            no_absorption=False,
        )

        result = self.solvent_oa_workflow.fit_automatic(params)

        if result.converged:
            self.solventOA_minimizerResult = result.minimizer_result
            self.solventOA_T = result.params.get("t", self.solventOA_T)

        return result

    def _fit_automatically_sample_ca(self):
        """Sample CA automatic fitting via workflow."""
        from workflows.sample_workflow import SampleCAFittingParams

        line = self.sampleCA_figure.axes.get_lines()[0]

        params = SampleCAFittingParams(
            dphi0=self.sampleCA_DPhi0,
            beamwaist=self.sampleCA_beamwaist,
            zero_level=self.sampleCA_zeroLevel,
            center_point=self.sampleCA_centerPoint,
            z_positions=self.sample_data_set[0],
            ca_data=line.get_ydata(),
            n2=getattr(self, "sample_n2", 1e-15),
            d0=self.d0,
            ra=self.ra,
            wavelength=self.lda,
            z_range=self.z_range,
        )

        result = self.sample_ca_workflow.fit_automatic(params)

        if result.converged:
            self.sampleCA_minimizerResult = result.minimizer_result

        return result

    def _fit_automatically_sample_oa(self):
        """Sample OA automatic fitting via workflow."""
        from workflows.sample_workflow import SampleOAFittingParams

        line = self.sampleOA_figure.axes.get_lines()[0]

        params = SampleOAFittingParams(
            transmittance=getattr(self, "sampleOA_transmittance", 0.5),
            beamwaist=self.sampleCA_beamwaist,
            zero_level=self.sampleOA_zeroLevel,
            center_point=self.sampleOA_centerPoint,
            z_positions=self.sample_data_set[0],
            oa_data=line.get_ydata(),
            beta=getattr(self, "sample_beta", 0),
            d0=self.d0,
            ra=self.ra,
            wavelength=self.lda,
            z_range=self.z_range,
        )

        result = self.sample_oa_workflow.fit_automatic(params)

        if result.converged:
            self.sampleOA_minimizerResult = result.minimizer_result

        return result

    def fit_manually(self, ftype: str, stype: str, activated_by=None) -> None:
        """Fit using sliders. Delegates to workflows."""
        logger.debug(
            f"Window.fit_manually: running with params: {ftype=}, {stype=}, {activated_by=}"
        )
        # Precondition check
        if ftype != "Silica" and not getattr(
            self, "silica_autofit_done", False
        ):
            self.showdialog("Info", "Fit silica first.")
            return

        if ftype == "Sample" and not getattr(
            self, "solventCA_autofit_done", False
        ):
            self.showdialog("Info", "Fit solvent first.")
            return

        logger.debug(
            "Window.fit_manually: about to Window.get_general_parameters"
        )
        self.get_general_parameters()
        logger.debug(
            "Window.fit_manually: about to Window.get_curve interpretation"
        )
        self.get_curve_interpretation(ftype, stype, "from_geometry")

        try:
            match (ftype, stype):
                case ("Silica", "CA"):
                    result = self._fit_manually_silica_ca()
                case ("Solvent", "CA"):
                    result = self._fit_manually_solvent_ca()
                case ("Solvent", "OA"):
                    result = self._fit_manually_solvent_oa()
                case ("Sample", "CA"):
                    result = self._fit_manually_sample_ca()
                case ("Sample", "OA"):
                    result = self._fit_manually_sample_oa()
                case _:
                    self.showdialog(
                        "Error", f"Unsupported fitting type: {ftype} {stype}"
                    )
                    return

            if result is None:
                return  # Fitting was skipped

            if result.converged:
                self.result = result.curve
                self.draw_fitting_line(ftype, stype)
            else:
                self.showdialog("Warning", result.message)

        except Exception as e:
            import logging

            logging.error(
                f"fit_manually failed for {ftype} {stype}: {e}", exc_info=True
            )
            self.showdialog("Error", f"Manual fitting failed: {str(e)}")

        self.calculate_derived_parameters(ftype)
        self.set_fit_summary(ftype, stype, caller="manual")

    def _fit_manually_silica_ca(self):
        """Silica CA manual fitting via workflow."""
        from workflows.silica_workflow import SilicaCAFittingParams

        line = self.silicaCA_figure.axes.get_lines()[0]

        params = SilicaCAFittingParams(
            dphi0=self.silicaCA_DPhi0,
            beamwaist=self.silicaCA_beamwaist,
            zero_level=self.silicaCA_zeroLevel,
            center_point=self.silicaCA_centerPoint,
            z_positions=self.silica_data_set[0],
            ca_data=line.get_ydata(),
            n2=self.silica_n2,
            d0=self.d0,
            ra=self.ra,
            wavelength=self.lda,
            z_range=self.z_range,
        )

        return self.silica_ca_workflow.fit_manual(params)

    def _fit_manually_solvent_ca(self):
        """Solvent CA manual fitting via workflow."""
        from workflows.solvent_workflow import SolventCAFittingParams

        line = self.solventCA_figure.axes.get_lines()[0]

        params = SolventCAFittingParams(
            dphi0=self.solventCA_DPhi0,
            beamwaist=self.solventCA_beamwaist,
            zero_level=self.solventCA_zeroLevel,
            center_point=self.solventCA_centerPoint,
            z_positions=self.solvent_data_set[0],
            ca_data=line.get_ydata(),
            n2=getattr(self, "solvent_n2", 0),
            d0=self.d0,
            ra=self.ra,
            wavelength=self.lda,
            z_range=self.z_range,
        )

        return self.solvent_ca_workflow.fit_manual(params)

    def _fit_manually_solvent_oa(self):
        """Solvent OA manual fitting via workflow."""
        from workflows.solvent_workflow import SolventOAFittingParams

        # Skip if no absorption assumed
        if self.solventOA_isAbsorption_checkBox.isChecked():
            return None

        line = self.solventOA_figure.axes.get_lines()[0]

        # Determine center point to use
        try:
            if (
                hasattr(self, "solventOA_customCenterPoint_checkBox")
                and self.solventOA_customCenterPoint_checkBox.isChecked()
            ):
                center = float(self.solventOA_centerPoint_doubleSpinBox.value())
            else:
                center = float(self.solventCA_centerPoint)
        except Exception:
            center = self.solventCA_centerPoint

        self.solventOA_centerPoint = center

        params = SolventOAFittingParams(
            t_value=getattr(self, "solventOA_T", 0.5),
            beamwaist=self.solventCA_beamwaist,
            zero_level=self.solventOA_zeroLevel,
            center_point=center,
            z_positions=self.solvent_data_set[0],
            oa_data=line.get_ydata(),
            n2=getattr(self, "solvent_n2", 0),
            d0=self.d0,
            ra=self.ra,
            wavelength=self.lda,
            z_range=self.z_range,
            no_absorption=self.solventOA_isAbsorption_checkBox.isChecked(),
        )

        return self.solvent_oa_workflow.fit_manual(params)

    def _fit_manually_sample_ca(self):
        """Sample CA manual fitting via workflow."""
        from workflows.sample_workflow import SampleCAFittingParams

        line = self.sampleCA_figure.axes.get_lines()[0]

        params = SampleCAFittingParams(
            dphi0=self.sampleCA_DPhi0,
            beamwaist=self.sampleCA_beamwaist,
            zero_level=self.sampleCA_zeroLevel,
            center_point=self.sampleCA_centerPoint,
            z_positions=self.sample_data_set[0],
            ca_data=line.get_ydata(),
            n2=getattr(self, "sample_n2", 1e-15),
            d0=self.d0,
            ra=self.ra,
            wavelength=self.lda,
            z_range=self.z_range,
        )

        return self.sample_ca_workflow.fit_manual(params)

    def _fit_manually_sample_oa(self):
        """Sample OA manual fitting via workflow."""
        from workflows.sample_workflow import SampleOAFittingParams

        line = self.sampleOA_figure.axes.get_lines()[0]

        params = SampleOAFittingParams(
            transmittance=getattr(self, "sampleOA_transmittance", 0.5),
            beamwaist=self.sampleCA_beamwaist,
            zero_level=self.sampleOA_zeroLevel,
            center_point=self.sampleOA_centerPoint,
            z_positions=self.sample_data_set[0],
            oa_data=line.get_ydata(),
            beta=getattr(self, "sample_beta", 0),
            d0=self.d0,
            ra=self.ra,
            wavelength=self.lda,
            z_range=self.z_range,
        )

        return self.sample_oa_workflow.fit_manual(params)

    def _fit_automatically_legacy(self, ftype: str, stype: str) -> None:
        self.get_general_parameters()
        self.get_curve_interpretation(ftype, stype, "from_geometry")

        match ftype:
            case "Silica":
                if stype == "CA":
                    # Datapoints for the curve to be fitted to
                    line = self.silicaCA_figure.axes.get_lines()[0]
                else:
                    return

                line_data = line.get_data()

                self.silica_calculation = Fitting(
                    self.silica_curves,
                    self.silicaCA_DPhi0,
                    self.silicaCA_beamwaist,
                    self.silicaCA_zeroLevel,
                    self.silicaCA_centerPoint,
                    self.silica_nop,
                    line_data[1],
                )
                minimizer_result, self.result = (
                    self.silica_calculation.automatic(
                        self.z_range, ftype, stype, line_data, window=self
                    )
                )
                self.silicaCA_minimizerResult = minimizer_result
                self.silica_autofit_done = True
                self.draw_fitting_line(ftype, stype)
                self.get_curve_interpretation(ftype, stype, "from_autofit")
                self.set_fit_summary(ftype, stype, caller="auto")
                # Use these exact number from error-corrected fit parameters to set sliders to their positions
                self.set_sliders_positions(ftype, stype)

            case "Solvent":
                if not self.silica_autofit_done:
                    self.showdialog("Info", "Fit silica first.")
                    return

                # If user chose 'Assume no absorption', skip OA fitting
                if stype == "OA":
                    try:
                        if self.solventOA_isAbsorption_checkBox.isChecked():
                            # Set OA curve to flat line (T = 1)
                            line = self.solventOA_figure.axes.get_lines()[0]
                            line.set_ydata(np.ones_like(line.get_xdata()))
                            self.solventOA_figure.draw_idle()
                            return
                    except Exception:
                        pass

                if stype == "CA":
                    # Data to be fitted
                    line = self.solventCA_figure.axes.get_lines()[0]
                elif stype == "OA":
                    line = self.solventOA_figure.axes.get_lines()[0]

                line_data = line.get_data()

                nop = len(self.solvent_data_set[0])
                self.solvent_data_set[0] = np.array(
                    [
                        (self.z_range * zz / nop - self.z_range / 2) * 1000
                        for zz in range(nop)
                    ]
                )  # [mm] update positions with newly-read Z-scan range value

                # Build Integration object appropriate for OA/CA automatic fitting
                if stype == "CA":
                    self.solvent_curves = Integration(
                        beta=0,
                        n2=getattr(self, "solvent_n2", 0),
                        DPhi0=self.solventCA_DPhi0,
                        positions=self.solvent_data_set[0],
                        d0=self.d0,
                        aperture_radius=self.ra,
                        wavelength=self.lda,
                        beamwaist=self.solventCA_beamwaist,
                        n_components=N_COMPONENTS,
                        integration_steps=INTEGRATION_STEPS,
                        stype="CA",
                    )
                    amplitude_init = self.solventCA_DPhi0
                else:  # OA
                    amplitude_init = as_float(
                        getattr(self, "solventOA_T", None)
                    )
                    if amplitude_init is None and hasattr(
                        self, "solventOA_T_doubleSpinBox"
                    ):
                        amplitude_init = self.solventOA_T_doubleSpinBox.value()
                    if amplitude_init is None:
                        amplitude_init = 0
                    self.solvent_curves = Integration(
                        beta=amplitude_init,
                        n2=getattr(self, "solvent_n2", 0),
                        DPhi0=0,
                        positions=self.solvent_data_set[0],
                        d0=self.d0,
                        aperture_radius=self.ra,
                        wavelength=self.lda,
                        beamwaist=self.solventCA_beamwaist,
                        n_components=N_COMPONENTS,
                        integration_steps=INTEGRATION_STEPS,
                        stype="OA",
                    )

                # Determine initial center for Fitting
                if stype == "OA":
                    # For OA: if custom center is enabled, read from spinbox; otherwise use CA center
                    try:
                        if (
                            hasattr(
                                self, "solventOA_customCenterPoint_checkBox"
                            )
                            and self.solventOA_customCenterPoint_checkBox.isChecked()
                        ):
                            # Custom enabled: read directly from spinbox to ensure latest value
                            initial_center = float(
                                self.solventOA_centerPoint_doubleSpinBox.value()
                            )
                        else:
                            # Custom disabled: use CA center
                            initial_center = float(self.solventCA_centerPoint)
                    except Exception:
                        initial_center = self.solventCA_centerPoint
                else:
                    initial_center = self.solventCA_centerPoint

                # Update internal state to match what will be fitted
                self.solventOA_centerPoint = initial_center

                self.solvent_calculation = Fitting(
                    self.solvent_curves,
                    ensure_float(amplitude_init),
                    ensure_float(self.solventCA_beamwaist),
                    ensure_float(self.solventCA_zeroLevel),
                    ensure_float(initial_center),
                    nop,
                    line_data[1],
                )

                try:
                    # Determine whether center should be varied based on 'Custom' checkbox
                    # If custom center is NOT checked, do not vary the center (use CA center)
                    vary_center = False
                    try:
                        vary_center = self.solventOA_customCenterPoint_checkBox.isChecked()
                    except Exception:
                        pass

                    # Decide whether beamwaist should be allowed to vary:
                    # - For CA: vary only if user enabled custom beamwaist
                    # - For OA: keep beamwaist fixed (OA fit uses CA beamwaist)
                    try:
                        if stype == "CA":
                            vary_bw = bool(
                                getattr(
                                    self,
                                    "solventCA_customBeamwaist_checkBox",
                                    None,
                                )
                                and self.solventCA_customBeamwaist_checkBox.isChecked()
                            )
                        else:
                            vary_bw = False
                    except Exception:
                        vary_bw = False

                    minimizer_result, self.result = (
                        self.solvent_calculation.automatic(
                            self.z_range,
                            ftype,
                            stype,
                            line_data,
                            window=self,
                            vary_beamwaist=vary_bw,
                            vary_centerpoint=vary_center,
                        )
                    )
                except Exception as e:
                    logging.error(f"Solvent automatic fitting error: {e}")
                    self.showdialog(
                        "Error", f"Solvent automatic fitting failed: {str(e)}"
                    )
                    return

                self.draw_fitting_line(ftype, stype)

                if stype == "CA":
                    self.solventCA_minimizerResult = minimizer_result
                    self.solventCA_DPhi0 = as_float(
                        minimizer_result.params["DPhi0"].value
                    )
                    self.solventCA_beamwaist = as_float(
                        minimizer_result.params["Beamwaist"].value
                    )
                    # UPDATE centerPoint from CA fit result (critical!)
                    self.solventCA_centerPoint = minimizer_result.params[
                        "Center"
                    ].value
                    # If OA custom center is disabled, update OA center to CA center now
                    try:
                        if (
                            hasattr(
                                self, "solventOA_customCenterPoint_checkBox"
                            )
                            and not self.solventOA_customCenterPoint_checkBox.isChecked()
                        ):
                            self.solventOA_centerPoint = (
                                self.solventCA_centerPoint
                            )
                            try:
                                self.solventOA_centerPoint_doubleSpinBox.blockSignals(
                                    True
                                )
                                self.solventOA_centerPoint_doubleSpinBox.setValue(
                                    self.solventOA_centerPoint
                                )
                            finally:
                                try:
                                    self.solventOA_centerPoint_doubleSpinBox.blockSignals(
                                        False
                                    )
                                except Exception:
                                    pass
                            try:
                                mapped = int(
                                    round(
                                        self.solventOA_centerPoint
                                        + SOLVENT_OA_SLIDER_CENTER_OFFSET
                                    )
                                )
                                self.solventOA_centerPoint_slider.blockSignals(
                                    True
                                )
                                self.solventOA_centerPoint_slider.setValue(
                                    mapped
                                )
                            finally:
                                try:
                                    self.solventOA_centerPoint_slider.blockSignals(
                                        False
                                    )
                                except Exception:
                                    pass
                    except Exception:
                        pass
                elif stype == "OA":
                    self.solventOA_minimizerResult = minimizer_result
                    # Extract fitted T and store for slider/spinbox mapping
                    try:
                        if "T" in minimizer_result.params:
                            self.solventOA_T = minimizer_result.params[
                                "T"
                            ].value
                        else:
                            # fallback: check for any param that looks like T
                            for k, v in minimizer_result.params.items():
                                if k.lower() == "t":
                                    self.solventOA_T = v.value
                                    break
                    except Exception:
                        pass

                self.calculate_derived_parameters(ftype)
                # Extract curve interpretation from autofit result
                self.get_curve_interpretation(ftype, stype, "from_autofit")

                # Use these exact number from error-corrected fit parameters to set sliders to their positions and display these values
                self.set_sliders_positions(ftype, stype)

                self.set_fit_summary(ftype, stype, caller="auto")

                self.solventCA_autofit_done = True

            case "Sample":
                if not self.silica_autofit_done:
                    self.showdialog("Info", "Fit silica first.")

                    if not self.solventCA_autofit_done:
                        self.showdialog("Info", "Fit solvent first.")
                        return
                    else:
                        return

                # Initialize sample fitter
                sample_fitter = SampleFitter()

                # Proceed with sample fitting
                if stype == "CA":
                    line = self.sampleCA_figure.axes.get_lines()[0]
                    line_data = line.get_data()

                    nop = len(self.sample_data_set[0])
                    self.sample_data_set[0] = np.array(
                        [
                            (self.z_range * zz / nop - self.z_range / 2) * 1000
                            for zz in range(nop)
                        ]
                    )

                    self.sample_curves = Integration(
                        beta=0,
                        n2=self.sample_n2
                        if hasattr(self, "sample_n2")
                        else 1e-15,
                        DPhi0=self.sampleCA_DPhi0,
                        positions=self.sample_data_set[0],
                        d0=self.d0,
                        aperture_radius=self.ra,
                        wavelength=self.lda,
                        beamwaist=self.sampleCA_beamwaist,
                        n_components=N_COMPONENTS,
                        integration_steps=INTEGRATION_STEPS,
                        stype="CA",
                    )

                    try:
                        minimizer_result, self.result = (
                            sample_fitter.fit_automatically_ca(
                                sample_curves=self.sample_curves,
                                initial_dphi0=self.sampleCA_DPhi0,
                                initial_beamwaist=ensure_float(
                                    self.sampleCA_beamwaist
                                ),
                                initial_zero_level=ensure_float(
                                    self.sampleCA_zeroLevel
                                ),
                                initial_centerpoint=ensure_float(
                                    self.sampleCA_centerPoint
                                ),
                                nop=nop,
                                z_range=self.z_range,
                                line_data=line_data,
                                window=self,
                            )
                        )

                        self.sampleCA_minimizerResult = minimizer_result
                        self.draw_fitting_line(ftype, stype)
                        self.sampleCA_fittingDone = True
                        self.get_curve_interpretation(
                            ftype, stype, "from_autofit"
                        )
                        self.set_fit_summary(ftype, stype, caller="auto")
                        self.set_sliders_positions(ftype, stype)

                    except Exception as e:
                        logging.error(f"Sample CA fitting error: {str(e)}")
                        self.showdialog(
                            "Error", f"Sample CA fitting failed: {str(e)}"
                        )
                        return

                elif stype == "OA":
                    line = self.sampleOA_figure.axes.get_lines()[0]
                    line_data = line.get_data()
                    nop = len(self.sample_data_set[0])

                    self.sample_curves_oa = Integration(
                        beta=self.sample_beta
                        if hasattr(self, "sample_beta")
                        else 0,
                        n2=1e-15,
                        DPhi0=0,
                        positions=self.sample_data_set[0],
                        d0=self.d0,
                        aperture_radius=self.ra,
                        wavelength=self.lda,
                        beamwaist=self.sampleCA_beamwaist,
                        n_components=N_COMPONENTS,
                        integration_steps=INTEGRATION_STEPS,
                        stype="OA",
                    )

                    try:
                        minimizer_result, self.result = (
                            sample_fitter.fit_automatically_oa(
                                sample_curves=self.sample_curves_oa,
                                initial_transmittance=0.5,
                                initial_beamwaist=self.sampleCA_beamwaist,
                                initial_zero_level=self.sampleOA_zeroLevel,
                                initial_centerpoint=self.sampleOA_centerPoint,
                                nop=nop,
                                z_range=self.z_range,
                                line_data=line_data,
                                window=self,
                            )
                        )

                        self.sampleOA_minimizerResult = minimizer_result
                        self.draw_fitting_line(ftype, stype)
                        self.sampleOA_fittingDone = True
                        self.get_curve_interpretation(
                            ftype, stype, "from_autofit"
                        )
                        self.set_fit_summary(ftype, stype, caller="auto")
                        self.set_sliders_positions(ftype, stype)

                    except Exception as e:
                        logging.error(f"Sample OA fitting error: {str(e)}")
                        self.showdialog(
                            "Error", f"Sample OA fitting failed: {str(e)}"
                        )
                        return

    # def fit_manually(self, ftype: str, stype: str) -> None:
    #     """Fit using sliders. Delegates to workflows."""

    #     # Precondition check
    #     if ftype != "Silica" and not getattr(self, 'silica_autofit_done', False):
    #         self.showdialog('Info', 'Fit silica first.')
    #         return

    #     self.get_general_parameters()
    #     self.get_curve_interpretation(ftype, stype, 'from_geometry')

    #     try:
    #         if ftype == "Solvent" and stype == "OA":
    #             # Use workflow for OA
    #             from workflows.solvent_workflow import SolventOAFittingParams

    #             params = SolventOAFittingParams(
    #                 t_value=self.solventOA_T,
    #                 beamwaist=self.solventCA_beamwaist,
    #                 zero_level=self.solventOA_zeroLevel,
    #                 center_point=self.solventOA_centerPoint,
    #                 z_positions=self.solvent_data_set[0],
    #                 oa_data=self.solvent_data_set[3],
    #                 n2=getattr(self, 'solvent_n2', 0),
    #                 d0=self.d0,
    #                 ra=self.ra,
    #                 wavelength=self.lda,
    #                 z_range=self.z_range,
    #                 no_absorption=self.solventOA_isAbsorption_checkBox.isChecked(),
    #             )

    #             result = self.solvent_oa_workflow.fit_manual(params)

    #             if result is None:
    #                 # Fitting was skipped
    #                 self.solventOA_fittingLine_drawn = False
    #                 return

    #             if result.converged:
    #                 self.result = result.curve
    #                 self.draw_fitting_line(ftype, stype)
    #             else:
    #                 self.showdialog('Warning', result.message)

    #         elif ftype == "Solvent" and stype == "CA":
    #             # Use workflow for CA
    #             from workflows.solvent_workflow import SolventCAFittingParams

    #             params = SolventCAFittingParams(
    #                 dphi0=self.solventCA_DPhi0,
    #                 beamwaist=self.solventCA_beamwaist,
    #                 zero_level=self.solventCA_zeroLevel,
    #                 center_point=self.solventCA_centerPoint,
    #                 z_positions=self.solvent_data_set[0],
    #                 ca_data=self.solvent_data_set[1],
    #                 n2=getattr(self, 'solvent_n2', 0),
    #                 d0=self.d0,
    #                 ra=self.ra,
    #                 wavelength=self.lda,
    #                 z_range=self.z_range,
    #             )

    #             result = self.solvent_ca_workflow.fit_manual(params)

    #             if result.converged:
    #                 self.result = result.curve
    #                 self.draw_fitting_line(ftype, stype)
    #                 self.solventCA_autofit_done = False
    #             else:
    #                 self.showdialog('Warning', result.message)

    #         else:
    #             # Keep existing logic for Silica and Sample (for now)
    #             print("Switching to legacy")
    #             self._fit_manually_legacy(ftype, stype)

    #     except Exception as e:
    #         import logging
    #         logging.error(f"fit_manually failed: {e}", exc_info=True)
    #         self.showdialog('Error', f'Fitting failed: {str(e)}')

    def _fit_manually_legacy(self, ftype: str, stype: str) -> None:
        """Keep old Silica/Sample logic here temporarily."""
        # Move existing fit_manually code here for Silica and Sample
        # This is temporary - you'll refactor these later

        """
        Perform a manual fit using the current GUI slider or spinbox values.

        This method is triggered whenever the user interacts with any fitting-related
        control (sliders, spinboxes, or programmatic updates). It reads the current
        parameter values from the UI, recomputes the Z-scan curve for the selected
        sample and aperture type, and updates the corresponding fitting plot.

        Args:
            ftype:
                The sample category being fitted. One of:
                    - "Silica"
                    - "Solvent"
                    - "Sample"

            stype:
                The aperture mode:
                    - "CA" for closed-aperture fitting
                    - "OA" for open-aperture fitting

            activated_by:
                Optional string describing what triggered the update.
                Typical values include:
                    - "slider" (user moved a slider)
                    - "spinbox" (user edited a spinbox)
                    - None (internal or programmatic update)
                This is used only for UI logic and does not affect the physics.

        Notes:
            • This function does not perform optimization; it simply evaluates the
            model with the current UI parameters.

            • It is called frequently (e.g., on every slider move), so it must remain
            lightweight.

            • All required fitting parameters (center point, zero level, DPhi0/T,
            beam waist, etc.) are taken directly from the GUI state.
        """

        # Preconditions: require silica to be fitted first for solvent/sample
        if ftype != "Silica" and not getattr(
            self, "silica_autofit_done", False
        ):
            self.showdialog("Info", "Fit silica first.")
            return

        # Retrieve "General parameters"
        self.get_general_parameters()
        self.get_curve_interpretation(ftype, stype, "from_geometry")

        match ftype:
            case "Silica":
                if stype == "CA":
                    # Datapoints to fit the curve to
                    line = self.silicaCA_figure.axes.get_lines()[0]
                    line_data = line.get_ydata()

                    self.silica_curves = Integration(
                        SILICA_BETA,
                        self.silica_n2,
                        self.silicaCA_DPhi0,
                        self.silica_data_set[0],
                        self.d0,
                        self.ra,
                        self.lda,
                        self.silicaCA_beamwaist,
                        N_COMPONENTS,
                        INTEGRATION_STEPS,
                    )
                    self.silica_calculation = Fitting(
                        self.silica_curves,
                        self.silicaCA_DPhi0,
                        self.silicaCA_beamwaist,
                        self.silicaCA_zeroLevel,
                        self.silicaCA_centerPoint,
                        len(self.silica_data_set[0]),
                        line_data,
                    )
                    self.result = self.silica_calculation.manual(
                        self.silicaCA_zeroLevel,
                        self.silicaCA_centerPoint,
                        self.silicaCA_DPhi0,
                        self.silicaCA_beamwaist,
                        self.z_range,
                        window=self,
                        stype=stype,
                    )
                    self.draw_fitting_line(ftype, stype)
                    self.silica_autofit_done = False
                    self.get_curve_interpretation(ftype, stype, "from_geometry")

            case "Solvent":
                # Data to be fitted
                if stype == "CA":
                    line = self.solventCA_figure.axes.get_lines()[0]
                    line_data = line.get_ydata()

                    self.solvent_curves = Integration(
                        0,
                        self.solvent_n2,
                        self.solventCA_DPhi0,
                        self.solvent_data_set[0],
                        self.d0,
                        self.ra,
                        self.lda,
                        self.solventCA_beamwaist,
                        N_COMPONENTS,
                        INTEGRATION_STEPS,
                    )  # solvent_beta = 0
                    self.solvent_calculation = Fitting(
                        self.solvent_curves,
                        self.solventCA_DPhi0,
                        self.solventCA_beamwaist,
                        self.solventCA_zeroLevel,
                        self.solventCA_centerPoint,
                        len(self.solvent_data_set[0]),
                        line_data,
                    )
                    self.result = self.solvent_calculation.manual(
                        self.solventCA_zeroLevel,
                        self.solventCA_centerPoint,
                        self.solventCA_DPhi0,
                        self.solventCA_beamwaist,
                        self.z_range,
                        window=self,
                        stype=stype,
                    )
                    self.draw_fitting_line(ftype, stype)

                    self.solventCA_autofit_done = False

                elif stype == "OA":
                    # Manual OA fitting: map slider/spinbox into physical T and run fitting
                    try:
                        # Respect 'Assume no absorption' setting
                        try:
                            if self.solventOA_isAbsorption_checkBox.isChecked():
                                self.showdialog(
                                    "Info",
                                    "Assume no absorption selected — OA fitting skipped.",
                                )
                                return
                        except Exception:
                            pass

                        line = self.solventOA_figure.axes.get_lines()[0]
                        line_data = line.get_ydata()

                        # Ensure we have a physical T value (prefer self.solventOA_T)
                        t_val = getattr(self, "solventOA_T", None)
                        if t_val is None:
                            # fallback to spinbox
                            try:
                                t_val = float(
                                    self.solventOA_T_doubleSpinBox.value()
                                )
                            except Exception:
                                t_val = 0.0

                        self.solvent_curves = Integration(
                            beta=t_val,
                            n2=getattr(self, "solvent_n2", 0),
                            DPhi0=0,
                            positions=self.solvent_data_set[0],
                            d0=self.d0,
                            aperture_radius=self.ra,
                            wavelength=self.lda,
                            beamwaist=self.solventCA_beamwaist,
                            n_components=N_COMPONENTS,
                            integration_steps=INTEGRATION_STEPS,
                            stype="OA",
                        )

                        self.solvent_calculation = Fitting(
                            self.solvent_curves,
                            t_val,
                            self.solventCA_beamwaist,
                            self.solventOA_zeroLevel,
                            self.solventOA_centerPoint,
                            len(self.solvent_data_set[0]),
                            line_data,
                        )

                        self.result = self.solvent_calculation.manual(
                            zero_level=self.solventOA_zeroLevel,
                            centerpoint=self.solventOA_centerPoint,
                            amplitude=t_val,
                            beamwaist=self.solventCA_beamwaist,
                            z_range=self.z_range,
                            window=self,
                            stype=stype,
                        )

                        self.draw_fitting_line(ftype, stype)
                        self.solventOA_fittingDone = False

                    except Exception as e:
                        logging.error(
                            f"Solvent OA manual fitting error: {str(e)}"
                        )
                        return

            case "Sample":
                sample_fitter = SampleFitter()

                if stype == "CA":
                    line = self.sampleCA_figure.axes.get_lines()[0]
                    line_data = line.get_ydata()

                    if not hasattr(self, "sample_curves"):
                        self.sample_curves = Integration(
                            beta=0,
                            n2=self.sample_n2
                            if hasattr(self, "sample_n2")
                            else 1e-15,
                            DPhi0=self.sampleCA_DPhi0,
                            positions=self.sample_data_set[0],
                            d0=self.d0,
                            aperture_radius=self.ra,
                            wavelength=self.lda,
                            beamwaist=self.sampleCA_beamwaist,
                            n_components=N_COMPONENTS,
                            integration_steps=INTEGRATION_STEPS,
                        )

                    self.sample_calculation = Fitting(
                        self.sample_curves,
                        self.sampleCA_DPhi0,
                        self.sampleCA_beamwaist,
                        self.sampleCA_zeroLevel,
                        self.sampleCA_centerPoint,
                        len(self.sample_data_set[0]),
                        line_data,
                    )

                    self.result = self.sample_calculation.manual(
                        zero_level=self.sampleCA_zeroLevel,
                        centerpoint=self.sampleCA_centerPoint,
                        amplitude=self.sampleCA_DPhi0,
                        beamwaist=self.sampleCA_beamwaist,
                        z_range=self.z_range,
                        window=self,
                        stype="CA",
                    )
                    self.draw_fitting_line(ftype, stype)
                    self.sampleCA_fittingDone = False

                elif stype == "OA":
                    line = self.sampleOA_figure.axes.get_lines()[0]
                    line_data = line.get_ydata()

                    if not hasattr(self, "sample_curves_oa"):
                        self.sample_curves_oa = Integration(
                            beta=self.sample_beta
                            if hasattr(self, "sample_beta")
                            else 0,
                            n2=1e-15,
                            DPhi0=0,
                            positions=self.sample_data_set[0],
                            d0=self.d0,
                            aperture_radius=self.ra,
                            wavelength=self.lda,
                            beamwaist=self.sampleCA_beamwaist,
                            n_components=N_COMPONENTS,
                            integration_steps=INTEGRATION_STEPS,
                            stype="OA",
                        )

                    self.sample_calculation_oa = Fitting(
                        self.sample_curves_oa,
                        0.5,
                        self.sampleCA_beamwaist,
                        self.sampleOA_zeroLevel,
                        self.sampleOA_centerPoint,
                        len(self.sample_data_set[0]),
                        line_data,
                    )

                    self.result = self.sample_calculation_oa.manual(
                        zero_level=self.sampleOA_zeroLevel,
                        centerpoint=self.sampleOA_centerPoint,
                        amplitude=0.5,
                        beamwaist=self.sampleCA_beamwaist,
                        z_range=self.z_range,
                        window=self,
                        stype="OA",
                    )
                    self.draw_fitting_line(ftype, stype)
                    self.sampleOA_fittingDone = False

        self.calculate_derived_parameters(ftype)

        caller = "manual"
        self.set_fit_summary(ftype, stype, caller)

    def read_header_params(self, caller: str, ftype: str):
        if caller == "Current Measurement":
            # General parameters
            if not self.customWavelength_checkBox.isChecked():
                self.wavelength_dataFittingTab_doubleSpinBox.setValue(
                    self.wavelength_dataSavingTab_doubleSpinBox.value()
                )
            if not self.customZscanRange_checkBox.isChecked():
                self.zscanRange_doubleSpinBox.setValue(
                    np.abs(
                        self.endPos_doubleSpinBox.value()
                        - self.startPos_doubleSpinBox.value()
                    )
                )
            if not self.customSilicaThickness_checkBox.isChecked():
                self.silicaThickness_dataFittingTab_doubleSpinBox.setValue(
                    self.silicaThickness_dataSavingTab_doubleSpinBox.value()
                )

            if ftype == "Sample":
                self.concentration_dataFittingTab_doubleSpinBox.setValue(
                    self.concentration_dataSavingTab_doubleSpinBox.value()
                )

        elif caller == "Load From File":
            if len(self.header) != 0:
                for hl in self.header:
                    # -------------------------
                    # Wavelength
                    # -------------------------
                    if "Wavelength" in hl:
                        wavelength_match = re.search(
                            r"(([1-9][0-9]*\.?[0-9]*)|(\.[0-9]+))([Ee][+-]?[0-9]+)?",
                            hl,
                        )
                        if wavelength_match is None:
                            continue  # or raise, depending on your logic

                        wavelength = float(wavelength_match.group(1))
                        self.wavelength_dataFittingTab_doubleSpinBox.setValue(
                            wavelength
                        )

                    # -------------------------
                    # Starting / ending position
                    # -------------------------
                    if (
                        not self.customZscanRange_checkBox.isChecked()
                        and "Starting pos" in hl
                    ):
                        starting_pos_match = re.search(
                            r"(([0-9][0-9]*\.?[0-9]*)|(\.[0-9]+))([Ee][+-]?[0-9]+)?",
                            hl,
                        )
                        if starting_pos_match is None:
                            continue

                        starting_pos = float(starting_pos_match.group(1))

                        next_line = self.header[self.header.index(hl) + 1]
                        end_pos_match = re.search(
                            r"(([0-9][0-9]*\.?[0-9]*)|(\.[0-9]+))([Ee][+-]?[0-9]+)?",
                            next_line,
                        )
                        if end_pos_match is None:
                            continue

                        end_pos = float(end_pos_match.group(1))

                        self.zscanRange_doubleSpinBox.setValue(
                            np.abs(end_pos - starting_pos)
                        )

                    # -------------------------
                    # Concentration
                    # -------------------------
                    if ftype == "Sample" and "Concentration" in hl:
                        concentration_match = re.search(
                            r"([0-9]+\.?[0-9]*\s*%)",
                            hl,
                        )
                        if concentration_match is None:
                            continue

                        # Extract "12.5%" → remove spaces → strip "%"
                        raw = concentration_match.group(1).replace(" ", "")
                        self.concentr_percent = float(raw[:-1])

                        self.concentration_dataFittingTab_doubleSpinBox.setValue(
                            self.concentr_percent
                        )

    def solvent_autocomplete(self):
        self.solventDensity_lineEdit.setText(
            str(
                self.solvents[self.solventName_comboBox.currentText()][
                    "density"
                ]
            )
            + " g/cm3"
        )
        self.solventRefrIdx_lineEdit.setText(
            str(self.solvents[self.solventName_comboBox.currentText()]["index"])
        )

    # THREAD CONTROLS
    def print_output(self, returned_value):
        print(returned_value)

    def thread_complete(self):
        print("THREAD COMPLETE!")

    def thread_it(self, func_to_execute):
        # Pass the function to execute
        worker = Worker(
            func_to_execute
        )  # Any other args, kwargs are passed to the run function
        worker.signals.result.connect(self.print_output)
        worker.signals.finished.connect(self.thread_complete)

        if func_to_execute == self.mpositioner.movetostart:
            worker.signals.progress.connect(self.create_raw_log_line)

        # Execute
        self.threadpool.start(worker)

        return worker

    # DIALOG BOXES
    def showdialog(self, msg_type: str, message: str):
        """
        Message type (msg_type) is one of these: 'Error', 'Warning', 'Info'
        """
        args = (msg_type, message)
        match msg_type:
            case "Error":
                button = QMessageBox.critical(self, *args)
            case "Warning":
                button = QMessageBox.warning(self, *args)
            case "Info":
                button = QMessageBox.information(self, *args)

        # if button == QMessageBox.Ok:
        #    pass

    # COLOR THEMES
    @QtCore.pyqtSlot()
    def changeSkinDark(self):
        # PALETTE
        dark_palette = QPalette()

        # Active
        dark_palette.setColor(QPalette.Window, QColor(35, 35, 40))
        dark_palette.setColor(QPalette.WindowText, QColor(200, 200, 200))
        dark_palette.setColor(QPalette.Base, QColor(60, 60, 65))
        dark_palette.setColor(QPalette.AlternateBase, QColor(35, 35, 40))
        # dark_palette.setColor(QPalette.ToolTipBase, QColor(255,255,255))
        # dark_palette.setColor(QPalette.ToolTipText, QColor(255,255,255))
        dark_palette.setColor(QPalette.Text, QColor(200, 200, 200))
        dark_palette.setColor(QPalette.Button, QColor(35, 35, 40))
        dark_palette.setColor(QPalette.ButtonText, QColor(200, 200, 200))
        # dark_palette.setColor(QPalette.BrightText, QColor(255,0,0))
        # dark_palette.setColor(QPalette.Link, QColor(42, 130, 218))
        dark_palette.setColor(QPalette.Highlight, QColor(42, 130, 218))

        # Disabled
        dark_palette.setColor(
            QPalette.Disabled, QPalette.ButtonText, QColor(100, 100, 100)
        )

        QtGui.QGuiApplication.setPalette(dark_palette)

        # Measurement Tab Charts
        self.rms_text.set_color("white")
        self.rms_text.set_bbox(
            dict(
                boxstyle="round",
                facecolor=(60 / 255, 60 / 255, 65 / 255),
                edgecolor=(200 / 255, 200 / 255, 200 / 255),
                alpha=1,
            )
        )

        for chart in self.charts.values():
            chart.fig.patch.set_facecolor((35 / 255, 35 / 255, 40 / 255, 1))
            chart.axes.set_facecolor((35 / 255, 35 / 255, 40 / 255, 1))

            for spine in chart.axes.spines.values():
                spine.set_color((200 / 255, 200 / 255, 200 / 255, 1))

            chart.axes.set_xlabel(
                chart.axes.get_xlabel(),
                fontdict={"color": (200 / 255, 200 / 255, 200 / 255, 1)},
            )
            chart.axes.set_ylabel(
                chart.axes.get_ylabel(),
                fontdict={"color": (200 / 255, 200 / 255, 200 / 255, 1)},
            )
            chart.axes.tick_params(
                axis="both",
                which="both",
                colors=(200 / 255, 200 / 255, 200 / 255, 1),
            )
            chart.axes.grid(
                which="both", color=(60 / 255, 60 / 255, 65 / 255, 1)
            )
            chart.draw_idle()

        # Fitting Tab Charts
        for chart_types in self.fitting_charts.values():
            for chart in chart_types.values():
                chart.fig.patch.set_facecolor((35 / 255, 35 / 255, 40 / 255, 1))
                chart.axes.set_facecolor((35 / 255, 35 / 255, 40 / 255, 1))

                for spine in chart.axes.spines.values():
                    spine.set_color((200 / 255, 200 / 255, 200 / 255, 1))

                chart.axes.set_title(
                    chart.axes.get_title(),
                    fontdict={"color": (200 / 255, 200 / 255, 200 / 255, 1)},
                )
                chart.axes.set_xlabel(
                    chart.axes.get_xlabel(),
                    fontdict={"color": (200 / 255, 200 / 255, 200 / 255, 1)},
                )
                chart.axes.set_ylabel(
                    chart.axes.get_ylabel(),
                    fontdict={"color": (200 / 255, 200 / 255, 200 / 255, 1)},
                )
                chart.axes.tick_params(
                    axis="both",
                    which="both",
                    colors=(200 / 255, 200 / 255, 200 / 255, 1),
                )
                chart.axes.grid(
                    which="both", color=(60 / 255, 60 / 255, 65 / 255, 1)
                )
                chart.draw_idle()

        # STYLESHEET
        stylesheet = """ 
            QTabBar::tab {height: 1em; margin: 0px; padding: 4px; padding-left: 1em; padding-right: 1em; /* height: 1em is expected to recover original setting */
                border-top: 1px solid rgba(60, 60, 65, 1); border-top-left-radius: 3px; border-top-right-radius: 3px;
                border-left: 1px solid rgba(60, 60, 65, 1);
                border-right: 1px solid rgba(60, 60, 65, 1);}
            QTabBar::tab:!selected {margin-top: 3px;}
            QTabBar::tab:selected {height: 15 px; border: 1px solid rgba(42, 130, 218, 1); border-top-left-radius: 3px; border-top-right-radius: 3px;
                background: rgba(42, 130, 218, 1); color: white} /* height: 15 px is expected to cover for slight change */
                                                                 /* in the height when only one tab is present in the tab bar. */
                                                                 /* It is prone to font size changes though. */
            """

        self.setStyleSheet(stylesheet)

    @QtCore.pyqtSlot()
    def changeSkinLight(self):
        # PALETTE
        QtGui.QGuiApplication.setPalette(self.default_palette)

        # Measurement Tab Charts
        self.rms_text.set_color("black")
        self.rms_text.set_bbox(
            dict(
                boxstyle="round", facecolor="white", edgecolor="black", alpha=1
            )
        )

        for chart in self.charts.values():
            chart.fig.patch.set_facecolor((255 / 255, 255 / 255, 255 / 255, 1))
            chart.axes.set_facecolor((255 / 255, 255 / 255, 255 / 255, 1))

            for spine in chart.axes.spines.values():
                spine.set_color((0, 0, 0, 1))

            chart.axes.set_xlabel(
                chart.axes.get_xlabel(), fontdict={"color": (0, 0, 0, 1)}
            )
            chart.axes.set_ylabel(
                chart.axes.get_ylabel(), fontdict={"color": (0, 0, 0, 1)}
            )
            chart.axes.tick_params(
                axis="both", which="both", colors=(0, 0, 0, 1)
            )
            chart.axes.grid(which="both", color="darkgrey")
            chart.draw_idle()

        # Fitting Tab Charts
        for chart_types in self.fitting_charts.values():
            for chart in chart_types.values():
                chart.fig.patch.set_facecolor(
                    (255 / 255, 255 / 255, 255 / 255, 1)
                )
                chart.axes.set_facecolor((255 / 255, 255 / 255, 255 / 255, 1))

                for spine in chart.axes.spines.values():
                    spine.set_color((0, 0, 0, 1))

                chart.axes.set_title(
                    chart.axes.get_title(), fontdict={"color": (0, 0, 0, 1)}
                )
                chart.axes.set_xlabel(
                    chart.axes.get_xlabel(), fontdict={"color": (0, 0, 0, 1)}
                )
                chart.axes.set_ylabel(
                    chart.axes.get_ylabel(), fontdict={"color": (0, 0, 0, 1)}
                )
                chart.axes.tick_params(
                    axis="both", which="both", colors=(0, 0, 0, 1)
                )
                chart.axes.grid(which="both", color="darkgrey")
                chart.draw_idle()

        # STYLESHEET
        stylesheet = """ 
            QTabBar::tab {height: 1em; margin: 0px; padding: 4px; padding-left: 1em; padding-right: 1em; /* height: 1em is expected to recover original setting */
                border-top: 1px solid rgba(200, 200, 205, 1); border-top-left-radius: 3px; border-top-right-radius: 3px;
                border-left: 1px solid rgba(200, 200, 205, 1);
                border-right: 1px solid rgba(200, 200, 205, 1);
                }
            QTabBar::tab:!selected {margin-top: 3px;}
            QTabBar::tab:selected {height: 15px; border: 1px solid rgba(42, 130, 218, 1); border-top-left-radius: 3px; border-top-right-radius: 3px;
                background: rgba(42, 130, 218, 1); color: white} /* height: 15 px is expected to cover for slight change */
                                                                 /* in the height when only one tab is present in the tab bar. */
                                                                 /* It is prone to font size changes though. */
            """

        self.setStyleSheet(stylesheet)


class MotorPositioner(QObject):
    def movetostart(self, progress_callback):
        start_pos = window.startPos_doubleSpinBox.value()
        end_pos = window.endPos_doubleSpinBox.value()
        if window.motor.position <= (start_pos + end_pos) / 2:
            window.motor.move_to(start_pos)
            window.where_to_start = "start"
        else:
            window.motor.move_to(end_pos)
            window.where_to_start = "end"

        if window.motor.is_in_motion:
            print("Moving to starting position")

        while window.motor.is_in_motion:
            continue

        self.run(progress_callback)

    def movetocustompos(self, *args, **kwargs):
        text_val = window.custom_pos_dialog.new_pos.text()
        target = float(text_val.replace(",", "."))

        if target > 100:
            return

        window.motor.move_to(target)
        while window.motor.is_in_motion:
            continue

        window.current_pos_chooser.setEnabled(True)

    def moveby(self, *args, **kwargs):
        move_step = (
            window.endPos_doubleSpinBox.value()
            - window.startPos_doubleSpinBox.value()
        ) / window.stepsScan_spinBox.value()
        if window.where_to_start == "end":
            window.motor.move_by(-move_step, blocking=True)
        else:
            window.motor.move_by(move_step, blocking=True)

        # while window.motor.is_in_motion == True:
        #    continue

    def movehome(self, *args, **kwargs):
        if not window.motor.has_homing_been_completed:
            window.motor.move_home()

            time.sleep(0.2)  # otherwise it may not move_home()
            if window.motor.is_in_motion:
                print("Homing now")
            while not window.motor.has_homing_been_completed:
                continue  # wait until homing is completed
            time.sleep(
                0.2
            )  # wait a little more (so the motor.position gets exactly "0" position)

        return "Homing performed!"

    def run(self, progress_callback):
        window.data_acquisition_complete = False
        window.data_reversed = False  # when backwards scan is performed, it later gets reversed (the data_reverse() method)

        time.sleep(
            0.2
        )  # Sometimes the first datapoint is collected before the motor has settled

        nos = window.stepsScan_spinBox.value()

        for step in range(nos + 1):
            if window.experiment_stopped:
                window.experiment_stopped = False
                window.running = False
                break
            self.step = step
            ### Create a task
            with nidaqmx.Task() as task:
                # Create MultiChannel "channel"
                task.ai_channels.add_ai_voltage_chan(
                    window.detector_core_name
                    + f"0:{window.number_of_channels_used}"
                )  # "Dev1/ai0:3"

                # Start Digital Edge
                task.triggers.start_trigger.dig_edge_src = "/Dev1/PFI0"
                task_trigger_src = task.triggers.start_trigger.dig_edge_src

                # Sample Clock
                rate0 = 1000  # 1 kHz (repetition rate of the laser)
                task.timing.cfg_samp_clk_timing(
                    rate0,
                    source=task_trigger_src,
                    active_edge=Edge.FALLING,
                    sample_mode=AcquisitionType.FINITE,
                    samps_per_chan=window.samplesStep_spinBox.value(),
                )

                ### Data acquisition
                reader = nidaqmx.stream_readers.AnalogMultiChannelReader(
                    task.in_stream
                )
                values_read = np.zeros(
                    (
                        window.number_of_channels_used + 1,
                        window.samplesStep_spinBox.value(),
                    ),
                    dtype=np.float64,
                )
                # Start Task
                task.start()

                time.sleep(
                    0.2
                )  # THIS IS THE TIME NEEDED TO COLLECT ALL SAMPLES

                # Acquire data
                window.data["positions"].append(window.motor.position)
                reader.read_many_sample(
                    values_read,
                    number_of_samples_per_channel=window.samplesStep_spinBox.value(),
                )

                # Take mean for each channel and append to "data" dictionary
                data_mean = np.mean(values_read, axis=1)

                for chan_no in range(window.number_of_channels_used):
                    window.data["absolute"][chan_no].append(data_mean[chan_no])
                # These loops need to be separated because "absolute" list at current step MUST BE defined at current step when accessed by "relative" list
                for chan_no in range(window.number_of_channels_used):
                    # Take last/current value (at index [step]) in "absolute" values list and divide by channel [1] (reference)
                    # And append to "relative" values list
                    window.data["relative"][chan_no].append(
                        data_mean[chan_no] / window.data["absolute"][1][step]
                    )

                # 2) display data

                for type in window.charts.keys():
                    for chan_no in range(window.number_of_channels_used):
                        window.measurement_lines[type][chan_no].set_xdata(
                            window.data["positions"]
                        )  # and set data of lines in "lines" dictionary to empty lists of positions and values
                        window.measurement_lines[type][chan_no].set_ydata(
                            window.data[type][chan_no]
                        )

                y = window.data["absolute"][1]
                window.rms_value = (
                    np.abs(np.sqrt(np.mean([yi**2 for yi in y])) - y[0]) / y[0]
                )
                window.rms_text.set_text(
                    f"RMS noise = {window.rms_value * 100:.3f}%"
                )

                window.measurement_plot_rescale(
                    window.focusAt_comboBox.currentText()
                )

                progress_callback.emit(step)

                # 3) move the motor
                if step == nos:  # prevent useless additional step
                    window.data_acquisition_complete = True
                    break
                else:
                    window.mpositioner.moveby()

                task.stop()

        # EMIT SOUND AT FINISH
        duration = 500  # milliseconds
        freq = 800  # Hz
        for _ in range(3):
            winsound.Beep(freq, duration)
            time.sleep(0.05)

        window.running = False


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)

    window = Window()
    app.setStyle("Fusion")
    window.default_palette = QtGui.QGuiApplication.palette()
    window.changeSkinDark()  # Make sure the additional changes are applied
    app.exec_()
