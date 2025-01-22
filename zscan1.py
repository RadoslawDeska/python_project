"""Z-scan measurement program allows for automated data collection
from simultaneous measurement of open- and closed-aperture Z-scan traces.

Additionally, the program allows for data analysis by fitting the Z-scan traces
to theoretical functions developed by Mansoor Sheik-Bahae and Eric van Stryland
published in \"Characterization Techniques and Tabulations for Organic Nonlinear Materials\"
M. G. Kuzyk and C. W. Dirk, Eds., page 655-692, Marcel Dekker, Inc., 1998
"""

__author__ = "Radosław Deska"
__version__ = "0.1.8"

# TRY INTRODUCING THE BATCH SCAN!!
# CHANGES AFTER 0.1.7
# 1) SCRIPT LOGICS:
#      - Motor now moves always in the direction expected based on initial position (even after Red Button)
#      - Red Button causes the MotorPositioner.run() to raise ExperimentStopped exception to break off the experiment completely
#      - Closing the program while the thread executing the MotorPositioner.run() method is now causing the thread to finish immediately
#      - Added settings entry to control the Beep at the end of the experiment
# 2) GUI:
#      - Added system of multiple scan runs.
#          -> Still requires saving to H5PY and txt temporary files, displaying in Saving Tab requires
#                   repair (prepend, append - this is stupidly done if number of runs is > 1)
# 3) SOLVENTS:
#      -
# 4) SAMPLES:
#      -

# TO DO LIST:
# 1) SCRIPT LOGICS:
#      - Threads should not modify main thread data (USE SIGNALS)
# 2) Use thorlabs_apt_device to shorten the initialization process (in other words: by-pass APT.dll). This won't work for simulation mode!

import argparse
import functools
import logging
import os
import re
import sys
import threading
import time
import traceback
from typing import Tuple

import matplotlib
import numpy as np
from numpy.typing import NDArray
from PyQt5 import QtGui, QtWidgets
from PyQt5.QtCore import QEvent, QObject, QSettings, QThreadPool, pyqtSignal
from PyQt5.QtWidgets import QApplication, QFileDialog, QSlider
from scipy.signal import medfilt
from sigfig import round as error_rounding

from classes.mainwindow.managers import FocusManager
from lib.cursors import BlittedCursor
from lib.gui_modifiers import GUI_Modifiers
from lib.settings import *
from lib.theory import Fitting, Integration

matplotlib.rcParams.update({"font.size": 7})

# CONSTANTS
SILICA_BETA = 0
N_COMPONENTS = 8  # number of electric field components (for Gaussian decomposition)
INTEGRATION_STEPS = 30  # accuracy of integration infinitesimal element, dx.
CUVETTE_PATH_LENGTH = 0.001  # [m] path length inside cuvette
SOLVENT_T_SLIDER_MAX = 1
MAX_DPHI0 = 3.142  # maximum DeltaPhi0 for silica (for sliders)
SIG_FIG_MAN_FIT = 3  # number of significant digits in rounding values while manually fitting the curves

# @Adam notes:
# 1. Separation of UI frrom computations / domain logic
#     Data representation (GUI / FE)   ---- API ---- Server logic / BE
#         REST API (webdev)
#           HTTP (Method, host, headers, body)
#               /measurements/{type}/{id}
#           m = Measurements("somejsonfile.json")
#            m.get_results(xyz)
# 2. Common utils
# 3. linters: black, mypy (typing), ruff
# 4. Flow -> think in terms of data flow (functional programming)


def time_it(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        print(f"Function {func.__name__} took {end_time - start_time:.6f} seconds to run")
        return result

    return wrapper

class Filter(QObject):
    clickedOutside = pyqtSignal(object)
    
    def __init__(self, target_widgets, parent=None):
        super().__init__(parent)
        self.parent_ = parent
        self.target_widgets = target_widgets

    def eventFilter(self, obj, event):
        if event.type() == QEvent.MouseButtonPress:
            for widget in self.target_widgets:
                if widget.rect().contains(widget.mapFromGlobal(event.globalPos())):
                    return super().eventFilter(obj, event)
            print(f"{FocusManager()._previously_focused=}")
            self.clickedOutside.emit(FocusManager()._previously_focused)
        
        return super().eventFilter(obj, event)


class Window(QtWidgets.QMainWindow, GUI_Modifiers):
    from classes.mainwindow.colorthemes import changeReadOnlyStyle, changeSkinDark, changeSkinLight
    from classes.mainwindow.dialogs import showdialog
    from classes.mainwindow.initialization_procedures import (
        initialize,
        modify_gui_to_initial_looks,
        set_additional_variables,
        set_initial_states,
        setup_triggers,
    )
    from classes.mainwindow.measurement_plotting import (
        clear_measurement,
        plot_updated_data,
        rescale_measurement_plots,
    )
    from classes.mainwindow.measurement_processing import (
        create_string_list_for_temp_file,
        create_string_list_to_display,
        data_save,
        process_current_datapoint,
        switch_widget_items,
        update_fullLog,
        write_raw_data_log,
        write_temp_file,
    )
    from classes.mainwindow.motor_control import set_motor_to_start, stop_experiment
    from classes.mainwindow.threadcontrols import (
        acquisition_complete,
        acquisition_starting,
        experiment_thread_finished,
        homing_complete,
        homing_start,
        print_output,
        thread_complete,
        thread_it,
    )
    from classes.mainwindow.user_triggered_signals import fitting_tab_signals, on_focus_out

    resume_scanning = pyqtSignal()

    # INITIALIZATION
    def __init__(self, settings):
        super(Window, self).__init__()
        self.settings = QSettings(settings, QSettings.IniFormat)

        self.setupUi(self)
        apply_settings(self, self.settings)

        self.experiment_stopped = threading.Event()  # a flag to control experiment run time
        self.threadpool = QThreadPool()
        # These have to be called prior to all else to provide names in namespace
        self.set_initial_states()
        self.set_additional_variables()
        # These have to be called after initial states and variables are set
        self.modify_gui_to_initial_looks()
        self.setup_triggers()
        
        '''Fix the issue of persistent focus by detecting clicks outside the specified widgets,
        even in areas that do not accept focus, and handling the event appropriately'''
        self.widgets_to_check = [self.codeOfSample_comboBox,
                                 self.codeOfSample_comboBox.view(),
                                self.silicaThickness_dataSavingTab_doubleSpinBox,
                                self.concentration_dataSavingTab_doubleSpinBox,
                                self.solventName_dataSavingTab_comboBox,
                                ]
        for widget in self.widgets_to_check:
            widget.elementsReceived.emit(self.widgets_to_check)
            
        self.filter = Filter(self.widgets_to_check, self)
        self.filter.clickedOutside.connect(self.on_focus_out)
        self.installEventFilter(self.filter)

    # SELECT DIRECTORY
    def choose_dir(self, caller=""):
        if caller == "DataSaving":
            path = self.mainDirectory_lineEdit.text()
        elif caller == "DataFitting":
            path = self.dataDirectory_lineEdit.text()

        p = QFileDialog.getExistingDirectory(self, "Select a directory", path)
        if p != "":  # This keeps old path in directory QLineEdit lines, if dialog is closed with Cancel
            if caller == "DataSaving":
                self.mainDirectory_lineEdit.setText(p.replace("/", "\\"))
            elif caller == "DataFitting":
                self.dataDirectory_lineEdit.setText(p.replace("/", "\\"))

    # DATA FITTING
    def data_loader(self, caller: str, ftype: str):
        """Load data either from file or from current measurement ('caller' argument).\n
        'ftype' is passed by proper button for loading from file, or by proper option in "Cuvette type" combobox in "Data Saving" tab"""

        def prepare_to_fit_manually(caller: str, ftype: str):
            # Read parameters
            self.read_header_params(caller=caller, ftype=ftype)
            self.data_display(self.data_set, ftype)

            self.switch_fitting_to_on_state(ftype)
            # Adjust centerPoint slider (prevent unexpected shifting of the slider and the curve to a side of the plot)
            match ftype:
                case "silica":
                    self.silica_centerPoint_slider.setMaximum(self.silica_nop - 1)
                    # Reset slider to center position
                    self.silica_centerPoint_slider.setValue(int(self.silica_nop / 2))

            if ftype == "silica":
                self.silica_autofit_done = False  # this is to start afresh with fitting

            self.fit_manually(ftype, stype="CA")
            self.fit_manually(ftype, stype="OA")

        if caller == "Current Measurement":
            # Fill in names in QLineEdits
            self.dataDirectory_lineEdit.setText(self.accurate_path + "\\")
            self.mainTabs.setCurrentIndex(2)

            match ftype:
                case "silica":
                    self.silicaFilename_lineEdit.setText(self.rawLogFilename_lineEdit.text())
                    self.fittingTabs.setCurrentIndex(0)
                    self.silicaAperture_tabWidget.setCurrentIndex(0)
                case "solvent":
                    self.solventFilename_lineEdit.setText(self.rawLogFilename_lineEdit.text())
                    self.fittingTabs.setCurrentIndex(1)
                    self.solventAperture_tabWidget.setCurrentIndex(0)
                case "sample":
                    self.sampleFilename_lineEdit.setText(self.rawLogFilename_lineEdit.text())
                    self.fittingTabs.setCurrentIndex(2)
                    self.sampleAperture_tabWidget.setCurrentIndex(0)

            prepare_to_fit_manually(caller=caller, ftype=ftype)

        elif caller == "Load From File":
            # Load data
            def load_data(caller, ftype):
                """Fills proper frame in GUI with file information, sets current tab to lead the user,\n
                loads data on screen, toggles to active state the fitting controls and calls for the initial fit\n
                for both CA and OA traces.
                """
                # if the header is correct, uncheck "set custom" checkboxes
                if self.header_correct:
                    for target in targets_checkboxes:
                        target.setChecked(False)
                # otherwise uncheck only those checkboxes that were found in the header
                else:
                    for i, matched in enumerate(header_matches[:-1]):  # apart from the last element which is the header beacon
                        targets_checkboxes[i].setChecked(not matched)

                # Fill in names in QLineEdits
                match ftype:
                    case "silica":
                        self.silicaFilename_lineEdit.setText(fname)
                        self.fittingTabs.setCurrentIndex(0)
                        self.silicaAperture_tabWidget.setCurrentIndex(0)
                    case "solvent":
                        self.solventFilename_lineEdit.setText(fname)
                        self.fittingTabs.setCurrentIndex(1)
                        self.solventAperture_tabWidget.setCurrentIndex(0)
                    case "sample":
                        self.sampleFilename_lineEdit.setText(fname)
                        self.fittingTabs.setCurrentIndex(2)
                        self.sampleAperture_tabWidget.setCurrentIndex(0)

                # Get data
                data = np.genfromtxt(p, skip_header=last_header_line)
                self.data_set = data[:, 0], data[:, 1], data[:, 2], data[:, 3]  # , data[:,4] not using the last column with zeros

                prepare_to_fit_manually(caller=caller, ftype=ftype)

            p = QFileDialog.getOpenFileName(self, "Select full description file", os.path.normpath(self.dataDirectory_lineEdit.text()))
            if p[0] != "":  # This keeps old filename in given file type QLineEdit lines, if dialog is closed with Cancel
                p = p[0]
                fname = p.split("/")[-1]
                self.dataDirectory_lineEdit.setText("\\".join(p.split("/")[:-1]))

                # Header check and manipulation
                with open(p, "r") as file:
                    # Correct header must include these and a beacon at the end in the form of "SNo" substring
                    required_matches = ["silica thickness", "Concentration", "Wavelength", "Starting pos", "Ending pos", "SNo"]
                    targets_checkboxes = [
                        self.customSilicaThickness_checkBox,
                        self.customConcentration_checkBox,
                        self.customWavelength_checkBox,
                        self.customZscanRange_checkBox,
                        self.customZscanRange_checkBox,
                    ]
                    missing_matches = []
                    header_matches = len(required_matches) * [False]
                    self.header = []
                    last_header_line = 0

                    for line_no, line in enumerate(file):
                        # first look for the header end beacon
                        if re.match(rf"\b{required_matches[-1]}\b", line):
                            last_header_line = line_no + 3
                            break

                    if last_header_line > 0:  # Header found
                        # look for other matches
                        for i, match in enumerate(required_matches[:-1]):
                            file.seek(0)
                            for line_no, line in enumerate(file):
                                if line_no == last_header_line:  # do not look further than header length
                                    break

                                if re.match(rf"\b{match}\b", line):  # lookup whole words
                                    header_matches[i] = True
                                    self.header.append(line.strip())

                            if not header_matches[i]:
                                missing_matches.append((i, match))

                        if missing_matches:
                            self.header_correct = False

                            # Make nicely looking response in the dialog box
                            if len(missing_matches) > 1:
                                missing = ", ".join(miss for (_, miss) in missing_matches[:-1])
                                missing = " and ".join((missing, missing_matches[-1][1]))
                            else:
                                missing = missing_matches[-1][1]

                            self.showdialog(
                                "Warning", ("Missing parameters!\n\n" f"{missing} not found in the header.\n\nSet measurement parameters manually.")
                            )

                        else:
                            self.header_correct = True

                        load_data(caller=caller, ftype=ftype)

                    else:  # Header not found
                        self.header_correct = False
                        self.showdialog("Warning", ("Header not found!\n\n" "Set measurement parameters manually or load a file with the header."))

                        load_data(caller=caller, ftype=ftype)

    def enable_custom(self, o: str):
        """Toggles readOnly parameter on `o` element from GUI. Uses match-case structure with `o` parameter to match\n
        to speed up processing at each call.

        Args:
            o (str): case for parameter to toggle
        """
        match o:
            case "ApertureDiameter":
                if self.customApertureDiameter_checkBox.isChecked() is True:
                    self.apertureDiameter_doubleSpinBox.setReadOnly(False)
                else:
                    self.apertureDiameter_doubleSpinBox.setReadOnly(True)
            case "ApertureDistance":
                if self.customApertureToFocusDistance_checkBox.isChecked() is True:
                    self.apertureToFocusDistance_doubleSpinBox.setReadOnly(False)
                else:
                    self.apertureToFocusDistance_doubleSpinBox.setReadOnly(True)
            case "SilicaThickness":
                if self.customSilicaThickness_checkBox.isChecked() is True:
                    self.silicaThickness_dataFittingTab_doubleSpinBox.setReadOnly(False)
                else:
                    self.silicaThickness_dataFittingTab_doubleSpinBox.setReadOnly(True)
            case "Wavelength":
                if self.customWavelength_checkBox.isChecked() is True:
                    self.wavelength_dataFittingTab_doubleSpinBox.setReadOnly(False)
                else:
                    self.wavelength_dataFittingTab_doubleSpinBox.setReadOnly(True)
            case "ZscanRange":
                if self.customZscanRange_checkBox.isChecked() is True:
                    self.zscanRange_doubleSpinBox.setReadOnly(False)
                else:
                    self.zscanRange_doubleSpinBox.setReadOnly(True)
            case "Concentration":
                if self.customConcentration_checkBox.isChecked() is True:
                    self.concentration_dataFittingTab_doubleSpinBox.setReadOnly(False)
                else:
                    self.concentration_dataFittingTab_doubleSpinBox.setReadOnly(True)
            case "SolventBeamwaist":
                if self.solventCA_customBeamwaist_checkBox.isChecked() is False:
                    self.solventCA_RayleighLength_slider.setEnabled(False)
                    if hasattr(window, "silica_beamwaist"):
                        self.solventCA_RayleighLength_slider.valueChanged.disconnect()
                        self.solventCA_RayleighLength_slider.setValue(int(round(self.silica_beamwaist * 1e6)))
                        self.solventCA_RayleighLength_slider.setValue(int(round(self.silica_beamwaist * 1e6)))
                        self.solventCA_RayleighLength_slider.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="CA"))
                        self.solventCA_beamwaistSummary_doubleSpinBox.setValue(self.silica_beamwaist * 1e6)
                        self.solventCA_beamwaistSummary_doubleSpinBox.setValue(self.silica_beamwaist * 1e6)

                else:
                    self.solventCA_RayleighLength_slider.setEnabled(True)
            case "SolventCenterPoint":  # OA case
                if self.solventOA_customCenterPoint_checkBox.isChecked() is False:
                    self.solventOA_centerPoint_slider.setEnabled(False)
                    # self.solventOA_centerPoint_doubleSpinBox.setValue(self.solventCA_centerPoint_doubleSpinBox.value())
                else:
                    self.solventOA_centerPoint_slider.setEnabled(True)
        
        self.changeReadOnlyStyle()

    def slider_fit_manually_connect(self, current_slider: QSlider, mode):
        """This function is intended to connect and disconnect other sliders on demand, to prevent them from updating the fitting line:
        1) "Disconnect" means current slider should be kept and all others should disconnect from their method
        2) "Connect" means reconnect all sliders to their method"""
        silica_sliders = [self.silica_RayleighLength_slider, self.silica_centerPoint_slider, self.silica_zeroLevel_slider, self.silica_DPhi0_slider]
        # solventCA_sliders = [self.solventCA_RayleighLength_slider, self.solventCA_centerPoint_slider, self.solventCA_zeroLevel_slider, self.solventCA_DPhi0_slider]
        # solventOA_sliders = [self.solventOA_deltaZpv_slider, self.solventOA_centerPoint_slider, self.solventOA_zeroLevel_slider, self.solventOA_deltaTpv_slider]
        # sampleCA_sliders = [self.sampleCA_deltaZpv_slider, self.sampleCA_centerPoint_slider, self.sampleCA_zeroLevel_slider, self.sampleCA_deltaTpv_slider]
        # sampleOA_sliders = [self.sampleOA_deltaZpv_slider, self.sampleOA_centerPoint_slider, self.sampleOA_zeroLevel_slider, self.sampleOA_deltaTpv_slider]

        available_sliders = [silica_sliders]  # , solventCA_sliders, solventOA_sliders, sampleCA_sliders, sampleOA_sliders]
        self.disconnected_sliders = "All"

        match mode:
            case "Disconnect":
                for sliders in available_sliders:
                    try:
                        sliders.remove(current_slider)
                        for s in sliders:
                            s.disconnect()
                            self.disconnected_sliders = "silica"

                    except ValueError:
                        pass  # If not found

            case "Connect":
                match self.disconnected_sliders:
                    case "silica":
                        for s in available_sliders[0]:
                            s.valueChanged.connect(lambda: self.fit_manually(ftype="silica", stype="CA", activated_by=self.disconnected_sliders))
                    case "silicaOA":
                        for s in available_sliders[1]:
                            s.valueChanged.connect(lambda: self.fit_manually(ftype="silica", stype="OA", activated_by=self.disconnected_sliders))
                    case "solventCA":
                        for s in available_sliders[2]:
                            s.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="CA", activated_by=self.disconnected_sliders))
                    case "solventOA":
                        for s in available_sliders[3]:
                            s.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="OA", activated_by=self.disconnected_sliders))
                    case "sampleCA":
                        for s in available_sliders[4]:
                            s.valueChanged.connect(lambda: self.fit_manually(ftype="sample", stype="CA", activated_by=self.disconnected_sliders))
                    case "sampleOA":
                        for s in available_sliders[5]:
                            s.valueChanged.connect(lambda: self.fit_manually(ftype="sample", stype="OA", activated_by=self.disconnected_sliders))
                    case "All":
                        for s in available_sliders[0]:
                            s.valueChanged.connect(lambda: self.fit_manually(ftype="silica", stype="CA"))
                        # for s in available_sliders[1]:
                        #     s.valueChanged.connect(lambda: self.fit_manually(ftype="silica", stype="OA"))
                        # for s in available_sliders[2]:
                        #     s.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="CA"))
                        # for s in available_sliders[3]:
                        #     s.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="OA"))
                        # for s in available_sliders[4]:
                        #     s.valueChanged.connect(lambda: self.fit_manually(ftype="sample", stype="CA"))
                        # for s in available_sliders[5]:
                        #     s.valueChanged.connect(lambda: self.fit_manually(ftype="sample", stype="OA"))

            case None:
                pass

    def switch_fitting_to_on_state(self, ftype: str):
        match ftype:
            case "silica":
                self.silica_RayleighLength_slider.setEnabled(True)
                self.silica_centerPoint_slider.setEnabled(True)
                self.silica_zeroLevel_slider.setEnabled(True)
                self.silica_DPhi0_slider.setEnabled(True)
                self.silica_fit_pushButton.setEnabled(True)
                self.silica_filterSize_slider.setEnabled(True)

            case "solvent":
                if self.solventCA_customBeamwaist_checkBox.isChecked() is True:
                    self.solventCA_RayleighLength_slider.setEnabled(True)
                else:
                    self.solventCA_RayleighLength_slider.setEnabled(False)
                self.solventCA_centerPoint_slider.setEnabled(True)
                self.solventCA_zeroLevel_slider.setEnabled(True)
                self.solventCA_DPhi0_slider.setEnabled(True)
                self.solventCA_fit_pushButton.setEnabled(True)
                self.solventCA_filterSize_slider.setEnabled(True)
                self.solventCA_customBeamwaist_checkBox.setEnabled(True)

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

            case "sample":
                pass

    def set_fit_summary(self, ftype: str, stype: str, caller=""):
        """Sets proper number of digits in summary display and shows parameters (with errors) after (automatic) fitting.

        Args:
            ftype (str): sample type ("silica","solvent","sample")
            stype (str): Z-scan mode ("CA", "OA")
            caller (str, optional): Which fit type has called this function. Defaults to "".
        """

        # DISPLAY VALUES
        match ftype:
            case "silica":
                if caller == "auto":
                    self.silica_deltaPhi0Summary_label.setText(str(error_rounding(self.silica_DPhi0, self.silica_DPhi0Error)))  # rad
                    self.silica_laserIntensitySummary_label.setText(str(error_rounding(self.laserI0 * 1e-13, self.laserI0Error * 1e-13)))  # [GW/cm2]
                    self.silica_beamwaistSummary_label.setText(
                        str(error_rounding(self.silica_beamwaist * 1e6, self.silica_beamwaistError * 1e6))
                    )  # [um] radius in focal point
                    self.silica_rayleighRangeSummary_label.setText(
                        str(error_rounding(self.silica_rayleighLength * 1e3, self.silica_rayleighLengthError * 1e3))
                    )  # [mm]
                    self.numericalApertureSummary_label.setText(str(error_rounding(self.numericalAperture, self.numericalApertureError)))
                else:
                    self.silica_deltaPhi0Summary_label.setText(
                        str(error_rounding(self.silica_DPhi0, sigfigs=SIG_FIG_MAN_FIT, warn=False)) + " ± #.##"
                    )  # rad
                    self.silica_laserIntensitySummary_label.setText(
                        str(error_rounding(self.laserI0 * 1e-13, sigfigs=SIG_FIG_MAN_FIT, warn=False)) + " ± #.##"
                    )  # [GW/cm2]
                    self.silica_beamwaistSummary_label.setText(
                        str(error_rounding(self.silica_beamwaist * 1e6, sigfigs=SIG_FIG_MAN_FIT, warn=False)) + " ± #.##"
                    )  # [um] radius in focal point
                    self.silica_rayleighRangeSummary_label.setText(
                        str(error_rounding(self.silica_rayleighLength * 1e3, sigfigs=SIG_FIG_MAN_FIT, warn=False)) + " ± #.##"
                    )  # [mm]
                    self.numericalApertureSummary_label.setText(
                        str(error_rounding(self.numericalAperture, sigfigs=SIG_FIG_MAN_FIT, warn=False)) + " ± #.##"
                    )

            case "solvent":
                if stype == "CA":
                    self.solventCA_deltaPhi0Summary_doubleSpinBox.setValue(self.solvent_DPhi0)
                    self.solventCA_n2Summary_doubleSpinBox.setValue(self.solvent_n2 * 1e22)  # display n2 in multiples of 10^-9 cm^2/GW
                    self.solventCA_beamwaistSummary_doubleSpinBox.setValue(self.solvent_beamwaist * 1e6)  # [um] radius in focal point
                    self.solventCA_deltaPhi0Summary_doubleSpinBox.setValue(self.solvent_DPhi0)
                    self.solventCA_n2Summary_doubleSpinBox.setValue(self.solvent_n2 * 1e22)  # display n2 in multiples of 10^-9 cm^2/GW
                    self.solventCA_beamwaistSummary_doubleSpinBox.setValue(self.solvent_beamwaist * 1e6)  # [um] radius in focal point
                    self.solventCA_rayleighRangeSummary_doubleSpinBox.setValue(self.solvent_rayleighLength * 1e3)  # [mm]

                elif stype == "OA":
                    self.solventOA_n2Summary_doubleSpinBox.setValue(self.solvent_n2 * 1e22)  # display n2 in multiples of 10^-9 cm^2/GW
                    self.solventOA_n2Summary_doubleSpinBox.setValue(self.solvent_n2 * 1e22)  # display n2 in multiples of 10^-9 cm^2/GW
                    # self.solventOA_TSummary_doubleSpinBox.setValue(self.solventOA_T)
                    # TEMPORARILY DISABLED
                    # self.solventOA_betaSummary_doubleSpinBox.setValue(self.solvent_beta) # CUVETTE_PATH_LENGTH is the solvent/sample thickness assuming no one-photon absorption

            case "sample":
                pass

    def get_general_parameters(self):
        """Gets the values from 'General Parameters' GUI frame and assigns them to variables with basic SI units (mm -> m):\n
        `self.l_silica`, `self.lda`, `self.z_range`, `self.ra`, `self.d0`\n
        and then calculate parameters that depend only on them and not on fitting process:\n
        `self.silica_n2`
        """
        self.l_silica = self.silicaThickness_dataFittingTab_doubleSpinBox.value() * 1e-3  # [m] silica thickness
        self.lda = self.wavelength_dataFittingTab_doubleSpinBox.value() * 1e-9  # [m] wavelength
        self.z_range = self.zscanRange_doubleSpinBox.value() * 1e-3  # [m] z-scan range
        self.ra = self.apertureDiameter_doubleSpinBox.value() / 2 * 1e-3  # [m] CA aperture radius
        self.d0 = self.apertureToFocusDistance_doubleSpinBox.value() * 1e-3  # [m] distance from focal point to aperture plane
        self.silica_n2 = (
            2.8203e-20 - 3e-27 / (self.lda) + 2e-33 / (self.lda) ** 2
        )  # [m2/W] Bandar A. Babgi formula (based on David Milam tables for n2)

    def get_curve_interpretation(self, ftype, stype, from_what: str, on_data_load=False):
        """Interprets the curve based on its geometry according to some approximated formulas (manual fitting)\n
        or based on automatic fitting parameters (automatic fitting) and calls for error estimation."""
        match ftype:
            case "silica":
                match from_what:
                    case "from_geometry":  # fit_manually
                        # read sliders
                        self.silica_zeroLevel = self.silica_zeroLevel_slider.value() / 100
                        self.silica_centerPoint = np.round(self.silica_centerPoint_slider.value() - self.silica_nop / 2)
                        if on_data_load:
                            self.silica_DPhi0, self.silica_beamwaist, self.silica_rayleighLength = self.params_from_geometry(ftype, "CA")
                        else:
                            self.silica_DPhi0 = self.silica_DPhi0_slider.value() * MAX_DPHI0 / self.silica_DPhi0_slider.maximum()
                            self.silica_rayleighLength = (
                                self.silica_RayleighLength_slider.value() * self.z_range / 2 / self.silica_RayleighLength_slider.maximum()
                            )
                            self.silica_beamwaist = np.sqrt(self.silica_rayleighLength * self.lda / np.pi)
                    case "from_autofit":  # fit_automatically
                        self.silica_zeroLevel, self.silica_zeroLevelError = (
                            self.silica_minimizerResult.params["Zero"].value,
                            self.silica_minimizerResult.params["Zero"].stderr,
                        )
                        self.silica_centerPoint, self.silica_centerPointError = (
                            self.silica_minimizerResult.params["Center"].value,
                            self.silica_minimizerResult.params["Center"].stderr,
                        )
                        self.silica_DPhi0, self.silica_DPhi0Error = (
                            self.silica_minimizerResult.params["DPhi0"].value,
                            self.silica_minimizerResult.params["DPhi0"].stderr,
                        )
                        self.silica_beamwaist, self.silica_beamwaistError = (
                            self.silica_minimizerResult.params["Beamwaist"].value,
                            self.silica_minimizerResult.params["Beamwaist"].stderr,
                        )

                        self.silica_rayleighLength = float(np.pi * self.silica_beamwaist**2 / self.lda)
                        self.silica_rayleighLengthError = (
                            np.pi / 2 / self.lda * self.silica_minimizerResult.params["Beamwaist"].value * self.silica_beamwaistError
                        )  # [m] Rayleigh length

                        self.laserI0Error = self.silica_DPhi0Error * self.lda / (2 * np.pi * self.l_silica * self.silica_n2)  # [W/m2]
                        self.numericalApertureError = (
                            self.silica_beamwaistError / self.silica_rayleighLength
                            + self.silica_beamwaist * self.silica_rayleighLengthError / self.silica_rayleighLength**2
                        )

                self.laserI0 = self.silica_DPhi0 * self.lda / (2 * np.pi * self.l_silica * self.silica_n2)  # [W/m2]
                self.numericalAperture = self.silica_beamwaist / self.silica_rayleighLength
                if np.isnan(self.numericalAperture):
                    self.numericalAperture = 0

            case "solvent":
                match stype:
                    case "CA":
                        match from_what:
                            case "from_geometry":
                                self.solventCA_zeroLevel = self.solventCA_zeroLevel_slider.value() / 100
                                self.solventCA_centerPoint = np.round(self.solventCA_centerPoint_slider.value() - self.solvent_nop / 2)
                                if on_data_load:
                                    self.solvent_DPhi0, self.solvent_beamwaist, self.solvent_rayleighLength = self.params_from_geometry(ftype, "CA")
                                else:
                                    self.solvent_DPhi0 = self.solventCA_DPhi0_slider.value() * MAX_DPHI0 / self.solventCA_DPhi0_slider.maximum()
                                    self.solvent_rayleighLength = (
                                        self.solventCA_RayleighLength_slider.value()
                                        * self.z_range
                                        / 2
                                        / self.solventCA_RayleighLength_slider.maximum()
                                    )
                                    self.solvent_beamwaist = np.sqrt(self.solvent_rayleighLength * self.lda / np.pi)

                                if ftype == "solvent" and self.solventCA_customBeamwaist_checkBox.isChecked() is False:
                                    self.solvent_beamwaist, self.solvent_rayleighLength = self.silica_beamwaist, self.silica_rayleighLength
                            case "from_autofit":
                                pass

                        self.solvent_n2 = (
                            self.solvent_DPhi0 / self.silica_DPhi0 * self.silica_n2 * self.l_silica / 0.001
                        )  # [m2/W] 0.001 stands for beam path in 1-mm cuvette

                        # calculate errors
                        match from_what:
                            case "from_geometry":  # fit_manually
                                pass  # don't calculate errors, they are not known
                            case "from_autofit":  # fit_automatically
                                try:
                                    self.solventCA_zeroLevel, self.solventCA_zeroLevelError, self.solventCA_zeroLevelPrecision = error_rounding(
                                        self.solventCA_minimizerResult.params["Zero"].value, self.solventCA_minimizerResult.params["Zero"].stderr
                                    )
                                    self.solventCA_centerPoint, self.solventCA_centerPointError, self.solventCA_centerPointPrecision = error_rounding(
                                        self.solventCA_minimizerResult.params["Center"].value, self.solventCA_minimizerResult.params["Center"].stderr
                                    )
                                    self.solvent_DPhi0, self.solvent_DPhi0Error, self.solvent_DPhi0Precision = error_rounding(
                                        self.solventCA_minimizerResult.params["DPhi0"].value, self.solventCA_minimizerResult.params["DPhi0"].stderr
                                    )
                                    if self.solventCA_minimizerResult.params["Beamwaist"].vary is True:
                                        self.solvent_beamwaist, self.solvent_beamwaistError, self.solvent_beamwaistPrecision = error_rounding(
                                            self.solventCA_minimizerResult.params["Beamwaist"].value,
                                            self.solventCA_minimizerResult.params["Beamwaist"].stderr,
                                        )
                                    else:
                                        self.solvent_beamwaist, self.solvent_beamwaistError, self.solvent_beamwaistPrecision = (
                                            self.silica_beamwaist,
                                            self.silica_beamwaistError,
                                            self.silica_beamwaistPrecision,
                                        )

                                    self.solvent_rayleighLengthError = (
                                        np.pi / 2 / self.lda * self.solventCA_minimizerResult.params["Beamwaist"].value * self.solvent_beamwaistError
                                    )  # [m] Rayleigh length
                                    self.solvent_rayleighLength, self.solvent_rayleighLengthError, self.solvent_rayleighLengthPrecision = (
                                        error_rounding(self.solvent_rayleighLength, self.solvent_rayleighLengthError)
                                    )

                                    self.solvent_n2Error = (
                                        self.solvent_DPhi0Error / self.silica_DPhi0 * self.silica_n2 * self.l_silica / 0.001
                                        + self.solvent_DPhi0 * self.silica_DPhi0Error / self.silica_DPhi0**2 * self.silica_n2 * self.l_silica / 0.001
                                    )
                                    self.solvent_n2, self.solvent_n2Error, self.solvent_n2Precision = error_rounding(
                                        self.solvent_n2 * 1e13, self.solvent_n2Error * 1e13
                                    )
                                    # recover original units of m2/W
                                    self.solvent_n2 *= 1e-13
                                    self.solvent_n2Error *= 1e-13
                                except Exception:
                                    logging.error(traceback.format_exc())
                                    self.showdialog("Error", "The fit didn't converge. Try using different initial parameters.")
                    case "OA":
                        pass
                match stype:
                    case "CA":
                        match from_what:
                            case "from_geometry":
                                self.solventCA_zeroLevel = self.solventCA_zeroLevel_slider.value() / 100
                                self.solventCA_centerPoint = np.round(self.solventCA_centerPoint_slider.value() - self.solvent_nop / 2)
                                if on_data_load:
                                    self.solvent_DPhi0, self.solvent_beamwaist, self.solvent_rayleighLength = self.params_from_geometry(ftype, "CA")
                                else:
                                    self.solvent_DPhi0 = self.solventCA_DPhi0_slider.value() * MAX_DPHI0 / self.solventCA_DPhi0_slider.maximum()
                                    self.solvent_rayleighLength = (
                                        self.solventCA_RayleighLength_slider.value()
                                        * self.z_range
                                        / 2
                                        / self.solventCA_RayleighLength_slider.maximum()
                                    )
                                    self.solvent_beamwaist = np.sqrt(self.solvent_rayleighLength * self.lda / np.pi)

                                if ftype == "solvent" and self.solventCA_customBeamwaist_checkBox.isChecked() is False:
                                    self.solvent_beamwaist, self.solvent_rayleighLength = self.silica_beamwaist, self.silica_rayleighLength
                            case "from_autofit":
                                pass

                        self.solvent_n2 = (
                            self.solvent_DPhi0 / self.silica_DPhi0 * self.silica_n2 * self.l_silica / 0.001
                        )  # [m2/W] 0.001 stands for beam path in 1-mm cuvette

                        # calculate errors
                        match from_what:
                            case "from_geometry":  # fit_manually
                                pass  # don't calculate errors, they are not known
                            case "from_autofit":  # fit_automatically
                                try:
                                    self.solventCA_zeroLevel, self.solventCA_zeroLevelError, self.solventCA_zeroLevelPrecision = error_rounding(
                                        self.solventCA_minimizerResult.params["Zero"].value, self.solventCA_minimizerResult.params["Zero"].stderr
                                    )
                                    self.solventCA_centerPoint, self.solventCA_centerPointError, self.solventCA_centerPointPrecision = error_rounding(
                                        self.solventCA_minimizerResult.params["Center"].value, self.solventCA_minimizerResult.params["Center"].stderr
                                    )
                                    self.solvent_DPhi0, self.solvent_DPhi0Error, self.solvent_DPhi0Precision = error_rounding(
                                        self.solventCA_minimizerResult.params["DPhi0"].value, self.solventCA_minimizerResult.params["DPhi0"].stderr
                                    )
                                    if self.solventCA_minimizerResult.params["Beamwaist"].vary is True:
                                        self.solvent_beamwaist, self.solvent_beamwaistError, self.solvent_beamwaistPrecision = error_rounding(
                                            self.solventCA_minimizerResult.params["Beamwaist"].value,
                                            self.solventCA_minimizerResult.params["Beamwaist"].stderr,
                                        )
                                    else:
                                        self.solvent_beamwaist, self.solvent_beamwaistError, self.solvent_beamwaistPrecision = (
                                            self.silica_beamwaist,
                                            self.silica_beamwaistError,
                                            self.silica_beamwaistPrecision,
                                        )

                                    self.solvent_rayleighLengthError = (
                                        np.pi / 2 / self.lda * self.solventCA_minimizerResult.params["Beamwaist"].value * self.solvent_beamwaistError
                                    )  # [m] Rayleigh length
                                    self.solvent_rayleighLength, self.solvent_rayleighLengthError, self.solvent_rayleighLengthPrecision = (
                                        error_rounding(self.solvent_rayleighLength, self.solvent_rayleighLengthError)
                                    )

                                    self.solvent_n2Error = (
                                        self.solvent_DPhi0Error / self.silica_DPhi0 * self.silica_n2 * self.l_silica / 0.001
                                        + self.solvent_DPhi0 * self.silica_DPhi0Error / self.silica_DPhi0**2 * self.silica_n2 * self.l_silica / 0.001
                                    )
                                    self.solvent_n2, self.solvent_n2Error, self.solvent_n2Precision = error_rounding(
                                        self.solvent_n2 * 1e13, self.solvent_n2Error * 1e13
                                    )
                                    # recover original units of m2/W
                                    self.solvent_n2 *= 1e-13
                                    self.solvent_n2Error *= 1e-13
                                except Exception:
                                    logging.error(traceback.format_exc())
                                    self.showdialog("Error", "The fit didn't converge. Try using different initial parameters.")
                    case "OA":
                        pass
            case "sample":
                if self.sampleCA_fittingLine_drawn is True:
                    pass

                self.sampleCA_zeroLevel = self.sampleCA_zeroLevel_slider.value() / 100
                self.sampleCA_centerPoint = self.sampleCA_centerPoint_slider.value() - 50

                if self.sampleOA_fittingLine_drawn is True:
                    pass

                self.sampleOA_zeroLevel = self.sampleOA_zeroLevel_slider.value() / 100
                self.sampleOA_centerPoint = self.sampleOA_centerPoint_slider.value() - 50

    def params_from_geometry(self, ftype, stype):
        """Using equations given by van Stryland and Sheik-Bahae in
        http://www.phys.unm.edu/msbahae/publications/z-scan.pdf
        return the fitting curve physical interpretation."""
        match ftype:
            case "silica":
                if stype == "CA":
                    dataset = self.silica_data_set[0], self.silica_data_set[1]
                else:
                    dataset = self.silica_data_set[0], self.silica_data_set[3]
            case "solvent":
                if stype == "CA":
                    dataset = self.solvent_data_set[0], self.solvent_data_set[1]
                else:
                    dataset = self.solvent_data_set[0], self.solvent_data_set[3]
            case "sample":
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
                deltaZpv = abs(fit_x[ymax_pos] - fit_x[ymin_pos]) * 1e-3  # [m]
                rayleighLength = deltaZpv / 1.7  # [m] Rayleigh length
                beamwaist = np.sqrt(rayleighLength * self.lda / np.pi)  # [m] beam radius in focal point
                return deltaPhi0, beamwaist, rayleighLength

            except TypeError:
                print("Something is wrong while interpreting the closed aperture curve.")

        elif stype == "OA":
            pass

    def set_sliders_positions(self, ftype, stype):
        """Called by `self.fit_automatically`.\n
        This method recalculates value to reflect given slider position.
        """
        match ftype:
            case "silica":
                # silica CA sliders
                # Because change in slider value triggers the 'fit_manually' method, it has to be disabled.
                self.silica_RayleighLength_slider.valueChanged.disconnect()
                self.silica_centerPoint_slider.valueChanged.disconnect()
                self.silica_zeroLevel_slider.valueChanged.disconnect()
                self.silica_DPhi0_slider.valueChanged.disconnect()

                self.silica_RayleighLength_slider.setValue(
                    int(round(self.silica_rayleighLength * self.silica_RayleighLength_slider.maximum() / (self.z_range / 2)))
                )
                self.silica_centerPoint_slider.setValue(int(round(self.silica_centerPoint + self.silica_centerPoint_slider.maximum() / 2)))
                self.silica_zeroLevel_slider.setValue(int(round(self.silica_zeroLevel * 100)))
                self.silica_DPhi0_slider.setValue(int(round(self.silica_DPhi0 / np.pi * self.silica_DPhi0_slider.maximum())))
                self.silica_DPhi0_slider.setValue(int(round(self.silica_DPhi0 / np.pi * self.silica_DPhi0_slider.maximum())))

                # And now when all is updated by the 'fit_automatically', reconnect the sliders to their slots
                self.fitting_tab_signals()

            case "solvent":
                if stype == "CA":
                    # solvent CA sliders
                    # Because change in slider value triggers the 'fit_manually' method, it has to be disabled.
                    self.solventCA_RayleighLength_slider.valueChanged.disconnect()
                    self.solventCA_centerPoint_slider.valueChanged.disconnect()
                    self.solventCA_zeroLevel_slider.valueChanged.disconnect()
                    self.solventCA_DPhi0_slider.valueChanged.disconnect()

                    if self.solventCA_customBeamwaist_checkBox.isChecked() is False:
                        self.solventCA_RayleighLength_slider.setValue(self.silica_RayleighLength_slider.value())
                    else:
                        self.solventCA_RayleighLength_slider.setValue(
                            int(round(self.solvent_rayleighLength * self.solventCA_RayleighLength_slider.maximum() / (self.z_range / 2)))
                        )
                    self.solventCA_centerPoint_slider.setValue(
                        int(round(self.solventCA_centerPoint + self.solventCA_centerPoint_slider.maximum() / 2))
                    )
                    self.solventCA_zeroLevel_slider.setValue(int(round(self.solventCA_zeroLevel * 100)))
                    self.solventCA_DPhi0_slider.setValue(int(round(self.solvent_DPhi0 / np.pi * self.solventCA_DPhi0_slider.maximum())))
                    self.solventCA_DPhi0_slider.setValue(int(round(self.solvent_DPhi0 / np.pi * self.solventCA_DPhi0_slider.maximum())))

                    # And now when all is updated by the 'fit_automatically', reconnect the sliders to their slots
                    self.solventCA_RayleighLength_slider.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="CA"))
                    self.solventCA_centerPoint_slider.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="CA"))
                    self.solventCA_zeroLevel_slider.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="CA"))
                    self.solventCA_DPhi0_slider.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="CA"))

                elif stype == "OA":
                    # solvent OA sliders
                    # Because change in slider value triggers the 'fit_manually' method, it has to be disabled.
                    self.solventOA_centerPoint_slider.valueChanged.disconnect()
                    self.solventOA_zeroLevel_slider.valueChanged.disconnect()
                    self.solventOA_T_slider.valueChanged.disconnect()

                    self.solventOA_centerPoint_slider.setValue(int(round(self.solventOA_centerPoint + 50)))
                    self.solventOA_zeroLevel_slider.setValue(int(round(self.solventOA_zeroLevel * 100)))
                    self.solventOA_T_slider.setValue(int(round(self.solventOA_T * self.solventOA_T_slider.maximum() / SOLVENT_T_SLIDER_MAX)))

                    # And now when all is updated by the 'fit_automatically', reconnect the sliders to their slots
                    self.solventOA_centerPoint_slider.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="OA"))
                    self.solventOA_zeroLevel_slider.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="OA"))
                    self.solventOA_T_slider.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="OA"))

            case "sample":
                pass

    def data_display(self, data_set, ftype: str):
        """Called by data_loader(). Displays:\n
        1) CA column divided by OA column to remove influence of OA on CA as 'CA'
        2) OA column as 'OA'\n
        to separate data processing for n2 and beta nonlinear coefficients."""
        self.get_general_parameters()

        def basic_data_manipulation(data) -> Tuple[NDArray, list, list, list]:
            """Reference division, separation of OA signal from CA fluctuation signal, normalization

            Args:
                data (any): four-column data containing datapoints, CA data, Reference data and OA data

            Returns:
                Tuple: Data ready for separated determination of n2 and beta nonlinear coefficients
            """
            positions, ca_data0, ref_data0, oa_data0 = data
            # Divide data by reference
            ca_data = [ca / ref for ca, ref in zip(ca_data0, ref_data0)]
            # ref_data = [ref/ref for ref in ref_data0]
            oa_data = [oa / ref for oa, ref in zip(oa_data0, ref_data0)]

            # centralize positions and normalize by z-scan range
            nop = len(positions)
            positions = np.array([(self.z_range * zz / nop - self.z_range / 2) * 1000 for zz in range(nop)])  # express in mm

            def normalize(data):  # Normalize to 1 (edge of the data as reference)
                zero_lvl = (np.mean(data[0:10]) + np.mean(data[len(data) - 10 :])) / 2
                return data / zero_lvl

            ca_data = normalize(ca_data)
            oa_data = normalize(oa_data)

            return positions, ca_data, ref_data0, oa_data

        self.positions, self.ca_data, self.ref_data, self.oa_data = basic_data_manipulation(data_set)

        # Create reference-corrected 'data_set' variable for each ftype for further reference
        match ftype:
            case "silica":
                self.silica_data_set = [self.positions, self.ca_data, self.ref_data, self.oa_data]
                self.silica_nop = len(self.positions)
            case "solvent":
                self.solvent_data_set = [self.positions, self.ca_data, self.ref_data, self.oa_data]
                self.solvent_nop = len(self.positions)
            case "sample":
                self.sample_data_set = [self.positions, self.ca_data, self.ref_data, self.oa_data]
                self.sample_nop = len(self.positions)

        # Display data on proper figures
        # Use reference-corrected data for fitting
        match ftype:
            case "silica":
                self.update_datafitting_plotlimits(self.silica_data_set, ftype=ftype)

            case "solvent":
                self.update_datafitting_plotlimits(self.solvent_data_set, ftype=ftype)

            case "sample":
                self.update_datafitting_plotlimits(self.sample_data_set, ftype=ftype)

    def set_new_positions(self):
        """Takes current values in 'General Parameters' from GUI and updates data_set[0] (positions) to meet new Z-scan range.\n
        Updates fitting plots limits to visualize the change.
        """
        self.get_general_parameters()

        if hasattr(self, "silica_data_set"):
            self.silica_data_set[0] = np.array(
                [(self.z_range * zz / self.silica_nop - self.z_range / 2) * 1000 for zz in range(self.silica_nop)]
            )  # [mm] update positions with newly-read Z-scan range value
            self.update_datafitting_plotlimits(self.silica_data_set, ftype="silica")

        if hasattr(self, "solvent_data_set"):
            self.solvent_data_set[0] = np.array(
                [(self.z_range * zz / self.solvent_nop - self.z_range / 2) * 1000 for zz in range(self.solvent_nop)]
            )  # [mm] update positions with newly-read Z-scan range value
            self.update_datafitting_plotlimits(self.solvent_data_set, ftype="solvent")

        if hasattr(self, "sample_data_set"):
            self.sample_data_set[0] = np.array(
                [(self.z_range * zz / self.sample_nop - self.z_range / 2) * 1000 for zz in range(self.sample_nop)]
            )  # [mm] update positions with newly-read Z-scan range value
            self.update_datafitting_plotlimits(self.sample_data_set, ftype="sample")

    def update_datafitting_plotlimits(self, data_set: list, ftype: str) -> None:
        """Updates limits of the data fitting plots for given `ftype` (silica, solvent, sample)

        Args:
            data_set (list): four-column reference-corrected data
            ftype (str): parameter holding information of sample type (silica, solvent, sample)
        """
        self.get_general_parameters()  # Ensures working with currently typed-in General Parameters values from GUI
        padding_vertical = 0.01

        def set_limits(axes, line, direction, padding):
            if direction == "vertical":
                axes.set_ylim(top=np.max(line.get_ydata()) * (1 + padding), bottom=np.min(line.get_ydata()) * (1 - padding))

        match ftype:
            case "silica":
                # Closed aperture
                line = self.silica_figure.axes.get_lines()[0]
                line.set_xdata(data_set[0])
                # line.set_ydata(ca_data)
                line.set_ydata([ca / oa for ca, oa in zip(data_set[1], data_set[3])])

                self.silica_figure.axes.set_xlim(left=-self.z_range / 2 * 1000, right=self.z_range / 2 * 1000)  # displayed in mm
                set_limits(self.silica_figure.axes, line, "vertical", padding_vertical)

                # Open aperture
                line = self.silicaOA_figure.axes.get_lines()[0]
                line.set_xdata(data_set[0])
                line.set_ydata(data_set[3])

                self.silicaOA_figure.axes.set_xlim(left=-self.z_range / 2 * 1000, right=self.z_range / 2 * 1000)  # displayed in mm
                set_limits(self.silicaOA_figure.axes, line, "vertical", padding_vertical)
                # self.silicaOA_figure.axes.set_ylim(top=np.max(line.get_ydata())*(1+margin_vertical),bottom=np.min(line.get_ydata())*(1-margin_vertical))

                # Update
                self.silica_figure.axes.relim()
                self.silicaOA_figure.axes.relim()

                self.silica_figure.axes.autoscale_view()
                self.silicaOA_figure.axes.autoscale_view()

                self.silica_figure.draw_idle()
                self.silicaOA_figure.draw_idle()

            case "solvent":
                # Closed aperture
                line = self.solventCA_figure.axes.get_lines()[0]
                line.set_xdata(data_set[0])
                # line.set_ydata(ca_data)
                line.set_ydata([ca / oa for ca, oa in zip(data_set[1], data_set[3])])

                self.solventCA_figure.axes.set_xlim(left=-self.z_range / 2 * 1000, right=self.z_range / 2 * 1000)  # displayed in mm
                set_limits(self.solventCA_figure.axes, line, "vertical", padding_vertical)
                # self.solventCA_figure.axes.set_ylim(top=np.max(line.get_ydata())*(1+margin_vertical),bottom=np.min(line.get_ydata())*(1-margin_vertical))

                # Open aperture
                line = self.solventOA_figure.axes.get_lines()[0]
                line.set_xdata(data_set[0])
                line.set_ydata(data_set[3])

                self.solventOA_figure.axes.set_xlim(left=-self.z_range / 2 * 1000, right=self.z_range / 2 * 1000)  # displayed in mm
                set_limits(self.solventOA_figure.axes, line, "vertical", padding_vertical)
                # self.solventOA_figure.axes.set_ylim(top=np.max(line.get_ydata())*(1+margin_vertical),bottom=np.min(line.get_ydata())*(1-margin_vertical))

                # Update
                self.solventCA_figure.axes.relim()
                self.solventOA_figure.axes.relim()

                self.solventCA_figure.axes.autoscale_view()
                self.solventOA_figure.axes.autoscale_view()

                self.solventCA_figure.draw_idle()
                self.solventOA_figure.draw_idle()

            case "sample":
                # Closed aperture
                line = self.sampleCA_figure.axes.get_lines()[0]
                line.set_xdata(data_set[0])
                # line.set_ydata(ca_data)
                line.set_ydata([ca / oa for ca, oa in zip(data_set[1], data_set[3])])

                self.sampleCA_figure.axes.set_xlim(left=-self.z_range / 2 * 1000, right=self.z_range / 2 * 1000)  # displayed in mm
                set_limits(self.sampleCA_figure.axes, line, "vertical", padding_vertical)
                # self.sampleCA_figure.axes.set_ylim(top=np.max(line.get_ydata())*1.04,bottom=np.min(line.get_ydata())*0.96)

                # Open aperture
                line = self.sampleOA_figure.axes.get_lines()[0]
                line.set_xdata(data_set[0])
                line.set_ydata(data_set[3])

                self.sampleOA_figure.axes.set_xlim(left=-self.z_range / 2 * 1000, right=self.z_range / 2 * 1000)  # displayed in mm
                set_limits(self.sampleOA_figure.axes, line, "vertical", padding_vertical)
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
            case "silica":
                if stype == "CA":
                    filter_size = self.silica_filterSize_slider.value()
                elif stype == "OA":
                    filter_size = self.silicaOA_filterSize_slider.value()
            case "solvent":
                if stype == "CA":
                    filter_size = self.solventCA_filterSize_slider.value()
                elif stype == "OA":
                    filter_size = self.solventOA_filterSize_slider.value()
            case "sample":
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
                    case "silica":
                        if stype == "CA":
                            line = self.silica_figure.axes.get_lines()[0]
                            line.set_ydata(filtered_y)

                            self.silica_figure.axes.relim()
                            self.silica_figure.axes.autoscale_view()
                            self.silica_figure.draw_idle()
                        elif stype == "OA":
                            pass

                    case "solvent":
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

                    case "sample":
                        pass

            except ValueError:
                self.showdialog("Error", "Possibly too few data points!\nMinimum required is 16 datapoints.")

    def enable_cursors(self, ftype: str, stype: str) -> None:
        match ftype:
            case "silica":
                if stype == "CA":
                    if self.silica_fixROI_checkBox.isChecked() is True:
                        self.silica_cursor_positioner = BlittedCursor(self.silica_figure.axes, color="magenta", linewidth=2)
                        self.on_mouse_move = self.silica_figure.mpl_connect("motion_notify_event", self.silica_cursor_positioner.on_mouse_move)
                        self.on_mouse_click = self.silica_figure.mpl_connect(
                            "button_press_event", lambda event: self.collect_cursor_clicks(event, ftype)
                        )
                    else:
                        try:
                            self.silica_cursor_positioner.vertical_line.remove()
                            self.silica_cursor_positioner.horizontal_line.remove()
                        except (AttributeError, ValueError):
                            print("Specified cross-hair doesn't exist.")
                        finally:
                            self.silica_figure.mpl_disconnect(self.on_mouse_move)
                            self.silica_figure.mpl_disconnect(self.on_mouse_click)
                            self.silica_cursorPositions = []
                        try:
                            self.silica_verline1.remove()
                            self.silica_verline2.remove()
                        except (AttributeError, ValueError):
                            print("Specified cross-hair doesn't exist.")
                        finally:
                            self.silica_figure.draw_idle()

                elif stype == "OA":
                    pass
            case "solvent":
                if stype == "CA":
                    if self.solventCA_fixROI_checkBox.isChecked() is True:
                        self.solventCA_cursor_positioner = BlittedCursor(self.solventCA_figure.axes, color="magenta", linewidth=2)
                        self.on_mouse_move = self.solventCA_figure.mpl_connect("motion_notify_event", self.solventCA_cursor_positioner.on_mouse_move)
                        self.on_mouse_click = self.solventCA_figure.mpl_connect(
                            "button_press_event", lambda event: self.collect_cursor_clicks(event, ftype, stype)
                        )
                    else:
                        try:
                            self.solventCA_cursor_positioner.vertical_line.remove()
                            self.solventCA_cursor_positioner.horizontal_line.remove()
                        except (AttributeError, ValueError):
                            print("Specified cross-hair doesn't exist.")
                        finally:
                            self.solventCA_figure.mpl_disconnect(self.on_mouse_move)
                            self.solventCA_figure.mpl_disconnect(self.on_mouse_click)
                            self.solventCA_cursorPositions = []
                        try:
                            self.solventCA_verline1.remove()
                            self.solventCA_verline2.remove()
                        except (AttributeError, ValueError):
                            print("Specified cross-hair doesn't exist.")
                        finally:
                            self.solventCA_figure.draw_idle()

                elif stype == "OA":
                    if self.solventOA_fixROI_checkBox.isChecked() is True:
                        self.solventOA_cursor_positioner = BlittedCursor(self.solventOA_figure.axes, color="magenta", linewidth=2)
                        self.on_mouse_move = self.solventOA_figure.mpl_connect("motion_notify_event", self.solventOA_cursor_positioner.on_mouse_move)
                        self.on_mouse_click = self.solventOA_figure.mpl_connect(
                            "button_press_event", lambda event: self.collect_cursor_clicks(event, ftype, stype)
                        )
                    else:
                        try:
                            self.solventOA_cursor_positioner.vertical_line.remove()
                            self.solventOA_cursor_positioner.horizontal_line.remove()
                        except (AttributeError, ValueError):
                            print("Specified cross-hair doesn't exist.")
                        finally:
                            self.solventOA_figure.mpl_disconnect(self.on_mouse_move)
                            self.solventOA_figure.mpl_disconnect(self.on_mouse_click)
                            self.solventOA_cursorPositions = []
                        try:
                            self.solventOA_verline1.remove()
                            self.solventOA_verline2.remove()
                        except (AttributeError, ValueError):
                            print("Specified cross-hair doesn't exist.")
                        finally:
                            self.solventOA_figure.draw_idle()

            case "sample":
                if stype == "CA":
                    self.sampleCA_cursor_positioner = BlittedCursor(self.sampleCA_figure.axes, color="magenta", linewidth=2)
                elif stype == "OA":
                    self.sampleOA_cursor_positioner = BlittedCursor(self.sampleOA_figure.axes, color="magenta", linewidth=2)

    def collect_cursor_clicks(self, event, ftype: str, stype="") -> None:
        x, y = event.xdata, event.ydata
        match ftype:
            case "silica":
                if not hasattr(self, "silica_cursorPositions"):
                    self.silica_cursorPositions = []
                if len(self.silica_cursorPositions) < 2:
                    self.silica_cursorPositions.append((x, y))
                    if len(self.silica_cursorPositions) == 1:
                        self.silica_verline1 = self.silica_figure.axes.axvline(x, color="orange", linewidth=2)
                    else:
                        self.silica_verline2 = self.silica_figure.axes.axvline(x, color="orange", linewidth=2)
                else:
                    self.silica_cursorPositions = []
                    self.silica_cursorPositions.append((x, y))
                    self.silica_verline1.set_xdata([x])
                    self.silica_verline2.set_xdata([None])

                self.silica_figure.draw_idle()

            case "solvent":
                if stype == "CA":
                    if not hasattr(self, "solventCA_cursorPositions"):
                        self.solventCA_cursorPositions = []
                    if len(self.solventCA_cursorPositions) < 2:
                        self.solventCA_cursorPositions.append((x, y))
                        if len(self.solventCA_cursorPositions) == 1:
                            self.solventCA_verline1 = self.solventCA_figure.axes.axvline(x, color="orange", linewidth=2)
                        else:
                            self.solventCA_verline2 = self.solventCA_figure.axes.axvline(x, color="orange", linewidth=2)
                    else:
                        self.solventCA_cursorPositions = []
                        self.solventCA_cursorPositions.append((x, y))
                        self.solventCA_verline1.set_xdata(x)
                        self.solventCA_verline2.set_xdata(None)

                    self.solventCA_figure.draw_idle()

                elif stype == "OA":
                    if not hasattr(self, "solventOA_cursorPositions"):
                        self.solventOA_cursorPositions = []
                    if len(self.solventOA_cursorPositions) < 2:
                        self.solventOA_cursorPositions.append((x, y))
                        if len(self.solventOA_cursorPositions) == 1:
                            self.solventOA_verline1 = self.solventOA_figure.axes.axvline(x, color="orange", linewidth=2)
                        else:
                            self.solventOA_verline2 = self.solventOA_figure.axes.axvline(x, color="orange", linewidth=2)
                    else:
                        self.solventOA_cursorPositions = []
                        self.solventOA_cursorPositions.append((x, y))
                        self.solventOA_verline1.set_xdata(x)
                        self.solventOA_verline2.set_xdata(None)

                    self.solventOA_figure.draw_idle()

            case "sample":
                pass

    def draw_fitting_line(self, ftype: str, stype: str) -> None:
        match ftype:
            case "silica":
                if self.silica_fittingLine_drawn is False:
                    (self.silica_fitting_line_ca,) = self.silica_figure.axes.plot(self.silica_data_set[0], self.result, "r")
                    self.silica_figure.draw_idle()

                    self.silica_fittingLine_drawn = True
                else:
                    if hasattr(self, "silica_data_set"):
                        self.silica_fitting_line_ca.set_xdata(self.silica_data_set[0])
                        self.silica_fitting_line_ca.set_ydata(self.result)

                        self.silica_figure.axes.set_xlim(left=-self.z_range / 2 * 1000, right=self.z_range / 2 * 1000)  # displayed in mm
                        self.silica_figure.axes.relim()
                        self.silica_figure.axes.autoscale_view()
                        self.silica_figure.draw_idle()

            case "solvent":
                if stype == "CA":
                    if self.solventCA_fittingLine_drawn is False:
                        (self.solvent_fitting_line_ca,) = self.solventCA_figure.axes.plot(self.solvent_data_set[0], self.result, "r")
                        self.solventCA_figure.draw_idle()

                        self.solventCA_fittingLine_drawn = True
                    else:
                        if hasattr(self, "solvent_data_set"):
                            self.solvent_fitting_line_ca.set_xdata(self.solvent_data_set[0])
                            self.solvent_fitting_line_ca.set_ydata(self.result)

                            self.solventCA_figure.axes.set_xlim(left=-self.z_range / 2 * 1000, right=self.z_range / 2 * 1000)  # displayed in mm
                            self.solventCA_figure.axes.relim()
                            self.solventCA_figure.axes.autoscale_view()
                            self.solventCA_figure.draw_idle()

                elif stype == "OA":
                    if self.solventOA_fittingLine_drawn is False:
                        (self.solvent_fitting_line_oa,) = self.solventOA_figure.axes.plot(self.solvent_data_set[0], self.result, "r")
                        self.solventOA_figure.draw_idle()

                        self.solventOA_fittingLine_drawn = True
                    else:
                        if hasattr(self, "solvent_data_set"):
                            self.solvent_fitting_line_oa.set_xdata(self.solvent_data_set[0])
                            self.solvent_fitting_line_oa.set_ydata(self.result)

                            self.solventOA_figure.axes.set_xlim(left=-self.z_range / 2 * 1000, right=self.z_range / 2 * 1000)  # displayed in mm
                            self.solventOA_figure.axes.relim()
                            self.solventOA_figure.axes.autoscale_view()
                            self.solventOA_figure.draw_idle()

            case "sample":
                if stype == "CA":
                    if self.sampleCA_fittingLine_drawn is False:
                        (self.sample_fitting_line_ca,) = self.sampleCA_figure.axes.plot(self.sample_data_set[0], self.result, "r")
                        self.sampleCA_figure.draw_idle()

                        self.sampleCA_fittingLine_drawn = True
                    else:
                        if hasattr(self, "sample_data_set"):
                            self.sample_fitting_line_ca.set_xdata(self.sample_data_set[0])
                            self.sample_fitting_line_ca.set_ydata(self.result)

                            self.sampleCA_figure.axes.set_xlim(left=-self.z_range / 2 * 1000, right=self.z_range / 2 * 1000)  # displayed in mm
                            self.sampleCA_figure.axes.relim()
                            self.sampleCA_figure.axes.autoscale_view()
                            self.sampleCA_figure.draw_idle()

                elif stype == "OA":
                    if self.sampleOA_fittingLine_drawn is False:
                        (self.sample_fitting_line_oa,) = self.sampleOA_figure.axes.plot(self.sample_data_set[0], self.result, "r")
                        self.sampleOA_figure.draw_idle()

                        self.sampleOA_fittingLine_drawn = True
                    else:
                        if hasattr(self, "sample_data_set"):
                            self.sample_fitting_line_oa.set_xdata(self.sample_data_set[0])
                            self.sample_fitting_line_oa.set_ydata(self.result)

                            self.sampleOA_figure.axes.set_xlim(left=-self.z_range / 2 * 1000, right=self.z_range / 2 * 1000)  # displayed in mm
                            self.sampleOA_figure.axes.relim()
                            self.sampleOA_figure.axes.autoscale_view()
                            self.sampleOA_figure.draw_idle()

    def fit_automatically(self, ftype: str, stype: str):
        self.get_general_parameters()
        self.get_curve_interpretation(ftype, stype, "from_geometry")

        match ftype:
            case "silica":
                if stype == "CA":
                    # Datapoints for the curve to be fitted to
                    line = self.silica_figure.axes.get_lines()[0]
                else:
                    return

                line_data = line.get_data()

                self.silica_calculation = Fitting(
                    self,
                    self.silica_curves,
                    self.silica_DPhi0,
                    self.silica_beamwaist,
                    self.silica_zeroLevel,
                    self.silica_centerPoint,
                    self.silica_nop,
                    line_data[1],
                )
                self.silica_calculation = Fitting(
                    self,
                    self.silica_curves,
                    self.silica_DPhi0,
                    self.silica_beamwaist,
                    self.silica_zeroLevel,
                    self.silica_centerPoint,
                    self.silica_nop,
                    line_data[1],
                )
                minimizer_result, self.result = self.silica_calculation.automatic(self.z_range, ftype, stype, line_data)
                # Make sure that the result doesn't contain NoneTypes
                for key in minimizer_result.params.keys():
                    if minimizer_result.params[key].stderr is None:
                        minimizer_result.params[key].stderr = 0
                self.silica_minimizerResult = minimizer_result
                self.silica_autofit_done = True
                self.draw_fitting_line(ftype, stype)
                self.get_curve_interpretation(ftype, stype, "from_autofit")
                self.set_fit_summary(ftype, stype, caller="auto")
                # Use these exact number from error-corrected fit parameters to set sliders to their positions
                self.set_sliders_positions(ftype, stype)

            case "solvent":
                if self.silica_autofit_done is False:
                    self.showdialog("Info", "Fit silica first.")
                    return

                if stype == "CA":
                    # Data to be fitted
                    line = self.solventCA_figure.axes.get_lines()[0]
                elif stype == "OA":
                    line = self.solventOA_figure.axes.get_lines()[0]

                line_data = line.get_data()

                self.solvent_nop = len(self.solvent_data_set[0])
                self.solvent_data_set[0] = np.array(
                    [(self.z_range * zz / self.solvent_nop - self.z_range / 2) * 1000 for zz in range(self.solvent_nop)]
                )  # [mm] update positions with newly-read Z-scan range value
                self.solvent_nop = len(self.solvent_data_set[0])
                self.solvent_data_set[0] = np.array(
                    [(self.z_range * zz / self.solvent_nop - self.z_range / 2) * 1000 for zz in range(self.solvent_nop)]
                )  # [mm] update positions with newly-read Z-scan range value

                self.solvent_calculation = Fitting(
                    self,
                    self.solvent_curves,
                    self.solvent_DPhi0,
                    self.solvent_beamwaist,
                    self.solventCA_zeroLevel,
                    self.solventCA_centerPoint,
                    self.solvent_nop,
                    line_data[1],
                )
                self.solvent_calculation = Fitting(
                    self,
                    self.solvent_curves,
                    self.solvent_DPhi0,
                    self.solvent_beamwaist,
                    self.solventCA_zeroLevel,
                    self.solventCA_centerPoint,
                    self.solvent_nop,
                    line_data[1],
                )
                minimizer_result, self.result = self.solvent_calculation.automatic(self.z_range, ftype, stype, line_data)
                match stype:
                    case "CA":
                        self.solventCA_minimizerResult = minimizer_result
                        self.solventCA_autofit_done = True
                    case "OA":
                        self.solventOA_minimizerResult = minimizer_result
                        self.solventOA_autofit_done = True
                self.draw_fitting_line(ftype, stype)
                self.get_curve_interpretation(ftype, stype, "from_autofit")
                match stype:
                    case "CA":
                        self.solventCA_minimizerResult = minimizer_result
                        self.solventCA_autofit_done = True
                    case "OA":
                        self.solventOA_minimizerResult = minimizer_result
                        self.solventOA_autofit_done = True
                self.draw_fitting_line(ftype, stype)
                self.get_curve_interpretation(ftype, stype, "from_autofit")
                self.set_fit_summary(ftype, stype, caller="auto")
                # Use these exact number from error-corrected fit parameters to set sliders to their positions
                self.set_sliders_positions(ftype, stype)
                # Use these exact number from error-corrected fit parameters to set sliders to their positions
                self.set_sliders_positions(ftype, stype)
            case "sample":
                if self.silica_autofit_done is False:
                    self.showdialog("Info", "Fit silica first.")

                    if self.solventCA_autofit_done is False:
                        self.showdialog("Info", "Fit solvent first.")
                        return
                    else:
                        return

                pass

    def fit_manually(self, ftype: str, stype: str) -> None:
        """Triggered by loading the data (from experiment or from file) or by fitting sliders value change.

        Args:
            ftype (str): `silica`, `solvent`, `sample`
            stype (str): `CA`, `OA`
        """
        # Retrieve "General parameters"
        self.get_general_parameters()
        self.get_curve_interpretation(ftype, stype, "from_geometry")

        match ftype:
            case "silica":
                if stype == "CA":
                    # Datapoints to fit the curve to
                    line = self.silica_figure.axes.get_lines()[0]
                    line_data = line.get_ydata()

                    self.silica_curves = Integration(
                        self,
                        SILICA_BETA,
                        self.silica_n2,
                        self.silica_DPhi0,
                        self.silica_data_set[0],
                        self.d0,
                        self.ra,
                        self.lda,
                        self.silica_beamwaist,
                        N_COMPONENTS,
                        INTEGRATION_STEPS,
                    )
                    self.silica_calculation = Fitting(
                        self,
                        self.silica_curves,
                        self.silica_DPhi0,
                        self.silica_beamwaist,
                        self.silica_zeroLevel,
                        self.silica_centerPoint,
                        len(self.silica_data_set[0]),
                        line_data,
                    )
                    self.result = self.silica_calculation.manual(
                        self.silica_zeroLevel, self.silica_centerPoint, self.silica_DPhi0, self.silica_beamwaist, self.z_range
                    )
                    self.silica_autofit_done = False
            case "solvent":
                # Data to be fitted
                if stype == "CA":
                    line = self.solventCA_figure.axes.get_lines()[0]
                    line_data = line.get_ydata()

                    self.solvent_curves = Integration(
                        self,
                        0,
                        self.solvent_n2,
                        self.solvent_DPhi0,
                        self.solvent_data_set[0],
                        self.d0,
                        self.ra,
                        self.lda,
                        self.solvent_beamwaist,
                        N_COMPONENTS,
                        INTEGRATION_STEPS,
                    )  # solvent_beta = 0
                    self.solvent_calculation = Fitting(
                        self,
                        self.solvent_curves,
                        self.solvent_DPhi0,
                        self.solvent_beamwaist,
                        self.solventCA_zeroLevel,
                        self.solventCA_centerPoint,
                        len(self.solvent_data_set[0]),
                        line_data,
                    )
                    self.result = self.solvent_calculation.manual(
                        self.solventCA_zeroLevel, self.solventCA_centerPoint, self.solvent_DPhi0, self.solvent_beamwaist, self.z_range
                    )
                    self.solvent_curves = Integration(
                        self,
                        0,
                        self.solvent_n2,
                        self.solvent_DPhi0,
                        self.solvent_data_set[0],
                        self.d0,
                        self.ra,
                        self.lda,
                        self.solvent_beamwaist,
                        N_COMPONENTS,
                        INTEGRATION_STEPS,
                    )  # solvent_beta = 0
                    self.solvent_calculation = Fitting(
                        self,
                        self.solvent_curves,
                        self.solvent_DPhi0,
                        self.solvent_beamwaist,
                        self.solventCA_zeroLevel,
                        self.solventCA_centerPoint,
                        len(self.solvent_data_set[0]),
                        line_data,
                    )
                    self.result = self.solvent_calculation.manual(
                        self.solventCA_zeroLevel, self.solventCA_centerPoint, self.solvent_DPhi0, self.solvent_beamwaist, self.z_range
                    )
                    self.solventCA_autofit_done = False
                # elif stype == "OA":
                #     line = self.solventOA_figure.axes.get_lines()[0]
                #     line_data = line.get_ydata()

                #     self.solvent_curves = Integration(self,self.solventOA_T_doubleSpinBox.value(),self.solvent_n2,0,self.solvent_data_set[0],self.d0,self.ra,self.lda,self.solvent_beamwaist,N_COMPONENTS,INTEGRATION_STEPS, stype="OA")
                #     self.solvent_calculation = Fitting(self,self.solvent_curves, 0, self.solvent_beamwaist, self.solventOA_zeroLevel, self.solventOA_centerPoint,len(self.solvent_data_set[0]),line_data)
                #     self.result = self.solvent_calculation.manual(self.solventOA_zeroLevel, self.solventOA_centerPoint, 0, self.solvent_beamwaist, self.z_range, stype)
                #     self.solvent_curves = Integration(self,self.solventOA_T_doubleSpinBox.value(),self.solvent_n2,0,self.solvent_data_set[0],self.d0,self.ra,self.lda,self.solvent_beamwaist,N_COMPONENTS,INTEGRATION_STEPS, stype="OA")
                #     self.solvent_calculation = Fitting(self,self.solvent_curves, 0, self.solvent_beamwaist, self.solventOA_zeroLevel, self.solventOA_centerPoint,len(self.solvent_data_set[0]),line_data)
                #     self.result = self.solvent_calculation.manual(self.solventOA_zeroLevel, self.solventOA_centerPoint, 0, self.solvent_beamwaist, self.z_range, stype)

                #     self.draw_fitting_line(ftype,stype)

                #     self.solventOA_fittingDone = False

            case "sample":
                pass

        self.draw_fitting_line(ftype, stype)
        self.get_curve_interpretation(ftype, stype, "from_geometry")

        caller = "manual"
        self.set_fit_summary(ftype, stype, caller)

    def read_header_params(self, caller: str, ftype: str):
        if caller == "Current Measurement":
            # General parameters
            # No longer do the "set custom" checkboxes affect these actions
            self.silicaThickness_dataFittingTab_doubleSpinBox.setValue(self.silicaThickness_dataSavingTab_doubleSpinBox.value())
            self.wavelength_dataFittingTab_doubleSpinBox.setValue(self.wavelength_dataSavingTab_doubleSpinBox.value())
            self.zscanRange_doubleSpinBox.setValue(np.abs(self.endPos_doubleSpinBox.value() - self.startPos_doubleSpinBox.value()))

            if ftype == "sample":
                self.concentration_dataFittingTab_doubleSpinBox.setValue(self.concentration_dataSavingTab_doubleSpinBox.value())

        elif caller == "Load From File":
            # Get parameters from the file header
            if len(self.header) != 0:
                for hl in self.header:
                    if hl.find("silica thickness") != -1:
                        silica_thickness_match = re.search(r"(?<=silica thickness:)(.*)(?=mm)", hl)
                        silica_thickness = float(silica_thickness_match.groups()[0])
                        self.silicaThickness_dataFittingTab_doubleSpinBox.setValue(silica_thickness)

                    if hl.find("Wavelength") != -1:
                        wavelength_match = re.search(r"(([1-9][0-9]*\.?[0-9]*)|(\.[0-9]+))([Ee][+-]?[0-9]+)?", hl)
                        wavelength = float(wavelength_match.groups()[0])
                        self.wavelength_dataFittingTab_doubleSpinBox.setValue(wavelength)

                    if self.customZscanRange_checkBox.isChecked() is False and hl.find("Starting pos") != -1:  # read zscanRange for any ftype
                        starting_pos_match = re.search(r"(([0-9][0-9]*\.?[0-9]*)|(\.[0-9]+))([Ee][+-]?[0-9]+)?", hl)
                        starting_pos = float(starting_pos_match.groups()[0])
                        next_line = self.header[self.header.index(hl) + 1]
                        end_pos_match = re.search(r"(([0-9][0-9]*\.?[0-9]*)|(\.[0-9]+))([Ee][+-]?[0-9]+)?", next_line)
                        end_pos = float(end_pos_match.groups()[0])
                        self.zscanRange_doubleSpinBox.setValue(np.abs(end_pos - starting_pos))

                    if ftype == "sample" and hl.find("Concentration") != -1:
                        concentration_match = re.search(
                            r"(([0-9][0-9]*\.?[0-9]*\s*%)|(\.[0-9]()\s*%+))([Ee][+-]?[0-9]()\s*%+)?", hl
                        )  # Here make sure there is % symbol
                        self.concentr_percent = float(
                            concentration_match.groups()[0].replace(" ", "")[:-1]
                        )  # remove redundant space and % symbol and change to float
                        self.concentration_dataFittingTab_doubleSpinBox.setValue(self.concentr_percent)

    # QUITTING THE PROGRAM
    def closeEvent(self, event):
        # reply = QMessageBox.question(self, 'Window Close', 'Are you sure you want to close the program?', QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

        # if reply == QMessageBox.Yes:
        #     self.stop_experiment()  # otherwise the thread with running experiment still runs
        #     save_settings(self,self.settings)  # noqa: F405
        #     # print('Program exited.')
        #     event.accept()
        # else:
        #     event.ignore()
        pass


if __name__ == "__main__":
    # ensure that every relative path down the script is referring to current directory of the script
    # wherever it is being run and in whatever environment
    if os.getcwd() != os.path.dirname(__file__):
        os.chdir(os.path.dirname(__file__))

    last_path = load_settings()  # noqa: F405
    if last_path is None:
        print("Try getting write access to the directory or move the program folder to another one, in which you have write access.")
    else:
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "-s", "--settings", dest="settings_ini", default=f"{last_path}", help="Path to settings file (with the filename.*ini)", type=str
        )
        args = parser.parse_args()
        
        app = QApplication(sys.argv)
        window = Window(args.settings_ini)
        
        # Apply temporary widgets
        GUI_Modifiers.cover_widget(window.silicaCA_fittingsummary_3)
        GUI_Modifiers.cover_widget(window.labels_frame_2)
        GUI_Modifiers.cover_widget(window.gridLayout_48)
        GUI_Modifiers.cover_widget(window.frame_3)
        GUI_Modifiers.cover_widget(window.silica_fit_pushButton_2)
        GUI_Modifiers.cover_widget(window.bottom_layout_4)

        window.show()
        
        # Apply styles
        '''The stylesheet is applied to the widgets only after the window is shown
        because the rendering process needs to be completed before the styles can
        be applied. This ensures that all the widgets are properly initialized and
        laid out before any visual changes are made.'''
        app.setStyle("Fusion")
        window.default_palette = QtGui.QGuiApplication.palette()
        window.changeSkinDark()
        window.changeReadOnlyStyle()
        
        app.exec_()
