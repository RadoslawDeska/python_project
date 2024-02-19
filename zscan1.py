"""Z-scan measurement program allows for automated data collection
from simultaneous measurement of open- and closed-aperture Z-scan traces.

Additionally, the program allows for data analysis by fitting the Z-scan traces
to theoretical functions developed by Mansoor Sheik-Bahae and Eric van Stryland
published in \"Characterization Techniques and Tabulations for Organic Nonlinear Materials\"
M. G. Kuzyk and C. W. Dirk, Eds., page 655-692, Marcel Dekker, Inc., 1998
"""

__author__ = "Radosław Deska"
__version__ = '0.1.1'

import json
import logging
import os
import re
import sys
import time
import traceback
import winsound
from datetime import datetime
from math import factorial
from typing import Tuple

import matplotlib
import nidaqmx
import numpy as np
from lmfit import Minimizer, Parameters  #, fit_report
from nidaqmx.constants import AcquisitionType, Edge
from numpy.typing import NDArray
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtCore import QFile, QFileInfo, QObject, QSettings, QThreadPool, QTimer
from PyQt5.QtGui import QColor, QPalette
from PyQt5.QtWidgets import QFileDialog, QMessageBox, QSlider
from scipy.signal import medfilt
from scipy.special import hyp2f1
from sigfig import round as error_rounding


from lib.cursors import BlittedCursor
from lib.figure import MplCanvas
from lib.mgmotor import MG17Motor
# from lib.scientific_rounding import error_rounding
from lib.worker import Worker

matplotlib.rcParams.update({'font.size': 7})
#from matplotlib.widgets import BlittedCursor

# CONSTANTS
SILICA_BETA = 0
N_COMPONENTS = 8 # number of electric field components (for Gaussian decomposition)
INTEGRATION_STEPS = 30 # accuracy of integration infinitesimal element, dx.
CUVETTE_PATH_LENGTH = 0.001 # [m] path length inside cuvette
SOLVENT_T_SLIDER_MAX = 1
MAX_DPHI0 = 3.142 # maximum DeltaPhi0 for silica (for sliders)

# @Adam notes:
# 0. stdlib Pathlib  a= Path("solent.json").read_text()
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


def some_foo(s: str) -> str:
    return s + "ss"

def apply_settings(w, settings):
    # Additions to UI design
    w.mainDirectory_lineEdit.setText(settings.value('SavingTab/main_directory').replace("/","\\"))
    w.dataDirectory_lineEdit.setText(settings.value('FittingTab/data_directory').replace("/","\\"))
class Window(QtWidgets.QMainWindow):

# @dataclass
# class Configuration:
#     path: Path
#     solvent_file: Path = "solvent.json"
    
    
#     @classmethod
#     def load_from_file(config_filename) -> Configuration:
#         with open(cnfig_filename) as f:
#             data = json.load(f)
            
#         self.path = data["path"]
#         self.solvent_flename =\ ///

# INITIALIZATION
    def __init__(self, settings):
        super(Window, self).__init__()
        
        info = QFileInfo(settings)
        if QFile(settings).exists() and bool(info.permissions() & QFile.WriteUser):
            # if file exists and is writable            
            self.settings = QSettings(settings, QSettings.IniFormat)
        else:
            settings_lines = [
                "[UI]",
                fr"ui_path={os.path.join(os.path.dirname(__file__), 'window.ui')}".replace("\\","/"),
                "[MeasurementTab]",
                "starting_position=35",
                "ending_position=75",
                "step_per_cycle=200",
                "samples_per_position=200",
                "[SavingTab]",
                "silica_thickness=4",
                fr"main_directory={os.path.join(os.path.dirname(__file__),'data')}".replace("\\","/"),
                "[FittingTab]",
                fr"data_directory={os.path.join(os.path.dirname(__file__),'data')}".replace("\\","/"),
                "silica_thickness=4",
                "wavelength=800",
                "zrange=40",
                "aperture_diameter=1",
                "distance_from_focus_to_CA=260",
                fr"solvents_path={os.path.join(os.path.dirname(__file__), 'solvents.json')}".replace("\\","/")
                ]
            default_settings_str = '\n'.join(settings_lines)
            
            if not QFile(os.path.join(os.path.dirname(__file__), "default_settings.ini")).exists():
                # create default settings
                with open(os.path.join(os.path.dirname(__file__), "default_settings.ini"), mode="w", encoding="utf-8") as fi:
                    fi.write(default_settings_str)
                os.chmod(os.path.join(os.path.dirname(__file__), "default_settings.ini"), 0o444) # prevent editing by setting permissions to read-only
            
            # create settings file that will be modified
            with open(os.path.join(os.path.dirname(__file__), "settings.ini"), mode="w", encoding="utf-8") as fi:
                fi.write(default_settings_str)
            self.settings = QSettings(os.path.join(os.path.dirname(__file__), "settings.ini"), QSettings.IniFormat)
                
        try:
            uic.loadUi(self.settings.value('UI/ui_path'), self)
        except FileNotFoundError:
            uic.loadUi('window.ui', self)
        apply_settings(self,self.settings)
        
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
        
        # SHOW THE APP WINDOW
        self.show()

    def timing_and_threading(self):
        self.timer=QTimer()
        self.threadpool = QThreadPool()
        print(f"Multithreading with maximum {self.threadpool.maxThreadCount()} threads")

    def states(self):
        self.clearing = False
        self.data_acquisition_complete = False
        self.experiment_stopped = False
        self.initialized = False
        self.initializing = False # variable for watching if initialize method has been called
        self.running = False # Variable for watching if run method has been called
        self.silica_fittingLine_drawn = False
        self.silica_autofit_done = False
        self.solventCA_fittingLine_drawn = False
        self.solventCA_autofit_done = False
        self.solventOA_fittingLine_drawn = False
        self.solventOA_fittingDone = False
        self.sampleCA_fittingLine_drawn = False
        self.sampleCA_fittingDone = False
        self.sampleOA_fittingLine_drawn = False
        self.sampleOA_fittingDone = False
    
    def additional_variables(self):
        self.motor_list = []
        self.offset = 0
        self.where_to_start = "start" # or "end" - whether the scan starts from start_pos or from end_pos
        
        self.previous_end_pos = self.endPos_doubleSpinBox.value()
        self.previous_start_pos = self.startPos_doubleSpinBox.value()
        self.previous_steps_scan = self.stepsScan_spinBox.value()

        self.header_correct = False
        self.rms_value = 0.0 # RMS of Reference signal

    def additional_objects(self):
        # Motor control
        self.ocx = MG17Motor()
        ocx_layout = self.mg17motor_control_vlayout
        ocx_layout.addWidget(self.ocx)
        
        # Charts
        self.initialize_measurement_charts()
        self.initialize_fitting_charts()
        
        # Solvents list
        self.load_solvents()
        # Others
        
    def load_solvents(self, caller=""):
        # Populate Solvent combobox with data from file
        if caller == "":
            try:
                with open(os.path.join(os.path.dirname(__file__), 'solvents.json'), mode="r", encoding="utf-8") as json_file:
                    self.load_and_autocomplete(json_file)
                
            except FileNotFoundError:
                self.showdialog('Warning','solvents.json not found in the default location. Select the file.')
                path = os.path.abspath(os.path.dirname(__file__)) # this is where solvents.json is expected to be
                file = QFileDialog.getOpenFileName(self, "Open File", path,filter="JSON file (*.json)")
                if file[0]!='': # if dialog was not cancelled
                    with open(file[0], mode="r", encoding="utf-8") as json_file:
                        self.load_and_autocomplete(json_file)
            
        elif caller == "LoadSolvents":
            path = os.path.abspath(os.path.dirname(__file__)) # this is where solvents.json is expected to be
            file = QFileDialog.getOpenFileName(self, "Open File", path,filter="JSON file (*.json)")
            if file[0]!='': # if dialog was not cancelled
                with open(file[0], mode="r", encoding="utf-8") as json_file:
                    self.load_and_autocomplete(json_file)
            else:
                return

    def load_and_autocomplete(self, json_file):
        """
        The function loads a JSON file containing solvents, populates a combo box with solvent names,
        and enables autocomplete functionality.
        
        :param json_file: The `json_file` parameter in the `load_and_autocomplete` method is expected to
        be a file object containing solvent data in JSON format. This method reads the JSON data from
        the file, populates the `self.solvents` attribute with the data, adds the solvent names to a
        combo box
        """
        self.solvents = json.load(json_file)
        self.solventName_comboBox.addItems([key for key in self.solvents.keys()])
        self.solvent_autocomplete()
            
    def value_change_triggers(self):
            # Measurement Tab
        self.endPos_doubleSpinBox.editingFinished.connect(lambda: self.measurement_plot_rescale("end"))
        self.startPos_doubleSpinBox.editingFinished.connect(lambda: self.measurement_plot_rescale("start"))
        self.stepsScan_spinBox.valueChanged.connect(self.measurement_plot_rescale)
            
            # Data saving Tab
        self.concentration_dataSavingTab_doubleSpinBox.editingFinished.connect(
            lambda: self.concentration_dataSavingTab_doubleSpinBox.setText(
                self.concentration_dataSavingTab_doubleSpinBox.text()+" %")
            if "%" not in self.concentration_dataSavingTab_doubleSpinBox.text()
            else self.concentration_dataSavingTab_doubleSpinBox.text())
        self.wavelength_dataSavingTab_doubleSpinBox.editingFinished.connect(
            lambda: self.wavelength_dataSavingTab_doubleSpinBox.setText(
                self.wavelength_dataSavingTab_doubleSpinBox.text()+" nm") 
            if "nm" not in self.wavelength_dataSavingTab_doubleSpinBox.text() 
            else self.wavelength_dataSavingTab_doubleSpinBox.text())
            
            # Data fitting Tab
        self.solventName_comboBox.currentIndexChanged.connect(self.solvent_autocomplete)
        self.zscanRange_doubleSpinBox.editingFinished.connect(self.set_new_positions)
    
    def slider_triggers(self):
        # Update display related to the sliders
            # Silica
        self.slider_fit_manually_connect(self.silica_RayleighLength_slider,"Connect")
        self.slider_fit_manually_connect(self.silica_centerPoint_slider, "Connect")
        self.slider_fit_manually_connect(self.silica_zeroLevel_slider, "Connect")
        self.slider_fit_manually_connect(self.silica_DPhi0_slider, "Connect")
        # self.silica_RayleighLength_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Silica", stype="CA"))
        # self.silica_centerPoint_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Silica", stype="CA"))
        # self.silica_zeroLevel_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Silica", stype="CA"))
        # self.silica_DPhi0_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Silica", stype="CA"))
        self.silica_filterSize_slider.valueChanged.connect(lambda: self.reduce_noise_in_data(self.silica_data_set, ftype="Silica", stype="CA"))

            # Solvent
        self.solventCA_RayleighLength_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="CA"))
        self.solventCA_centerPoint_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="CA"))
        self.solventCA_zeroLevel_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="CA"))
        self.solventCA_DPhi0_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="CA"))
        self.solventCA_filterSize_slider.valueChanged.connect(lambda: self.reduce_noise_in_data(self.solvent_data_set, ftype="Solvent", stype="CA"))

        self.solventOA_centerPoint_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="OA"))
        self.solventOA_zeroLevel_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="OA"))
        self.solventOA_T_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="OA"))
        self.solventOA_filterSize_slider.valueChanged.connect(lambda: self.reduce_noise_in_data(self.solvent_data_set, ftype="Solvent", stype="OA"))

    def clicker_triggers(self):
        # Menu triggers
        self.actionExit.triggered.connect(self.shutdown)
        self.actionLoadSolvents.triggered.connect(lambda: self.load_solvents(caller="LoadSolvents"))
        self.actionLight.triggered.connect(self.changeSkinLight)
        self.actionDark.triggered.connect(self.changeSkinDark)

        # Top bar
        self.exitRed_pushButton.clicked.connect(self.stop_experiment)

        # Measurement Tab
        self.clear_pushButton.clicked.connect(self.measurement_clear)
        self.focusAt_comboBox.currentTextChanged.connect(lambda: self.measurement_plot_rescale(self.focusAt_comboBox.currentText()))
        self.initialize_pushButton.clicked.connect(self.initialize)
        self.run_pushButton.clicked.connect(self.set_to_start)

        # Save data Tab
        self.chooseDirectory_pushButton.clicked.connect(lambda: self.choose_dir(caller="DataSaving"))
        self.saveData_pushButton.clicked.connect(self.data_save)
        self.sendToFit_pushButton.clicked.connect(lambda: self.data_loader(caller="Current Measurement", ftype=self.cuvetteType_comboBox.currentText()))
        self.update_pushButton.clicked.connect(self.update_data_and_filenames)

        # Data fitting Tab
            # Files
        self.customDirectory_pushButton.clicked.connect(lambda: self.choose_dir(caller="DataFitting"))
        self.customSilicaFile_pushButton.clicked.connect(lambda: self.data_loader(caller="Load From File", ftype="Silica"))
        self.customSolventFile_pushButton.clicked.connect(lambda: self.data_loader(caller="Load From File", ftype="Solvent"))
        self.customSampleFile_pushButton.clicked.connect(lambda: self.data_loader(caller="Load From File", ftype="Sample"))
            # General parameters
        self.customApertureDiameter_checkBox.stateChanged.connect(lambda: self.enable_custom('ApertureDiameter'))
        self.customApertureToFocusDistance_checkBox.stateChanged.connect(lambda: self.enable_custom('ApertureDistance'))
        self.customSilicaThickness_checkBox.stateChanged.connect(lambda: self.enable_custom('SilicaThickness'))
        self.customWavelength_checkBox.stateChanged.connect(lambda: self.enable_custom('Wavelength'))
        self.customZscanRange_checkBox.stateChanged.connect(lambda: self.enable_custom('ZscanRange'))
            # Sample properties
        self.customConcentration_checkBox.stateChanged.connect(lambda: self.enable_custom('Concentration'))
            
            # Silica CA tab
        self.silica_fit_pushButton.clicked.connect(lambda: self.fit_automatically(ftype="Silica", stype="CA"))
        self.silica_fixROI_checkBox.stateChanged.connect(lambda: self.enable_cursors(ftype="Silica", stype="CA"))

            # Solvent CA tab
        self.solventCA_fit_pushButton.clicked.connect(lambda: self.fit_automatically(ftype="Solvent", stype="CA"))
        self.solventCA_fixROI_checkBox.stateChanged.connect(lambda: self.enable_cursors(ftype="Solvent", stype="CA"))
        self.solventCA_customBeamwaist_checkBox.stateChanged.connect(lambda: self.enable_custom('SolventBeamwaist'))

            # Solvent OA tab
        self.solventOA_fit_pushButton.clicked.connect(lambda: self.fit_automatically(ftype="Solvent", stype="OA"))
        self.solventOA_fixROI_checkBox.stateChanged.connect(lambda: self.enable_cursors(ftype="Solvent", stype="OA"))
        self.solventOA_isAbsorption_checkBox.stateChanged.connect(lambda: self.toggle_absorption_model(ftype="Solvent"))
        self.solventOA_absorptionModel_comboBox.currentIndexChanged.connect(lambda: self.toggle_saturation_model(ftype="Solvent"))
        self.solventOA_customCenterPoint_checkBox.stateChanged.connect(lambda: self.enable_custom('SolventCenterPoint'))
        
    def timer_triggers(self): #(started with Initalize button click)
        self.timer.timeout.connect(self.motion_detection)
        
        self.start_timer()
    
    def initialize(self, *args, **kwargs):
        self.initializing = True

        # Initialize detectors
        device_name = "/Dev1"
        self.detector_core_name = device_name+"/ai"
        self.number_of_channels_used = 3
        
        print('Detectors initialized')

        # Initialize data dictionaries and apply empty data to lines
        for type, chart in self.charts.items():             # e.g.: take the tuple ("relative", "self.rel_chart")
            for chan_no in range(self.number_of_channels_used):
                self.data[type].update({chan_no: []})       # then fill "relative" dictionary in "self.data" dictionary with pairs of channel number and list of values.
                line, = chart.axes.plot(self.data["positions"],self.data[type][chan_no], marker='.') # add empty line for each channel on the "relative" chart
                self.measurement_lines[type].update({chan_no: line})    # update the "lines" dictionary with line for each channel
        
        # Initialize translation stage motor
        self.motor_detection_and_homing()

        self.initialized = True
        self.initializing = False

# INITIALIZE CHARTS
    def initialize_measurement_charts(self):
        '''Initializes charts in the 'Measurement' Tab
        and sets their scales using 'measurement_plot_rescale()' method.'''
        # "Measurement charts"
        self.rel_chart = MplCanvas(self)
        layout_rel = self.relative_layout
        layout_rel.addWidget(self.rel_chart)

        self.abs_chart = MplCanvas(self)
        self.abs_chart.axes.set_ylabel('Amplitude (V)')
        layout_abs = self.absolute_layout
        layout_abs.addWidget(self.abs_chart)

        self.rms_text = self.abs_chart.axes.text(0.87,0.9, f"RMS noise = {self.rms_value*100:.3f}%", transform=self.abs_chart.axes.transAxes,
                                                 bbox = dict(boxstyle='round', facecolor='white', alpha=1))
        
        self.charts = {"relative": self.rel_chart, "absolute": self.abs_chart}
        # initialize empty lines and data dictionaries
        self.measurement_lines = {"relative": {}, "absolute": {}}
        self.data = {"positions": [], "relative": {}, "absolute": {}}

        self.measurement_plot_rescale()

        print('Canvas loaded')

    def initialize_fitting_charts(self):
        # Silica chart
        self.silica_figure = MplCanvas(self)
        self.silica_figure.axes.plot([],[], marker='o', ms=5, linestyle='', color='tab:blue')
        self.silica_figure.axes.set_title("Silica - closed aperture")
        self.silica_figure.axes.set_ylabel("Norm. trasmittance")
        self.silica_figure.axes.set_position([0.135,0.105,.825,.825])
        layout_ca_silica = self.silicaCA_layout
        layout_ca_silica.addWidget(self.silica_figure)
        
        self.silicaOA_figure = MplCanvas(self)
        self.silicaOA_figure.axes.plot([],[], marker='o', ms=5, linestyle='', color='tab:blue')
        self.silicaOA_figure.axes.set_title("Silica - open aperture")
        self.silicaOA_figure.axes.set_ylabel("Norm. trasmittance")
        self.silicaOA_figure.axes.set_position([0.135,0.105,.825,.825])
        layout_oa_silica = self.silicaOA_layout
        layout_oa_silica.addWidget(self.silicaOA_figure)

        # Solvent charts
        self.solventCA_figure = MplCanvas(self)
        self.solventCA_figure.axes.plot([],[], marker='o', ms=5, linestyle='', color='tab:blue')
        self.solventCA_figure.axes.set_title("Solvent - closed aperture")
        self.solventCA_figure.axes.set_ylabel("Norm. trasmittance")
        self.solventCA_figure.axes.set_position([0.135,0.105,.825,.825])
        layout_ca_solvent = self.solventCA_layout
        layout_ca_solvent.addWidget(self.solventCA_figure)
        
        self.solventOA_figure = MplCanvas(self)
        self.solventOA_figure.axes.plot([],[], marker='o', ms=5, linestyle='', color='tab:blue')
        self.solventOA_figure.axes.set_title("Solvent - open aperture")
        self.solventOA_figure.axes.set_ylabel("Norm. trasmittance")
        self.solventOA_figure.axes.set_position([0.135,0.105,.825,.825])
        layout_oa_solvent = self.solventOA_layout
        layout_oa_solvent.addWidget(self.solventOA_figure)

        # Sample charts
        self.sampleCA_figure = MplCanvas(self)
        self.sampleCA_figure.axes.plot([],[], marker='o', ms=5, linestyle='', color='tab:blue')
        self.sampleCA_figure.axes.set_title("Sample - closed aperture")
        self.sampleCA_figure.axes.set_ylabel("Norm. trasmittance")
        self.sampleCA_figure.axes.set_position([0.135,0.105,.825,.825]) # [left, bottom, width, height]
        layout_ca_sample = self.sampleCA_layout
        layout_ca_sample.addWidget(self.sampleCA_figure)
        
        self.sampleOA_figure = MplCanvas(self)
        self.sampleOA_figure.axes.plot([],[], marker='o', ms=5, linestyle='', color='tab:blue')
        self.sampleOA_figure.axes.set_title("Sample - open aperture")
        self.sampleOA_figure.axes.set_ylabel("Norm. trasmittance")
        self.sampleOA_figure.axes.set_position([0.135,0.105,.825,.825])
        layout_oa_sample = self.sampleOA_layout
        layout_oa_sample.addWidget(self.sampleOA_figure)

        self.fitting_charts = {"Silica": {"CA": self.silica_figure, "OA": self.silicaOA_figure},
                               "Solvent": {"CA": self.solventCA_figure, "OA": self.solventOA_figure},
                               "Sample": {"CA": self.sampleCA_figure, "OA": self.sampleOA_figure}}
        
# MOTOR NAVIGATION
    def motion_detection(self):
        if self.initializing is True:
            self.clearLED_pushButton.setEnabled(False)
            self.initLED_pushButton.setEnabled(True)
            self.initialize_pushButton.setEnabled(False)
            self.runLED_pushButton.setEnabled(False)
        elif self.running is True:
            self.clearLED_pushButton.setEnabled(False)
            self.initLED_pushButton.setEnabled(False)
            self.runLED_pushButton.setEnabled(True)
        elif self.clearing is True:
            self.clearLED_pushButton.setEnabled(True)
            self.initLED_pushButton.setEnabled(False)
            self.runLED_pushButton.setEnabled(False)
        else:
            self.clearLED_pushButton.setEnabled(False)
            self.initLED_pushButton.setEnabled(False)
            self.runLED_pushButton.setEnabled(False)

        if hasattr(self,'motor'):
            if self.motor.is_in_motion is True:
                self.clear_pushButton.setEnabled(False)
                self.run_pushButton.setEnabled(False)
                self.waitLED_pushButton.setEnabled(True)
                
                self.stepsScan_spinBox.setEnabled(False)
                self.samplesStep_spinBox.setEnabled(False)

            else:
                if self.initializing is True:
                    self.clear_pushButton.setEnabled(False)
                    self.run_pushButton.setEnabled(False)
                    self.waitLED_pushButton.setEnabled(True)
                    
                    self.startPos_doubleSpinBox.setEnabled(False)
                    self.endPos_doubleSpinBox.setEnabled(False)
                    self.stepsScan_spinBox.setEnabled(False)
                    self.samplesStep_spinBox.setEnabled(False)

                elif self.running is True:
                    self.clear_pushButton.setEnabled(False)
                    self.run_pushButton.setEnabled(False)
                    self.waitLED_pushButton.setEnabled(True)

                    self.startPos_doubleSpinBox.setEnabled(False)
                    self.endPos_doubleSpinBox.setEnabled(False)
                    self.stepsScan_spinBox.setEnabled(False)
                    self.samplesStep_spinBox.setEnabled(False)

                elif self.clearing is True:
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
        
        if self.data_acquisition_complete is False:
            self.update_pushButton.setEnabled(False)
            self.sendToFit_pushButton.setEnabled(False)
        else:
            self.update_pushButton.setEnabled(True)

    def motor_detection_and_homing(self, *args, **kwargs):
        motor_id = 40180184

        from packages import thorlabs_apt as apt
        # importing it from custom location allows to place APT.dll
        # in the package directory to read it (it is more user-friendly)
        
        try:
            self.motor = apt.Motor(motor_id)
            self.ocx.configure(motor_id)
            print(f'Motor {motor_id} connected')
        
        except Exception:
            logging.error(traceback.format_exc())
            print('Try closing other programs that may be accessing the device.')
            print('Try re-plugging the USB device and then run this program.')
            print(f'If you are simulating the motion controller, make sure Thorlabs APT Server has been configured with the HWSerialNumber {motor_id}')
            self.showdialog('Error', 'Motion controller not found! For more information, see the console output.')
        
        self.mpositioner = MotorPositioner()
        self.motor.backlash_distance = 0
        self.thread_it(self.mpositioner.movehome)

    def set_custom_pos(self):
        self.thread_it(self.mpositioner.movetocustompos)

    def set_to_start(self):
        self.running = True
        self.thread_it(self.mpositioner.movetostart)

# GUI EVENTS TIMING
    def start_timer(self):
        self.timer.start(100) # How often various actions or states are triggered or checked for
    
    #def stop_timer(self):
    #    self.timer.stop()

    def stop_experiment(self):
        if hasattr(self,'motor'):
            self.motor.stop_profiled()
            self.experiment_stopped = True
                        
            self.initializing = False
            self.running = False
            self.clearing = False
            
            self.measurement_clear()

# QUITTING THE PROGRAM
    def shutdown(self):
        """
        """        
        reply = QMessageBox.question(self, 'Window Close', 'Are you sure you want to close the window?', QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            # new_settings = dict()
            # self.mainDirectory_lineEdit.setText(self.path.replace("/","\\"))
            # self.dataDirectory_lineEdit.setText(self.path.replace("/","\\"))
            # self.settings.setValue()
            print('Program exited.')
            sys.exit()

# DATA ACQUISITION AND DISPLAY
    def measurement_clear(self): # clears all in the first two tabs (Measurement and Data Saving)
        self.clearing = True
        self.data_acquisition_complete = False

        self.data["positions"] = []     # reset positions to empty list
        
        self.rms_value = 0.0
        self.rms_text.set_text(f"RMS noise = {self.rms_value*100:.3f}%")

        for type, chart in self.charts.items():     # e.g.: take the tuple ("relative", "self.rel_chart")
            for chan_no in range(self.number_of_channels_used):
                self.data[type][chan_no] = []       # then fill "relative" dictionary in "self.data" dictionary empty list of values per channel
                
                # UPDATE LINES INSTEAD OF DELETING AND REINSTANTIATING
                self.measurement_lines[type][chan_no].set_xdata(self.data["positions"]) # and set data of lines in "lines" dictionary to empty lists of positions and values
                self.measurement_lines[type][chan_no].set_ydata(self.data[type][chan_no])
            
            chart.axes.relim()
            chart.axes.autoscale_view()
            chart.draw_idle()

        self.rawLogData_textBrowser.clear()
        self.fullLogData_textBrowser.clear()

        self.saveData_pushButton.setEnabled(False)
        
        self.clearing = False

    def measurement_plot_rescale(self, who_called=""):
        '''Rescales plots in the 'Measurement' Tab based on
        values given in the 'Measurement control' panel\n
        If user changes the focus of chart at specific line, rescaling fits the charts
        to display the line in its min-max y-range.'''
        
        match who_called:
            case "start":
                if self.startPos_doubleSpinBox.value() >= self.endPos_doubleSpinBox.value():
                    self.startPos_doubleSpinBox.setValue(self.previous_start_pos)
            case "end":
                if self.endPos_doubleSpinBox.value() <= self.startPos_doubleSpinBox.value():
                    self.endPos_doubleSpinBox.setValue(self.previous_end_pos)
            
            case _: # if the who_called value is any of the values in the list of "Focus At" dropdown OR "ANYTHING ELSE"
                
                if self.initialized is False: # prevent the problem of accessing data for rescale when there is no data created yet.
                    return
                
                called_by = self.focusAt_comboBox.currentText()
                for type, chart in self.charts.items():
                    match called_by:
                        case "All":
                            chart.axes.relim()
                            chart.axes.autoscale()
                        case "Closed":
                            try:
                                chart.axes.set_ylim(min(self.data[type][0]),max(self.data[type][0]))
                                chart.axes.relim()
                            except ValueError:
                                #print(f'\nValueError: min() arg is an empty sequence.\nThis message was initiated by user clicking on {called_by}.\nNothing has to be done.')
                                return
                        case "Reference":
                            try:
                                chart.axes.set_ylim(min(self.data[type][1]),max(self.data[type][1]))
                                chart.axes.relim()
                            except ValueError:
                                #print(f'\nValueError: min() arg is an empty sequence.\nThis message was initiated by user clicking on {called_by}.\nNothing has to be done.')
                                return
                        case "Open":
                            try:
                                chart.axes.set_ylim(min(self.data[type][2]),max(self.data[type][2]))
                                chart.axes.relim()
                            except ValueError:
                                #print(f'\nValueError: min() arg is an empty sequence.\nThis message was initiated by user clicking on {called_by}.\nNothing has to be done.')
                                return
                        case _: # HERE GOES "ANYTHING ELSE"
                            return

                    chart.draw_idle()
        
        self.focalPoint_doubleSpinBox.setValue((self.endPos_doubleSpinBox.value()+self.startPos_doubleSpinBox.value())/2)
        self.zscanRange_measurementTab_doubleSpinBox.setValue(np.abs(self.endPos_doubleSpinBox.value()-self.startPos_doubleSpinBox.value()))

        self.offset = (self.endPos_doubleSpinBox.value()-self.startPos_doubleSpinBox.value())/self.stepsScan_spinBox.value()
        for chart in self.charts.values():
            chart.axes.set_xlim(self.startPos_doubleSpinBox.value()-self.offset, self.endPos_doubleSpinBox.value()+self.offset)
            
            chart.axes.relim()
            chart.draw_idle()

    def create_raw_log_line(self, step):
        if window.where_to_start == "end":
            new_step = np.abs(step-200)
            line = f"{new_step:4d}{' '}"
        else:
            line = f"{step:4d}{' '}"

        for chan_no in range(self.number_of_channels_used):
            line += f"{2*' '}{self.data['absolute'][chan_no][step]:12.8f}{9*' '}"
            #line += f"{2*' '}{abs_chart_data[chan_no][step]:12.8f}{9*' '}"

            if chan_no == self.number_of_channels_used-1:
                line += f"{2*' '}{0:12.8f}"

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
        
        p = QFileDialog.getExistingDirectory(self, 'Select a directory', path)
        if p != "": # This keeps old path in directory QLineEdit lines, if dialog is closed with Cancel
            if caller == "DataSaving":
                self.mainDirectory_lineEdit.setText(p.replace("/","\\"))
            elif caller == "DataFitting":
                self.dataDirectory_lineEdit.setText(p.replace("/","\\"))

# DATA SAVING
    def data_reverse(self):
        if self.where_to_start == "end":# and self.data_reversed is False:
            for chan_no in range(self.number_of_channels_used):
                self.data['absolute'][chan_no] = np.flip(self.data['absolute'][0],axis=0)
            #self.data_reversed = True
    
    def update_data_and_filenames(self):
        self.fullLogData_textBrowser.clear()

    # DATA PREVIEW
        # Full description header
        #if self.experimentDescription_plainTextEdit.toPlainText() != "": # this is for full log to look nicer
        #    self.experimentDescription_plainTextEdit.appendPlainText("")
        
        # full log header
        header = ("Z-scan Measurement\n"                                                                             # line 0
                  #f"Sample type: {self.cuvetteType_comboBox.currentText()}\n"                                           
                  f"Code: {self.codeOfSample_lineEdit.text()}\n"                                                     # line 1
                  f"Silica thickness: {self.silicaThickness_dataSavingTab_doubleSpinBox.text()}\n"                   # line 2
                  f"Concentration: {self.concentration_dataSavingTab_doubleSpinBox.text()}\n"                        # line 3
                  f"Wavelength: {self.wavelength_dataSavingTab_doubleSpinBox.text()}\n"                              # line 4
                  f"{self.experimentDescription_plainTextEdit.toPlainText()}\n"                                      # line 5
                  "--------------------------------------------------------------------------------------------\n\n" # line 6
                  f"Starting pos: {self.startPos_doubleSpinBox.value()}\n"                                           # line 7
                  f"Ending pos: {self.endPos_doubleSpinBox.value()}\n"                                               # line 8
                  "CH1:   Closed aperture\n"                                                                         # line 9
                  "CH2:   Reference\n"                                                                               # line 10
                  "CH3:   Open aperture\n"                                                                           # line 11
                  "CH4:   Empty channel\n\n"                                                                         # line 12
                  "--------------------------------------------------------------------------------------------\n\n" # line 13
                  "SNo.  [V] Voltage Max        [V] Voltage Max        [V] Voltage Max        [V] Voltage Max\n\n"   # line 14
                  "--------------------------------------------------------------------------------------------\n")  # line 15

        raw_log = self.rawLogData_textBrowser.toPlainText()
        self.fullLogData_textBrowser.append(header)
        self.fullLogData_textBrowser.append(raw_log)

        # log data
        if self.data_acquisition_complete is True:
            raw_log_data = np.genfromtxt(raw_log.split('\n'))
            self.data['absolute'] = {0: list(raw_log_data[:,1]), 1: list(raw_log_data[:,2]), 2: list(raw_log_data[:,3])}
            self.data_set = list(raw_log_data[:,0]), self.data["absolute"][0], self.data["absolute"][1], self.data["absolute"][2]#, data[:,3] not using the last column with zeros
            
            self.saveData_pushButton.setEnabled(True)
        else:
            self.saveData_pushButton.setEnabled(False)

    # FILENAMES
        now = datetime.now()
        self.cur_date = now.strftime("%Y_%m_%d")
        self.cur_time = now.strftime("%H_%M")
        sample_type = self.codeOfSample_lineEdit.text()
        
        concentration = self.concentration_dataSavingTab_doubleSpinBox.value()
        conc_hyphen = str(f"{concentration:.2f}").replace(".","-")

        wavelength = self.wavelength_dataSavingTab_doubleSpinBox.value()
        wavel_hyphen = str(f"{wavelength:.1f}").replace(".","-")
        
        self.rawLogFilename_lineEdit.setText(f"{self.cur_date}__{self.cur_time}__{sample_type}_{conc_hyphen}_{wavel_hyphen}.txt")
        self.fullLogFilename_lineEdit.setText(f"{self.cur_date}__{self.cur_time}__{sample_type}_{conc_hyphen}_{wavel_hyphen}_2.txt")

        self.files = (self.rawLogFilename_lineEdit.text(), self.fullLogFilename_lineEdit.text())
    
    def data_save(self):
        self.accurate_path = os.path.join(self.mainDirectory_lineEdit.text(),self.cur_date)
        
        try:
            if not os.path.exists(self.accurate_path):
                os.mkdir(self.accurate_path)
            
            for file in self.files:
                with open(os.path.join(self.accurate_path,file), 'w') as f:
                    if self.files.index(file) == 0:
                        f.write(self.rawLogData_textBrowser.toPlainText())
                    else:
                        f.write(self.fullLogData_textBrowser.toPlainText())
            
            self.saveData_pushButton.setEnabled(False)
            self.sendToFit_pushButton.setEnabled(True)
            self.showdialog('Info', "Files successfully written!")
        
        except PermissionError:
            self.showdialog('Error', 'Permission denied!\nCannot write the file in this directory.')

# DATA FITTING
    def data_loader(self, caller:str, ftype:str):
        '''Load data either from file or from current measurement ('caller' argument).\n
        'ftype' is passed by proper button for loading from file, or by proper option in "Cuvette type" combobox in "Data Saving" tab'''
        if caller == "Current Measurement":
            # Fill in names in QLineEdits
            self.dataDirectory_lineEdit.setText(self.accurate_path+"\\")
            self.mainTabs.setCurrentIndex(2)

            match ftype:
                case "Silica":
                    self.silicaFilename_lineEdit.setText(self.rawLogFilename_lineEdit.text())
                    self.fittingTabs.setCurrentIndex(0)
                    self.silicaAperture_tabWidget.setCurrentIndex(0)
                case "Solvent":
                    self.solventFilename_lineEdit.setText(self.rawLogFilename_lineEdit.text())
                    self.fittingTabs.setCurrentIndex(1)
                    self.solventAperture_tabWidget.setCurrentIndex(0)
                case "Sample":
                    self.sampleFilename_lineEdit.setText(self.rawLogFilename_lineEdit.text())
                    self.fittingTabs.setCurrentIndex(2)
                    self.sampleAperture_tabWidget.setCurrentIndex(0)
            
            # Get data
            self.data_set = range(self.stepsScan_spinBox.value()+1), self.data["absolute"][0], self.data["absolute"][1], self.data["absolute"][2]#, data[:,3] not using the last column with zeros
            
            # Read parameters
            self.read_header_params(caller = "Current Measurement", ftype=ftype)
            self.data_display(self.data_set, ftype)

            self.switch_fitting_to_on_state(ftype)
            if ftype == 'Silica':
                self.silica_autofit_done = False # this is to start afresh with fitting

            self.fit_manually(ftype,stype="CA")
            self.fit_manually(ftype,stype="OA")
        
        elif caller == "Load From File":
            # Load data
            def load_data():
                """Fills proper frame in GUI with file information, sets current tab to lead the user,\n                
                loads data on screen, toggles to active state the fitting controls and calls for the initial fit\n
                for both CA and OA traces.
                """
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
                data = np.genfromtxt(p, skip_header=last_header_line)
                self.data_set = data[:,0], data[:,1], data[:,2], data[:,3]#, data[:,4] not using the last column with zeros
                # Read parameters
                self.read_header_params(caller = "Load From File", ftype=ftype)
                self.data_display(self.data_set, ftype)

                self.switch_fitting_to_on_state(ftype)
                if ftype == 'Silica':
                    self.silica_autofit_done = False # this is to start afresh with fitting

                self.fit_manually(ftype, stype="CA")
                self.fit_manually(ftype, stype="OA")
            
            p = QFileDialog.getOpenFileName(self, 'Select full description file', os.path.normpath(self.dataDirectory_lineEdit.text()))
            if p[0] != "": # This keeps old filename in given file type QLineEdit lines, if dialog is closed with Cancel
                p = p[0]
                fname = p.split("/")[-1]
                self.dataDirectory_lineEdit.setText("\\".join(p.split("/")[:-1]))
                
                # Header check and manipulation
                with open(p, 'r') as file:
                    # Correct header must include these and a beacon at the end in the form of "SNo" substring
                    required_matches = ["Concentration","Wavelength","Starting pos","Ending pos", "SNo"]
                    optional_matches = ["Silica thickness"]
                    header_matches = len(required_matches)*[False]
                    self.header = []
                    last_header_line = 0

                    for line_no, line in enumerate(file):
                        for match in required_matches:
                            if re.match(r"\b%s\b" % match, line): # lookup whole words
                                header_matches[required_matches.index(match)] = True
                                self.header.append(line.strip())

                                if match == "SNo": # This is the header end beacon
                                    last_header_line = line_no+3
                        
                        opt_matched = 0
                        for opt_match in optional_matches:
                            if re.match(r"\b%s\b" % opt_match, line): # lookup whole words
                                header_matches.append(True)
                                self.header.append(line.strip())
                                opt_matched += 1
                    
                    if len(self.header) != 0:
                        if last_header_line != 0: # if the header end beacon was found
                            self.header_correct = all(header_matches[:-(1+opt_matched)]) # check if required matches are satisfied
                            
                            if self.header_correct is False:
                                self.showdialog('Warning',
                                ('Header corrupted!\n\n'
                                'Some parameters might have been loaded improperly.\n\nPut measurement parameters manually.'))

                                self.customSilicaThickness_checkBox.setChecked(True)
                                self.customWavelength_checkBox.setChecked(True)
                                self.customZscanRange_checkBox.setChecked(True)

                            load_data()

                        else:
                            self.header_correct = False
                            self.showdialog('Warning',
                                ('Header corrupted!\n\n'
                                'Load a file without the header.'))
                            
                    else:
                        self.header_correct = False
                        self.showdialog('Warning',
                        ('Header not found!\n\n'
                        'Set measurement parameters by yourself or load a file with the header.'))

                        self.customSilicaThickness_checkBox.setChecked(True)
                        self.customWavelength_checkBox.setChecked(True)
                        self.customZscanRange_checkBox.setChecked(True)

                        load_data()

    def enable_custom(self, o:str):
        """Toggles readOnly parameter on `o` element from GUI. Uses match-case structure with `o` parameter to match\n
        to speed up processing at each call.

        Args:
            o (str): case for parameter to toggle
        """        
        match o:
            case 'ApertureDiameter':
                if self.customApertureDiameter_checkBox.isChecked() is True:
                    self.apertureDiameter_doubleSpinBox.setReadOnly(False)
                else:
                    self.apertureDiameter_doubleSpinBox.setReadOnly(True)
            case 'ApertureDistance':
                if self.customApertureToFocusDistance_checkBox.isChecked() is True:
                    self.apertureToFocusDistance_doubleSpinBox.setReadOnly(False)
                else:
                    self.apertureToFocusDistance_doubleSpinBox.setReadOnly(True)
            case 'SilicaThickness': 
                if self.customSilicaThickness_checkBox.isChecked() is True:
                    self.silicaThickness_dataFittingTab_doubleSpinBox.setReadOnly(False)
                else:
                    self.silicaThickness_dataFittingTab_doubleSpinBox.setReadOnly(True)
            case 'Wavelength':
                if self.customWavelength_checkBox.isChecked() is True:
                    self.wavelength_dataFittingTab_doubleSpinBox.setReadOnly(False)
                else:
                    self.wavelength_dataFittingTab_doubleSpinBox.setReadOnly(True)
            case 'ZscanRange':
                if self.customZscanRange_checkBox.isChecked() is True:
                    self.zscanRange_doubleSpinBox.setReadOnly(False)
                else:
                    self.zscanRange_doubleSpinBox.setReadOnly(True)
            case 'Concentration':
                if self.customConcentration_checkBox.isChecked() is True:
                    self.concentration_dataFittingTab_doubleSpinBox.setReadOnly(False)
                else:
                    self.concentration_dataFittingTab_doubleSpinBox.setReadOnly(True)
            case 'SolventBeamwaist':
                if self.solventCA_customBeamwaist_checkBox.isChecked() is False:
                    self.solventCA_RayleighLength_slider.setEnabled(False)
                    if hasattr(window, 'silica_beamwaist'):
                        self.solventCA_RayleighLength_slider.valueChanged.disconnect()
                        self.solventCA_RayleighLength_slider.setValue(int(round(self.silica_beamwaist*1E6)))
                        self.solventCA_RayleighLength_slider.setValue(int(round(self.silica_beamwaist*1E6)))
                        self.solventCA_RayleighLength_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="CA"))
                        self.solventCA_beamwaistSummary_doubleSpinBox.setValue(self.silica_beamwaist*1E6)
                        self.solventCA_beamwaistSummary_doubleSpinBox.setValue(self.silica_beamwaist*1E6)
                        
                else:
                    self.solventCA_RayleighLength_slider.setEnabled(True)
            case 'SolventCenterPoint': # OA case
                if self.solventOA_customCenterPoint_checkBox.isChecked() is False:
                    self.solventOA_centerPoint_slider.setEnabled(False)
                    #self.solventOA_centerPoint_doubleSpinBox.setValue(self.solventCA_centerPoint_doubleSpinBox.value())
                else:
                    self.solventOA_centerPoint_slider.setEnabled(True)

    def slider_fit_manually_connect(self, current_slider:QSlider, mode):
        '''This function is intended to connect and disconnect other sliders on demand, to prevent them from updating the fitting line:
        1) "Disconnect" means current slider should be kept and all others should disconnect from their method
        2) "Connect" means reconnect all sliders to their method'''
        
        silica_sliders = [self.silica_RayleighLength_slider, self.silica_centerPoint_slider, self.silica_zeroLevel_slider, self.silica_DPhi0_slider]
        #solventCA_sliders = [self.solventCA_RayleighLength_slider, self.solventCA_centerPoint_slider, self.solventCA_zeroLevel_slider, self.solventCA_DPhi0_slider]
        #solventOA_sliders = [self.solventOA_deltaZpv_slider, self.solventOA_centerPoint_slider, self.solventOA_zeroLevel_slider, self.solventOA_deltaTpv_slider]
        #sampleCA_sliders = [self.sampleCA_deltaZpv_slider, self.sampleCA_centerPoint_slider, self.sampleCA_zeroLevel_slider, self.sampleCA_deltaTpv_slider]
        #sampleOA_sliders = [self.sampleOA_deltaZpv_slider, self.sampleOA_centerPoint_slider, self.sampleOA_zeroLevel_slider, self.sampleOA_deltaTpv_slider]
        
        available_sliders = [silica_sliders]#, solventCA_sliders, solventOA_sliders, sampleCA_sliders, sampleOA_sliders]
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
                        pass # If not found
            
            case "Connect":
                match self.disconnected_sliders:
                    case "silica":
                        for s in available_sliders[0]:
                            s.valueChanged.connect(lambda: self.fit_manually(ftype="Silica", stype="CA", activated_by=self.disconnected_sliders))
                    case "silicaOA":
                        for s in available_sliders[1]:
                            s.valueChanged.connect(lambda: self.fit_manually(ftype="Silica", stype="OA", activated_by=self.disconnected_sliders))
                    case "solventCA":
                        for s in available_sliders[2]:
                            s.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="CA", activated_by=self.disconnected_sliders))
                    case "solventOA":
                        for s in available_sliders[3]:
                            s.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="OA", activated_by=self.disconnected_sliders))
                    case "sampleCA":
                        for s in available_sliders[4]:
                            s.valueChanged.connect(lambda: self.fit_manually(ftype="Sample", stype="CA", activated_by=self.disconnected_sliders))
                    case "sampleOA":
                        for s in available_sliders[5]:
                            s.valueChanged.connect(lambda: self.fit_manually(ftype="Sample", stype="OA", activated_by=self.disconnected_sliders))
                    case "All":
                        for s in available_sliders[0]:
                            s.valueChanged.connect(lambda: self.fit_manually(ftype="Silica", stype="CA"))
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
                if self.solventOA_isAbsorption_checkBox.isChecked() is False:
                    self.solventOA_absorptionModel_label.setVisible(True)
                    self.solventOA_absorptionModel_comboBox.setVisible(True)
                    self.solventOA_fixROI_checkBox.setVisible(True)
                else:
                    self.solventOA_absorptionModel_label.setVisible(False)
                    self.solventOA_absorptionModel_comboBox.setVisible(False)
                    self.solventOA_fixROI_checkBox.setVisible(False)
        
            case "Sample":
                if self.sampleOA_isAbsorption_checkBox.isChecked() is False:
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
                if self.solventOA_isAbsorption_checkBox.isChecked() is False:
                    models = ["SA", "2PA+SA"] #, "RSA"]
                    if self.solventOA_absorptionModel_comboBox.currentText() in models:
                        self.solventOA_saturationModel_label.setVisible(True)
                        self.solventOA_saturationModel_comboBox.setVisible(True)
                    else:
                        self.solventOA_saturationModel_label.setVisible(False)
                        self.solventOA_saturationModel_comboBox.setVisible(False)
            case "Sample":
                if self.sampleOA_isAbsorption_checkBox.isChecked() is False:
                    models = ["SA", "2PA+SA"] #, "RSA"]
                    if self.sampleOA_absorptionModel_comboBox.currentText() in models:
                        self.sampleOA_saturationModel_label.setVisible(True)
                        self.sampleOA_saturationModel_comboBox.setVisible(True)
                    else:
                        self.sampleOA_saturationModel_label.setVisible(False)
                        self.sampleOA_saturationModel_comboBox.setVisible(False)

    def switch_fitting_to_on_state(self, ftype:str):
        match ftype:
            case "Silica":
                self.silica_RayleighLength_slider.setEnabled(True)
                self.silica_centerPoint_slider.setEnabled(True)
                self.silica_zeroLevel_slider.setEnabled(True)
                self.silica_DPhi0_slider.setEnabled(True)
                self.silica_fit_pushButton.setEnabled(True)
                self.silica_filterSize_slider.setEnabled(True)

            case "Solvent":
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
        
            case "Sample":
                pass
    
    def set_fit_summary(self, ftype:str, stype:str, caller=""):
        """Sets proper number of digits in summary display and shows parameters (with errors) after (automatic) fitting.

        Args:
            ftype (str): Sample type ("Silica","Solvent","Sample")
            stype (str): Z-scan mode ("CA", "OA")
            caller (str, optional): Which fit type has called this function. Defaults to "".
        """

        # DISPLAY VALUES (ALREADY ROUNDED)
        match ftype:
            case "Silica":
                if caller == "auto":
                    self.silica_deltaPhi0Summary_label.setText(str(self.silica_DPhi0)+' ± '+str(self.silica_DPhi0Error)) # rad
                    self.silica_laserIntensitySummary_label.setText(str(self.laserI0*1E-13)+' ± '+str(self.laserI0Error*1E-13)) # [GW/cm2]
                    self.silica_beamwaistSummary_label.setText(str(self.silica_beamwaist*1E6)+' ± '+str(self.silica_beamwaistError*1E6)) # [um] radius in focal point
                    self.silica_rayleighRangeSummary_label.setText(str(self.silica_rayleighLength*1E3)+' ± '+str(self.silica_rayleighLengthError*1E3)) # [mm]
                    self.numericalApertureSummary_label.setText(str(self.numericalAperture)+' ± '+str(self.numericalApertureError))
                else:
                    self.silica_deltaPhi0Summary_label.setText(str(self.silica_DPhi0)+' ± #.##') # rad
                    self.silica_laserIntensitySummary_label.setText(str(self.laserI0*1E-13)+' ± #.##') # [GW/cm2]
                    self.silica_beamwaistSummary_label.setText(str(error_rounding(self.silica_beamwaist*1E6, sigfigs=3, warn=False))+' ± #.##') # [um] radius in focal point
                    self.silica_rayleighRangeSummary_label.setText(str(self.silica_rayleighLength)+' ± #.##') # [mm]
                    self.numericalApertureSummary_label.setText(str(self.numericalAperture)+' ± #.##')
                
            case "Solvent":
                if stype == "CA":
                    self.solventCA_deltaPhi0Summary_doubleSpinBox.setValue(self.solvent_DPhi0)
                    self.solventCA_n2Summary_doubleSpinBox.setValue(self.solvent_n2*1E22) # display n2 in multiples of 10^-9 cm^2/GW
                    self.solventCA_beamwaistSummary_doubleSpinBox.setValue(self.solvent_beamwaist*1E6) # [um] radius in focal point
                    self.solventCA_deltaPhi0Summary_doubleSpinBox.setValue(self.solvent_DPhi0)
                    self.solventCA_n2Summary_doubleSpinBox.setValue(self.solvent_n2*1E22) # display n2 in multiples of 10^-9 cm^2/GW
                    self.solventCA_beamwaistSummary_doubleSpinBox.setValue(self.solvent_beamwaist*1E6) # [um] radius in focal point
                    self.solventCA_rayleighRangeSummary_doubleSpinBox.setValue(self.solvent_rayleighLength*1E3) # [mm]
                
                elif stype == "OA":
                    self.solventOA_n2Summary_doubleSpinBox.setValue(self.solvent_n2*1E22) # display n2 in multiples of 10^-9 cm^2/GW 
                    self.solventOA_n2Summary_doubleSpinBox.setValue(self.solvent_n2*1E22) # display n2 in multiples of 10^-9 cm^2/GW 
                    #self.solventOA_TSummary_doubleSpinBox.setValue(self.solventOA_T)            
                    # TEMPORARILY DISABLED
                    #self.solventOA_betaSummary_doubleSpinBox.setValue(self.solvent_beta) # CUVETTE_PATH_LENGTH is the solvent/sample thickness assuming no one-photon absorption
        
            case "Sample":
                pass

    def get_general_parameters(self):
        """Gets the values from 'General Parameters' GUI frame and assigns them to variables with basic SI units (mm -> m):\n
        `self.l_silica`, `self.lda`, `self.z_range`, `self.ra`, `self.d0`\n
        and then calculate parameters that depend only on them and not on fitting process:\n
        `self.silica_n2`
        """
        self.l_silica = self.silicaThickness_dataFittingTab_doubleSpinBox.value()*1E-3 # [m] silica thickness
        self.lda = self.wavelength_dataFittingTab_doubleSpinBox.value()*1E-9 # [m] wavelength
        self.z_range = self.zscanRange_doubleSpinBox.value()*1E-3 # [m] z-scan range
        self.ra = self.apertureDiameter_doubleSpinBox.value()/2*1E-3 # [m] CA aperture radius
        self.d0 = self.apertureToFocusDistance_doubleSpinBox.value()*1E-3 # [m] distance from focal point to aperture plane
        self.silica_n2 = 2.8203E-20 - 3E-27/(self.lda) + 2E-33/(self.lda)**2 # [m2/W] Bandar A. Babgi formula (based on David Milam tables for n2)
    
    def get_curve_interpretation(self, ftype, stype, from_what:str, on_data_load=False):
        '''Interprets the curve based on its geometry according to some approximated formulas (manual fitting)\n
        or based on automatic fitting parameters (automatic fitting) and calls for error estimation.'''
        match ftype:
            case "Silica":
                match from_what:
                    case "from_geometry": # fit_manually
                        # Adjust centerPoint slider (prevent unexpected shifting of the slider and the curve to a side of the plot)
                        self.silica_centerPoint_slider.setMaximum(self.silica_nop-1)
                        self.silica_centerPoint_slider.setValue(int(self.silica_nop/2))
                        self.silica_centerPoint_slider.maximum()
                        
                        # read sliders
                        self.silica_zeroLevel = self.silica_zeroLevel_slider.value()/100
                        self.silica_centerPoint = np.round(self.silica_centerPoint_slider.value()-self.silica_nop/2)
                        if on_data_load:
                            self.silica_DPhi0, self.silica_beamwaist, self.silica_rayleighLength = \
                                self.params_from_geometry(ftype, "CA")
                        else:
                            self.silica_DPhi0 = self.silica_DPhi0_slider.value()*MAX_DPHI0/self.silica_DPhi0_slider.maximum()
                            self.silica_rayleighLength = self.silica_RayleighLength_slider.value()*self.z_range/2/self.silica_RayleighLength_slider.maximum()
                            self.silica_beamwaist = np.sqrt(self.silica_rayleighLength*self.lda/np.pi)
                    case "from_autofit": # fit_automatically
                        self.silica_zeroLevel = self.silica_minimizerResult.params['Zero'].value
                        self.silica_centerPoint = self.silica_minimizerResult.params['Center'].value
                        self.silica_DPhi0 = self.silica_minimizerResult.params['DPhi0'].value
                        self.silica_beamwaist = self.silica_minimizerResult.params['Beamwaist'].value
                        self.silica_rayleighLength = float(np.pi*self.silica_beamwaist**2/self.lda)
                
                self.laserI0 = self.silica_DPhi0*self.lda/(2*np.pi*self.l_silica*self.silica_n2) # [W/m2]
                
                try:
                    self.numericalAperture = self.silica_beamwaist/self.silica_rayleighLength
                except ZeroDivisionError:
                    self.numericalAperture = 0
                
                # calculate errors
                match from_what:
                    case "from_geometry": # fit_manually
                        self.silica_zeroLevel = error_rounding(
                                self.silica_zeroLevel,
                                sigfigs=3, warn=False
                                )
                        self.silica_centerPoint = error_rounding(
                                self.silica_centerPoint,
                                sigfigs=3, warn=False
                                )
                        self.silica_DPhi0 = error_rounding(self.silica_DPhi0, sigfigs=3, warn=False)
                        
                        self.silica_beamwaist = error_rounding(self.silica_beamwaist, sigfigs=3, warn=False)                            
                        
                        # self.silica_rayleighLengthError = \
                        #     (np.pi/2/self.lda*self.silica_minimizerResult.params['Beamwaist'].value*self.silica_beamwaistError) # [m] Rayleigh length
                        self.silica_rayleighLength = error_rounding(self.silica_rayleighLength, sigfigs=3, warn=False)

                        # self.laserI0Error = self.silica_DPhi0Error*self.lda/(2*np.pi*self.l_silica*self.silica_n2) # [W/m2]
                        self.laserI0 = error_rounding(self.laserI0,sigfigs=3, warn=False) # W/m2
                        
                        # self.numericalApertureError = self.silica_beamwaistError/self.silica_rayleighLength + \
                        #     self.silica_beamwaist*self.silica_rayleighLengthError/self.silica_rayleighLength**2
                        
                        try:
                            self.numericalAperture = error_rounding(self.numericalAperture, sigfigs=3, warn=False)
                        except ValueError:
                            self.numericalAperture = 0
                    
                    case "from_autofit": # fit_automatically
                        self.silica_zeroLevel, self.silica_zeroLevelError = error_rounding(
                                self.silica_minimizerResult.params['Zero'].value,
                                uncertainty=self.silica_minimizerResult.params['Zero'].stderr,
                                separation=tuple
                                )
                        self.silica_centerPoint, self.silica_centerPointError = error_rounding(
                                self.silica_minimizerResult.params['Center'].value,
                                uncertainty=self.silica_minimizerResult.params['Center'].stderr,
                                separation=tuple
                                )
                        self.silica_DPhi0, self.silica_DPhi0Error = error_rounding(
                                self.silica_minimizerResult.params['DPhi0'].value,
                                uncertainty=self.silica_minimizerResult.params['DPhi0'].stderr,
                                separation=tuple
                                )
                        self.silica_beamwaist, self.silica_beamwaistError = error_rounding(
                                self.silica_minimizerResult.params['Beamwaist'].value,
                                uncertainty=self.silica_minimizerResult.params['Beamwaist'].stderr,
                                separation=tuple
                                )
                        
                        self.silica_rayleighLengthError = \
                            (np.pi/2/self.lda*self.silica_minimizerResult.params['Beamwaist'].value*self.silica_beamwaistError) # [m] Rayleigh length
                        self.silica_rayleighLength, self.silica_rayleighLengthError = error_rounding(
                                self.silica_rayleighLength,
                                uncertainty=self.silica_rayleighLengthError,
                                separation=tuple
                                )

                        self.laserI0Error = self.silica_DPhi0Error*self.lda/(2*np.pi*self.l_silica*self.silica_n2) # [W/m2]
                        
                        self.laserI0, self.laserI0Error = error_rounding(
                            self.laserI0,
                            uncertainty=self.laserI0Error,
                            separation=tuple
                            ) # W/m2
                        

                        self.numericalApertureError = self.silica_beamwaistError/self.silica_rayleighLength + \
                            self.silica_beamwaist*self.silica_rayleighLengthError/self.silica_rayleighLength**2
                            
                        if not np.isnan(self.numericalApertureError):
                            self.numericalAperture, self.numericalApertureError = error_rounding(
                                self.numericalAperture,
                                uncertainty=self.numericalApertureError,
                                separation=tuple
                                )

            case "Solvent":
                match stype:
                    case "CA":
                        match from_what:
                            case "from_geometry":
                                self.solventCA_zeroLevel = self.solventCA_zeroLevel_slider.value()/100
                                self.solventCA_centerPoint = np.round(self.solventCA_centerPoint_slider.value()-self.solvent_nop/2)
                                if on_data_load:
                                    self.solvent_DPhi0, self.solvent_beamwaist, self.solvent_rayleighLength = \
                                        self.params_from_geometry(ftype, "CA")
                                else:
                                    self.solvent_DPhi0 = self.solventCA_DPhi0_slider.value()*MAX_DPHI0/self.solventCA_DPhi0_slider.maximum()
                                    self.solvent_rayleighLength = self.solventCA_RayleighLength_slider.value()*self.z_range/2/self.solventCA_RayleighLength_slider.maximum()
                                    self.solvent_beamwaist = np.sqrt(self.solvent_rayleighLength*self.lda/np.pi)
                                
                                if ftype == "Solvent" and self.solventCA_customBeamwaist_checkBox.isChecked() is False:
                                    self.solvent_beamwaist, self.solvent_rayleighLength = self.silica_beamwaist, self.silica_rayleighLength
                            case "from_autofit":
                                pass
                            
                        self.solvent_n2 = self.solvent_DPhi0/self.silica_DPhi0*self.silica_n2*self.l_silica/0.001 # [m2/W] 0.001 stands for beam path in 1-mm cuvette
                        
                        # calculate errors
                        match from_what:
                            case "from_geometry": # fit_manually
                                pass # don't calculate errors, they are not known
                            case "from_autofit": # fit_automatically
                                try:
                                    self.solventCA_zeroLevel, self.solventCA_zeroLevelError, self.solventCA_zeroLevelPrecision = \
                                        error_rounding(self.solventCA_minimizerResult.params['Zero'].value, self.solventCA_minimizerResult.params['Zero'].stderr)
                                    self.solventCA_centerPoint, self.solventCA_centerPointError, self.solventCA_centerPointPrecision = \
                                        error_rounding(self.solventCA_minimizerResult.params['Center'].value, self.solventCA_minimizerResult.params['Center'].stderr)
                                    self.solvent_DPhi0, self.solvent_DPhi0Error, self.solvent_DPhi0Precision = \
                                        error_rounding(self.solventCA_minimizerResult.params['DPhi0'].value, self.solventCA_minimizerResult.params['DPhi0'].stderr)
                                    if self.solventCA_minimizerResult.params['Beamwaist'].vary is True:
                                        self.solvent_beamwaist, self.solvent_beamwaistError, self.solvent_beamwaistPrecision = \
                                            error_rounding(self.solventCA_minimizerResult.params['Beamwaist'].value, self.solventCA_minimizerResult.params['Beamwaist'].stderr)
                                    else:
                                        self.solvent_beamwaist, self.solvent_beamwaistError, self.solvent_beamwaistPrecision = \
                                            self.silica_beamwaist, self.silica_beamwaistError, self.silica_beamwaistPrecision
                                    
                                    self.solvent_rayleighLengthError = (np.pi/2/self.lda*self.solventCA_minimizerResult.params['Beamwaist'].value*self.solvent_beamwaistError) # [m] Rayleigh length
                                    self.solvent_rayleighLength, self.solvent_rayleighLengthError, self.solvent_rayleighLengthPrecision = \
                                        error_rounding(self.solvent_rayleighLength, self.solvent_rayleighLengthError)
                                    
                                    self.solvent_n2Error = self.solvent_DPhi0Error/self.silica_DPhi0*self.silica_n2*self.l_silica/0.001 + \
                                        self.solvent_DPhi0*self.silica_DPhi0Error/self.silica_DPhi0**2*self.silica_n2*self.l_silica/0.001
                                    self.solvent_n2, self.solvent_n2Error, self.solvent_n2Precision = error_rounding(self.solvent_n2*1E13, self.solvent_n2Error*1E13)
                                    # recover original units of m2/W
                                    self.solvent_n2 *= 1E-13
                                    self.solvent_n2Error *= 1E-13
                                except Exception:
                                    logging.error(traceback.format_exc())
                                    self.showdialog('Error', 'The fit didn\'t converge. Try using different initial parameters.')
                    case "OA":
                        pass
                match stype:
                    case "CA":
                        match from_what:
                            case "from_geometry":
                                self.solventCA_zeroLevel = self.solventCA_zeroLevel_slider.value()/100
                                self.solventCA_centerPoint = np.round(self.solventCA_centerPoint_slider.value()-self.solvent_nop/2)
                                if on_data_load:
                                    self.solvent_DPhi0, self.solvent_beamwaist, self.solvent_rayleighLength = \
                                        self.params_from_geometry(ftype, "CA")
                                else:
                                    self.solvent_DPhi0 = self.solventCA_DPhi0_slider.value()*MAX_DPHI0/self.solventCA_DPhi0_slider.maximum()
                                    self.solvent_rayleighLength = self.solventCA_RayleighLength_slider.value()*self.z_range/2/self.solventCA_RayleighLength_slider.maximum()
                                    self.solvent_beamwaist = np.sqrt(self.solvent_rayleighLength*self.lda/np.pi)
                                
                                if ftype == "Solvent" and self.solventCA_customBeamwaist_checkBox.isChecked() is False:
                                    self.solvent_beamwaist, self.solvent_rayleighLength = self.silica_beamwaist, self.silica_rayleighLength
                            case "from_autofit":
                                pass
                            
                        self.solvent_n2 = self.solvent_DPhi0/self.silica_DPhi0*self.silica_n2*self.l_silica/0.001 # [m2/W] 0.001 stands for beam path in 1-mm cuvette
                        
                        # calculate errors
                        match from_what:
                            case "from_geometry": # fit_manually
                                pass # don't calculate errors, they are not known
                            case "from_autofit": # fit_automatically
                                try:
                                    self.solventCA_zeroLevel, self.solventCA_zeroLevelError, self.solventCA_zeroLevelPrecision = \
                                        error_rounding(self.solventCA_minimizerResult.params['Zero'].value, self.solventCA_minimizerResult.params['Zero'].stderr)
                                    self.solventCA_centerPoint, self.solventCA_centerPointError, self.solventCA_centerPointPrecision = \
                                        error_rounding(self.solventCA_minimizerResult.params['Center'].value, self.solventCA_minimizerResult.params['Center'].stderr)
                                    self.solvent_DPhi0, self.solvent_DPhi0Error, self.solvent_DPhi0Precision = \
                                        error_rounding(self.solventCA_minimizerResult.params['DPhi0'].value, self.solventCA_minimizerResult.params['DPhi0'].stderr)
                                    if self.solventCA_minimizerResult.params['Beamwaist'].vary is True:
                                        self.solvent_beamwaist, self.solvent_beamwaistError, self.solvent_beamwaistPrecision = \
                                            error_rounding(self.solventCA_minimizerResult.params['Beamwaist'].value, self.solventCA_minimizerResult.params['Beamwaist'].stderr)
                                    else:
                                        self.solvent_beamwaist, self.solvent_beamwaistError, self.solvent_beamwaistPrecision = \
                                            self.silica_beamwaist, self.silica_beamwaistError, self.silica_beamwaistPrecision
                                    
                                    self.solvent_rayleighLengthError = (np.pi/2/self.lda*self.solventCA_minimizerResult.params['Beamwaist'].value*self.solvent_beamwaistError) # [m] Rayleigh length
                                    self.solvent_rayleighLength, self.solvent_rayleighLengthError, self.solvent_rayleighLengthPrecision = \
                                        error_rounding(self.solvent_rayleighLength, self.solvent_rayleighLengthError)
                                    
                                    self.solvent_n2Error = self.solvent_DPhi0Error/self.silica_DPhi0*self.silica_n2*self.l_silica/0.001 + \
                                        self.solvent_DPhi0*self.silica_DPhi0Error/self.silica_DPhi0**2*self.silica_n2*self.l_silica/0.001
                                    self.solvent_n2, self.solvent_n2Error, self.solvent_n2Precision = error_rounding(self.solvent_n2*1E13, self.solvent_n2Error*1E13)
                                    # recover original units of m2/W
                                    self.solvent_n2 *= 1E-13
                                    self.solvent_n2Error *= 1E-13
                                except Exception:
                                    logging.error(traceback.format_exc())
                                    self.showdialog('Error', 'The fit didn\'t converge. Try using different initial parameters.')
                    case "OA":
                        pass
            case "Sample":
                if self.sampleCA_fittingLine_drawn is True:
                    pass
                
                self.sampleCA_zeroLevel = self.sampleCA_zeroLevel_slider.value()/100
                self.sampleCA_centerPoint = self.sampleCA_centerPoint_slider.value()-50
                
                if self.sampleOA_fittingLine_drawn is True:
                    pass
                
                self.sampleOA_zeroLevel = self.sampleOA_zeroLevel_slider.value()/100
                self.sampleOA_centerPoint = self.sampleOA_centerPoint_slider.value()-50
        
    def params_from_geometry(self, ftype, stype):
        '''Using equations given by van Stryland and Sheik-Bahae in
        http://www.phys.unm.edu/msbahae/publications/z-scan.pdf
        return the fitting curve physical interpretation.'''
        match ftype:
            case "Silica":
                if stype == "CA":
                    dataset = self.silica_data_set[0],self.silica_data_set[1]
                else:
                    dataset = self.silica_data_set[0],self.silica_data_set[3]
            case "Solvent":
                if stype == "CA":
                    dataset = self.solvent_data_set[0],self.solvent_data_set[1]
                else:
                    dataset = self.solvent_data_set[0],self.solvent_data_set[3]
            case "Sample":
                if stype == "CA":
                    dataset = self.sample_data_set[0],self.sample_data_set[1]
                else:
                    dataset = self.sample_data_set[0],self.sample_data_set[3]
            
        fit_x, fit_y = dataset[0],dataset[1]
        
        if stype == "CA":
            try:
                ymax_pos = np.where(fit_y==np.max(fit_y))[0][0]
                ymin_pos = np.where(fit_y==np.min(fit_y))[0][0]
                
                if ymax_pos > ymin_pos:
                    deltaTpv_sign = 1
                elif ymax_pos < ymin_pos:
                    deltaTpv_sign = -1
                else:
                    deltaTpv_sign = 0
                
                deltaTpv = deltaTpv_sign*abs(np.max(fit_y)-np.min(fit_y))
                deltaPhi0 = deltaTpv/0.406
                deltaZpv = abs(fit_x[ymax_pos]-fit_x[ymin_pos])*1E-3 # [m]
                rayleighLength = deltaZpv/1.7 # [m] Rayleigh length
                beamwaist = np.sqrt(rayleighLength*self.lda/np.pi) # [m] beam radius in focal point
                return deltaPhi0, beamwaist, rayleighLength
            
            except TypeError:
                print('Something is wrong while interpreting the closed aperture curve.')
        
        elif stype == "OA":
            pass
            
    def set_sliders_positions(self, ftype, stype):
        ''' Called by `self.fit_automatically`.\n
        This method recalculates value to reflect given slider position.
        '''
        match ftype:
            case "Silica":
                # Silica CA sliders
                # Because change in slider value triggers the 'fit_manually' method, it has to be disabled.
                self.silica_RayleighLength_slider.valueChanged.disconnect()
                self.silica_centerPoint_slider.valueChanged.disconnect()
                self.silica_zeroLevel_slider.valueChanged.disconnect()
                self.silica_DPhi0_slider.valueChanged.disconnect()
                
                self.silica_RayleighLength_slider.setValue(int(round(self.silica_rayleighLength*self.silica_RayleighLength_slider.maximum()/(self.z_range/2))))
                self.silica_centerPoint_slider.setValue(int(round(self.silica_centerPoint+self.silica_centerPoint_slider.maximum()/2)))
                self.silica_zeroLevel_slider.setValue(int(round(self.silica_zeroLevel*100)))
                self.silica_DPhi0_slider.setValue(int(round(self.silica_DPhi0/np.pi*self.silica_DPhi0_slider.maximum())))
                self.silica_DPhi0_slider.setValue(int(round(self.silica_DPhi0/np.pi*self.silica_DPhi0_slider.maximum())))

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

                    if self.solventCA_customBeamwaist_checkBox.isChecked() is False:
                        self.solventCA_RayleighLength_slider.setValue(self.silica_RayleighLength_slider.value())
                    else:
                        self.solventCA_RayleighLength_slider.setValue(int(round(self.solvent_rayleighLength*self.solventCA_RayleighLength_slider.maximum()/(self.z_range/2))))
                    self.solventCA_centerPoint_slider.setValue(int(round(self.solventCA_centerPoint+self.solventCA_centerPoint_slider.maximum()/2)))
                    self.solventCA_zeroLevel_slider.setValue(int(round(self.solventCA_zeroLevel*100)))
                    self.solventCA_DPhi0_slider.setValue(int(round(self.solvent_DPhi0/np.pi*self.solventCA_DPhi0_slider.maximum())))
                    self.solventCA_DPhi0_slider.setValue(int(round(self.solvent_DPhi0/np.pi*self.solventCA_DPhi0_slider.maximum())))

                    # And now when all is updated by the 'fit_automatically', reconnect the sliders to their slots
                    self.solventCA_RayleighLength_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="CA"))
                    self.solventCA_centerPoint_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="CA"))
                    self.solventCA_zeroLevel_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="CA"))
                    self.solventCA_DPhi0_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="CA"))
                
                elif stype == "OA":
                    # Solvent OA sliders
                    # Because change in slider value triggers the 'fit_manually' method, it has to be disabled.
                    self.solventOA_centerPoint_slider.valueChanged.disconnect()
                    self.solventOA_zeroLevel_slider.valueChanged.disconnect()
                    self.solventOA_T_slider.valueChanged.disconnect()

                    self.solventOA_centerPoint_slider.setValue(int(round(self.solventOA_centerPoint+50)))
                    self.solventOA_zeroLevel_slider.setValue(int(round(self.solventOA_zeroLevel*100)))
                    self.solventOA_T_slider.setValue(int(round(self.solventOA_T*self.solventOA_T_slider.maximum()/SOLVENT_T_SLIDER_MAX)))

                    # And now when all is updated by the 'fit_automatically', reconnect the sliders to their slots
                    self.solventOA_centerPoint_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="OA"))
                    self.solventOA_zeroLevel_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="OA"))
                    self.solventOA_T_slider.valueChanged.connect(lambda: self.fit_manually(ftype="Solvent", stype="OA"))

            case "Sample":
                pass
        
    def data_display(self, data_set, ftype:str):
        '''Called by data_loader(). Displays:\n
        1) CA column divided by OA column to remove influence of OA on CA as 'CA'
        2) OA column as 'OA'\n
        to separate data processing for n2 and beta nonlinear coefficients.'''
        self.get_general_parameters()

        def basic_data_manipulation(data) -> Tuple[NDArray,list,list,list]:
            """Reference division, separation of OA signal from CA fluctuation signal, normalization

            Args:
                data (any): four-column data containing datapoints, CA data, Reference data and OA data

            Returns:
                Tuple: Data ready for separated determination of n2 and beta nonlinear coefficients
            """            
            positions, ca_data0, ref_data0, oa_data0 = data
            # Divide data by reference
            ca_data = [ca/ref for ca, ref in zip(ca_data0, ref_data0)]
            #ref_data = [ref/ref for ref in ref_data0]
            oa_data = [oa/ref for oa, ref in zip(oa_data0, ref_data0)]

            # centralize positions and normalize by z-scan range
            nop = len(positions)
            positions = np.array([(self.z_range*zz/nop-self.z_range/2)*1000 for zz in range(nop)]) # express in mm
            
            def normalize(data):# Normalize to 1 (edge of the data as reference)
                zero_lvl = (np.mean(data[0:10])+np.mean(data[len(data)-10:]))/2
                return data/zero_lvl
            
            ca_data = normalize(ca_data)
            oa_data = normalize(oa_data)
            
            return positions, ca_data, ref_data0, oa_data
        
        self.positions, self.ca_data, self.ref_data, self.oa_data = basic_data_manipulation(data_set)

        # Create reference-corrected 'data_set' variable for each ftype for further reference
        match ftype:
            case "Silica":
                self.silica_data_set = [self.positions, self.ca_data, self.ref_data, self.oa_data]
                self.silica_nop = len(self.positions)
            case "Solvent":
                self.solvent_data_set = [self.positions, self.ca_data, self.ref_data, self.oa_data]
                self.solvent_nop = len(self.positions)
            case "Sample":
                self.sample_data_set = [self.positions, self.ca_data, self.ref_data, self.oa_data]
                self.sample_nop = len(self.positions)
                
        # Display data on proper figures
        # Use reference-corrected data for fitting
        match ftype:
            case "Silica":
                self.update_datafitting_plotlimits(self.silica_data_set,ftype=ftype)

            case "Solvent":
                self.update_datafitting_plotlimits(self.solvent_data_set,ftype=ftype)

            case "Sample":
                self.update_datafitting_plotlimits(self.sample_data_set,ftype=ftype)
    
    def set_new_positions(self):
        """Takes current values in 'General Parameters' from GUI and updates data_set[0] (positions) to meet new Z-scan range.\n
        Updates fitting plots limits to visualize the change.
        """
        self.get_general_parameters()

        if hasattr(self,'silica_data_set'):
            self.silica_data_set[0] = np.array([(self.z_range*zz/self.silica_nop-self.z_range/2)*1000 for zz in range(self.silica_nop)]) # [mm] update positions with newly-read Z-scan range value
            self.update_datafitting_plotlimits(self.silica_data_set,ftype='Silica')

        if hasattr(self,'solvent_data_set'):
            self.solvent_data_set[0] = np.array([(self.z_range*zz/self.solvent_nop-self.z_range/2)*1000 for zz in range(self.solvent_nop)]) # [mm] update positions with newly-read Z-scan range value
            self.update_datafitting_plotlimits(self.solvent_data_set,ftype='Solvent')

        if hasattr(self,'sample_data_set'):
            self.sample_data_set[0] = np.array([(self.z_range*zz/self.sample_nop-self.z_range/2)*1000 for zz in range(self.sample_nop)]) # [mm] update positions with newly-read Z-scan range value
            self.update_datafitting_plotlimits(self.sample_data_set,ftype='Sample')

    def update_datafitting_plotlimits(self, data_set:list, ftype:str) -> None:
        """Updates limits of the data fitting plots for given `ftype` (Silica, Solvent, Sample)

        Args:
            data_set (list): four-column reference-corrected data
            ftype (str): parameter holding information of sample type (Silica, Solvent, Sample)
        """        
        self.get_general_parameters() # Ensures working with currently typed-in General Parameters values from GUI
        padding_vertical = 0.01
        
        def set_limits(axes, line, direction, padding):
            if direction == "vertical":
                axes.set_ylim(top=np.max(line.get_ydata())*(1+padding),bottom=np.min(line.get_ydata())*(1-padding))
        
        match ftype:
            case "Silica":
                # Closed aperture
                line = self.silica_figure.axes.get_lines()[0]
                line.set_xdata(data_set[0])
                #line.set_ydata(ca_data)
                line.set_ydata([ca/oa for ca, oa in zip(data_set[1], data_set[3])])
                
                self.silica_figure.axes.set_xlim(left=-self.z_range/2*1000, right=self.z_range/2*1000) # displayed in mm
                set_limits(self.silica_figure.axes, line, "vertical", padding_vertical)

                # Open aperture
                line = self.silicaOA_figure.axes.get_lines()[0]
                line.set_xdata(data_set[0])
                line.set_ydata(data_set[3])
                
                self.silicaOA_figure.axes.set_xlim(left=-self.z_range/2*1000, right=self.z_range/2*1000) # displayed in mm
                set_limits(self.silicaOA_figure.axes, line, "vertical", padding_vertical)
                #self.silicaOA_figure.axes.set_ylim(top=np.max(line.get_ydata())*(1+margin_vertical),bottom=np.min(line.get_ydata())*(1-margin_vertical))
                
                # Update
                self.silica_figure.axes.relim()
                self.silicaOA_figure.axes.relim()
                
                self.silica_figure.axes.autoscale_view()
                self.silicaOA_figure.axes.autoscale_view()
                
                self.silica_figure.draw_idle()
                self.silicaOA_figure.draw_idle()
            
            case "Solvent":
                # Closed aperture
                line = self.solventCA_figure.axes.get_lines()[0]
                line.set_xdata(data_set[0])
                #line.set_ydata(ca_data)
                line.set_ydata([ca/oa for ca, oa in zip(data_set[1], data_set[3])])

                self.solventCA_figure.axes.set_xlim(left=-self.z_range/2*1000, right=self.z_range/2*1000) # displayed in mm
                set_limits(self.solventCA_figure.axes, line, "vertical", padding_vertical)
                #self.solventCA_figure.axes.set_ylim(top=np.max(line.get_ydata())*(1+margin_vertical),bottom=np.min(line.get_ydata())*(1-margin_vertical))

                # Open aperture
                line = self.solventOA_figure.axes.get_lines()[0]
                line.set_xdata(data_set[0])
                line.set_ydata(data_set[3])

                self.solventOA_figure.axes.set_xlim(left=-self.z_range/2*1000, right=self.z_range/2*1000) # displayed in mm
                set_limits(self.solventOA_figure.axes, line, "vertical", padding_vertical)
                #self.solventOA_figure.axes.set_ylim(top=np.max(line.get_ydata())*(1+margin_vertical),bottom=np.min(line.get_ydata())*(1-margin_vertical))

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
                line.set_xdata(data_set[0])
                #line.set_ydata(ca_data)
                line.set_ydata([ca/oa for ca, oa in zip(data_set[1], data_set[3])])

                self.sampleCA_figure.axes.set_xlim(left=-self.z_range/2*1000, right=self.z_range/2*1000) # displayed in mm
                set_limits(self.sampleCA_figure.axes, line, "vertical", padding_vertical)
                #self.sampleCA_figure.axes.set_ylim(top=np.max(line.get_ydata())*1.04,bottom=np.min(line.get_ydata())*0.96)
            
                # Open aperture
                line = self.sampleOA_figure.axes.get_lines()[0]
                line.set_xdata(data_set[0])
                line.set_ydata(data_set[3])

                self.sampleOA_figure.axes.set_xlim(left=-self.z_range/2*1000, right=self.z_range/2*1000) # displayed in mm
                set_limits(self.sampleOA_figure.axes, line, "vertical", padding_vertical)
                #self.sampleOA_figure.axes.set_ylim(top=np.max(line.get_ydata())*1.04,bottom=np.min(line.get_ydata())*0.96)
            
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
                    filter_size = self.silica_filterSize_slider.value()
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
                            line = self.silica_figure.axes.get_lines()[0]
                            line.set_ydata(filtered_y)
                            
                            self.silica_figure.axes.relim()
                            self.silica_figure.axes.autoscale_view()
                            self.silica_figure.draw_idle()
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
                self.showdialog('Error',
                'Possibly too few data points!\nMinimum required is 16 datapoints.')

    def enable_cursors(self, ftype:str, stype:str) -> None:
        match ftype:
            case "Silica":
                if stype == "CA":
                    if self.silica_fixROI_checkBox.isChecked() is True:
                        self.silica_cursor_positioner = BlittedCursor(self.silica_figure.axes, color = 'magenta', linewidth = 2)
                        self.on_mouse_move = self.silica_figure.mpl_connect('motion_notify_event', self.silica_cursor_positioner.on_mouse_move)
                        self.on_mouse_click = self.silica_figure.mpl_connect('button_press_event', lambda event: self.collect_cursor_clicks(event,ftype))
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
            case "Solvent":
                if stype == "CA":
                    if self.solventCA_fixROI_checkBox.isChecked() is True:
                        self.solventCA_cursor_positioner = BlittedCursor(self.solventCA_figure.axes, color = 'magenta', linewidth = 2)
                        self.on_mouse_move = self.solventCA_figure.mpl_connect('motion_notify_event', self.solventCA_cursor_positioner.on_mouse_move)
                        self.on_mouse_click = self.solventCA_figure.mpl_connect('button_press_event', lambda event: self.collect_cursor_clicks(event,ftype,stype))
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
                        self.solventOA_cursor_positioner = BlittedCursor(self.solventOA_figure.axes, color = 'magenta', linewidth = 2)
                        self.on_mouse_move = self.solventOA_figure.mpl_connect('motion_notify_event', self.solventOA_cursor_positioner.on_mouse_move)
                        self.on_mouse_click = self.solventOA_figure.mpl_connect('button_press_event', lambda event: self.collect_cursor_clicks(event,ftype,stype))
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
            
            case "Sample":
                if stype == "CA":
                    self.sampleCA_cursor_positioner = BlittedCursor(self.sampleCA_figure.axes, color = 'magenta', linewidth = 2)
                elif stype == "OA":
                    self.sampleOA_cursor_positioner = BlittedCursor(self.sampleOA_figure.axes, color = 'magenta', linewidth = 2)
    
    def collect_cursor_clicks(self, event, ftype:str, stype="") -> None:
        x, y = event.xdata, event.ydata
        match ftype:
            case "Silica":
                if not hasattr(self, "silica_cursorPositions"):
                    self.silica_cursorPositions = []
                if len(self.silica_cursorPositions)<2:
                    self.silica_cursorPositions.append((x, y))
                    if len(self.silica_cursorPositions) == 1:
                        self.silica_verline1 = self.silica_figure.axes.axvline(x, color="orange", linewidth=2)
                    else:
                        self.silica_verline2 = self.silica_figure.axes.axvline(x, color="orange", linewidth=2)
                else:
                    self.silica_cursorPositions = []
                    self.silica_cursorPositions.append((x, y))
                    self.silica_verline1.set_xdata(x)
                    self.silica_verline2.set_xdata(None)
                
                self.silica_figure.draw_idle()
            
            case "Solvent":
                if stype == "CA":
                    if not hasattr(self, "solventCA_cursorPositions"):
                        self.solventCA_cursorPositions = []
                    if len(self.solventCA_cursorPositions)<2:
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
                    if len(self.solventOA_cursorPositions)<2:
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
            
            case "Sample":
                pass

    def draw_fitting_line(self, ftype:str, stype:str) -> None:
        match ftype:
            case 'Silica':
                if self.silica_fittingLine_drawn is False:
                    self.silica_fitting_line_ca, = self.silica_figure.axes.plot(self.silica_data_set[0],self.result,'r')
                    self.silica_figure.draw_idle()

                    self.silica_fittingLine_drawn = True
                else:
                    if hasattr(self,"silica_data_set"):
                        self.silica_fitting_line_ca.set_xdata(self.silica_data_set[0])
                        self.silica_fitting_line_ca.set_ydata(self.result)
                        
                        self.silica_figure.axes.set_xlim(left=-self.z_range/2*1000, right=self.z_range/2*1000) # displayed in mm
                        self.silica_figure.axes.relim()
                        self.silica_figure.axes.autoscale_view()
                        self.silica_figure.draw_idle()
            
            case 'Solvent':
                if stype == "CA":
                    if self.solventCA_fittingLine_drawn is False:
                        self.solvent_fitting_line_ca, = self.solventCA_figure.axes.plot(self.solvent_data_set[0],self.result,'r')
                        self.solventCA_figure.draw_idle()

                        self.solventCA_fittingLine_drawn = True
                    else:
                        if hasattr(self,"solvent_data_set"):
                            self.solvent_fitting_line_ca.set_xdata(self.solvent_data_set[0])
                            self.solvent_fitting_line_ca.set_ydata(self.result)
                            
                            self.solventCA_figure.axes.set_xlim(left=-self.z_range/2*1000, right=self.z_range/2*1000) # displayed in mm
                            self.solventCA_figure.axes.relim()
                            self.solventCA_figure.axes.autoscale_view()
                            self.solventCA_figure.draw_idle()
                
                elif stype == "OA":
                    if self.solventOA_fittingLine_drawn is False:
                        self.solvent_fitting_line_oa, = self.solventOA_figure.axes.plot(self.solvent_data_set[0],self.result,'r')
                        self.solventOA_figure.draw_idle()

                        self.solventOA_fittingLine_drawn = True
                    else:
                        if hasattr(self,"solvent_data_set"):
                            self.solvent_fitting_line_oa.set_xdata(self.solvent_data_set[0])
                            self.solvent_fitting_line_oa.set_ydata(self.result)
                            
                            self.solventOA_figure.axes.set_xlim(left=-self.z_range/2*1000, right=self.z_range/2*1000) # displayed in mm
                            self.solventOA_figure.axes.relim()
                            self.solventOA_figure.axes.autoscale_view()
                            self.solventOA_figure.draw_idle()
            
            case 'Sample':
                if stype == "CA":
                    if self.sampleCA_fittingLine_drawn is False:
                        self.sample_fitting_line_ca, = self.sampleCA_figure.axes.plot(self.sample_data_set[0],self.result,'r')
                        self.sampleCA_figure.draw_idle()

                        self.sampleCA_fittingLine_drawn = True
                    else:
                        if hasattr(self,"sample_data_set"):
                            self.sample_fitting_line_ca.set_xdata(self.sample_data_set[0])
                            self.sample_fitting_line_ca.set_ydata(self.result)
                            
                            self.sampleCA_figure.axes.set_xlim(left=-self.z_range/2*1000, right=self.z_range/2*1000) # displayed in mm
                            self.sampleCA_figure.axes.relim()
                            self.sampleCA_figure.axes.autoscale_view()
                            self.sampleCA_figure.draw_idle()
                
                elif stype == "OA":
                    if self.sampleOA_fittingLine_drawn is False:
                        self.sample_fitting_line_oa, = self.sampleOA_figure.axes.plot(self.sample_data_set[0],self.result,'r')
                        self.sampleOA_figure.draw_idle()

                        self.sampleOA_fittingLine_drawn = True
                    else:
                        if hasattr(self,"sample_data_set"):
                            self.sample_fitting_line_oa.set_xdata(self.sample_data_set[0])
                            self.sample_fitting_line_oa.set_ydata(self.result)
                            
                            self.sampleOA_figure.axes.set_xlim(left=-self.z_range/2*1000, right=self.z_range/2*1000) # displayed in mm
                            self.sampleOA_figure.axes.relim()
                            self.sampleOA_figure.axes.autoscale_view()
                            self.sampleOA_figure.draw_idle()

    def fit_automatically(self, ftype:str, stype:str):
        self.get_general_parameters()
        self.get_curve_interpretation(ftype,stype,'from_geometry')
        
        match ftype:
            case "Silica":
                if stype == "CA":
                    # Datapoints for the curve to be fitted to
                    line = self.silica_figure.axes.get_lines()[0]
                else:
                    return

                line_data = line.get_data()
                
                self.silica_calculation = Fitting(self.silica_curves, self.silica_DPhi0, self.silica_beamwaist, self.silica_zeroLevel, self.silica_centerPoint,self.silica_nop,line_data[1])
                self.silica_calculation = Fitting(self.silica_curves, self.silica_DPhi0, self.silica_beamwaist, self.silica_zeroLevel, self.silica_centerPoint,self.silica_nop,line_data[1])
                minimizer_result, self.result = self.silica_calculation.automatic(self.z_range, ftype, stype, line_data)
                # Make sure that the result doesn't contain NoneTypes
                for key in minimizer_result.params.keys():
                    if minimizer_result.params[key].stderr is None:
                        minimizer_result.params[key].stderr = 0
                self.silica_minimizerResult = minimizer_result
                self.silica_autofit_done = True
                self.draw_fitting_line(ftype, stype)
                self.get_curve_interpretation(ftype,stype,'from_autofit')
                self.set_fit_summary(ftype, stype, caller="auto")
                # Use these exact number from error-corrected fit parameters to set sliders to their positions
                self.set_sliders_positions(ftype, stype)
                
            case "Solvent":
                if self.silica_autofit_done is False:
                    self.showdialog('Info','Fit silica first.')
                    return
                
                if stype == "CA":
                    # Data to be fitted
                    line = self.solventCA_figure.axes.get_lines()[0]
                elif stype == "OA":
                    line = self.solventOA_figure.axes.get_lines()[0]

                line_data = line.get_data()
                
                self.solvent_nop = len(self.solvent_data_set[0])
                self.solvent_data_set[0] = np.array([(self.z_range*zz/self.solvent_nop-self.z_range/2)*1000 for zz in range(self.solvent_nop)]) # [mm] update positions with newly-read Z-scan range value
                self.solvent_nop = len(self.solvent_data_set[0])
                self.solvent_data_set[0] = np.array([(self.z_range*zz/self.solvent_nop-self.z_range/2)*1000 for zz in range(self.solvent_nop)]) # [mm] update positions with newly-read Z-scan range value
                
                self.solvent_calculation = Fitting(self.solvent_curves, self.solvent_DPhi0, self.solvent_beamwaist, self.solventCA_zeroLevel, self.solventCA_centerPoint,self.solvent_nop,line_data[1])
                self.solvent_calculation = Fitting(self.solvent_curves, self.solvent_DPhi0, self.solvent_beamwaist, self.solventCA_zeroLevel, self.solventCA_centerPoint,self.solvent_nop,line_data[1])
                minimizer_result, self.result = self.solvent_calculation.automatic(self.z_range, ftype, stype, line_data)
                match stype:
                    case "CA":
                        self.solventCA_minimizerResult = minimizer_result
                        self.solventCA_autofit_done = True
                    case "OA":
                        self.solventOA_minimizerResult = minimizer_result
                        self.solventOA_autofit_done = True
                self.draw_fitting_line(ftype, stype)
                self.get_curve_interpretation(ftype,stype,'from_autofit')
                match stype:
                    case "CA":
                        self.solventCA_minimizerResult = minimizer_result
                        self.solventCA_autofit_done = True
                    case "OA":
                        self.solventOA_minimizerResult = minimizer_result
                        self.solventOA_autofit_done = True
                self.draw_fitting_line(ftype, stype)
                self.get_curve_interpretation(ftype,stype,'from_autofit')
                self.set_fit_summary(ftype, stype, caller="auto")
                # Use these exact number from error-corrected fit parameters to set sliders to their positions
                self.set_sliders_positions(ftype, stype)
                # Use these exact number from error-corrected fit parameters to set sliders to their positions
                self.set_sliders_positions(ftype, stype)
            case "Sample":
                if self.silica_autofit_done is False:
                    self.showdialog('Info','Fit silica first.')
                    
                    if self.solventCA_autofit_done is False:
                        self.showdialog('Info','Fit solvent first.')
                        return
                    else:
                        return

                pass

    def fit_manually(self, ftype:str, stype:str) -> None:
        """Triggered by loading the data (from experiment or from file) or by fitting sliders value change.

        Args:
            ftype (str): `Silica`, `Solvent`, `Sample`
            stype (str): `CA`, `OA`
        """
        #Retrieve "General parameters"
        self.get_general_parameters()
        self.get_curve_interpretation(ftype,stype,'from_geometry')
        
        match ftype:
            case "Silica":
                if stype == "CA":
                    # Datapoints to fit the curve to
                    line = self.silica_figure.axes.get_lines()[0]
                    line_data = line.get_ydata()
                    
                    self.silica_curves = Integration(SILICA_BETA,self.silica_n2,self.silica_DPhi0,self.silica_data_set[0],self.d0,self.ra,self.lda,self.silica_beamwaist,N_COMPONENTS,INTEGRATION_STEPS)
                    self.silica_calculation = Fitting(self.silica_curves, self.silica_DPhi0, self.silica_beamwaist, self.silica_zeroLevel, self.silica_centerPoint,len(self.silica_data_set[0]),line_data)
                    self.result = self.silica_calculation.manual(self.silica_zeroLevel, self.silica_centerPoint, self.silica_DPhi0, self.silica_beamwaist, self.z_range)
                    self.silica_autofit_done = False
            case "Solvent":
                # Data to be fitted
                if stype == "CA":
                    line = self.solventCA_figure.axes.get_lines()[0]
                    line_data = line.get_ydata()
                                    
                    self.solvent_curves = Integration(0,self.solvent_n2,self.solvent_DPhi0,self.solvent_data_set[0],self.d0,self.ra,self.lda,self.solvent_beamwaist,N_COMPONENTS,INTEGRATION_STEPS) # solvent_beta = 0
                    self.solvent_calculation = Fitting(self.solvent_curves, self.solvent_DPhi0, self.solvent_beamwaist, self.solventCA_zeroLevel, self.solventCA_centerPoint,len(self.solvent_data_set[0]),line_data)
                    self.result = self.solvent_calculation.manual(self.solventCA_zeroLevel, self.solventCA_centerPoint, self.solvent_DPhi0, self.solvent_beamwaist, self.z_range)
                    self.solvent_curves = Integration(0,self.solvent_n2,self.solvent_DPhi0,self.solvent_data_set[0],self.d0,self.ra,self.lda,self.solvent_beamwaist,N_COMPONENTS,INTEGRATION_STEPS) # solvent_beta = 0
                    self.solvent_calculation = Fitting(self.solvent_curves, self.solvent_DPhi0, self.solvent_beamwaist, self.solventCA_zeroLevel, self.solventCA_centerPoint,len(self.solvent_data_set[0]),line_data)
                    self.result = self.solvent_calculation.manual(self.solventCA_zeroLevel, self.solventCA_centerPoint, self.solvent_DPhi0, self.solvent_beamwaist, self.z_range)
                    self.solventCA_autofit_done = False
                # elif stype == "OA":
                #     line = self.solventOA_figure.axes.get_lines()[0]
                #     line_data = line.get_ydata()

                #     self.solvent_curves = Integration(self.solventOA_T_doubleSpinBox.value(),self.solvent_n2,0,self.solvent_data_set[0],self.d0,self.ra,self.lda,self.solvent_beamwaist,N_COMPONENTS,INTEGRATION_STEPS, stype="OA")
                #     self.solvent_calculation = Fitting(self.solvent_curves, 0, self.solvent_beamwaist, self.solventOA_zeroLevel, self.solventOA_centerPoint,len(self.solvent_data_set[0]),line_data)
                #     self.result = self.solvent_calculation.manual(self.solventOA_zeroLevel, self.solventOA_centerPoint, 0, self.solvent_beamwaist, self.z_range, stype)
                #     self.solvent_curves = Integration(self.solventOA_T_doubleSpinBox.value(),self.solvent_n2,0,self.solvent_data_set[0],self.d0,self.ra,self.lda,self.solvent_beamwaist,N_COMPONENTS,INTEGRATION_STEPS, stype="OA")
                #     self.solvent_calculation = Fitting(self.solvent_curves, 0, self.solvent_beamwaist, self.solventOA_zeroLevel, self.solventOA_centerPoint,len(self.solvent_data_set[0]),line_data)
                #     self.result = self.solvent_calculation.manual(self.solventOA_zeroLevel, self.solventOA_centerPoint, 0, self.solvent_beamwaist, self.z_range, stype)
                
                #     self.draw_fitting_line(ftype,stype)

                #     self.solventOA_fittingDone = False

            case "Sample":
                pass
        
        self.draw_fitting_line(ftype,stype)
        self.get_curve_interpretation(ftype,stype,'from_geometry')
        
        caller = "manual"
        self.set_fit_summary(ftype, stype, caller)
        
    def read_header_params(self, caller:str, ftype:str):
        if caller == "Current Measurement":
            # General parameters
            if self.customWavelength_checkBox.isChecked() is False:
                self.wavelength_dataFittingTab_doubleSpinBox.setValue(self.wavelength_dataSavingTab_doubleSpinBox.value())
            if self.customZscanRange_checkBox.isChecked() is False:
                self.zscanRange_doubleSpinBox.setValue(np.abs(self.endPos_doubleSpinBox.value()-self.startPos_doubleSpinBox.value()))
            if self.customSilicaThickness_checkBox.isChecked() is False:
                self.silicaThickness_dataFittingTab_doubleSpinBox.setValue(self.silicaThickness_dataSavingTab_doubleSpinBox.value())
            
            if ftype == "Sample":
                self.concentration_dataFittingTab_doubleSpinBox.setValue(self.concentration_dataSavingTab_doubleSpinBox.value())

        elif caller == "Load From File":
            # Get parameters from the file header
            if len(self.header) != 0:
                for hl in self.header:
                    if hl.find("Wavelength") != -1:
                        wavelength_match = re.search(r'(([1-9][0-9]*\.?[0-9]*)|(\.[0-9]+))([Ee][+-]?[0-9]+)?',hl)
                        wavelength = float(wavelength_match.groups()[0])
                        self.wavelength_dataFittingTab_doubleSpinBox.setValue(wavelength)
                    
                    #if ftype == "Silica" and self.customZscanRange_checkBox.isChecked() is False and hl.find("Starting pos") != -1:
                    if self.customZscanRange_checkBox.isChecked() is False and hl.find("Starting pos") != -1: # read zscanRange for any ftype
                        starting_pos_match = re.search(r'(([0-9][0-9]*\.?[0-9]*)|(\.[0-9]+))([Ee][+-]?[0-9]+)?',hl)
                        starting_pos = float(starting_pos_match.groups()[0])
                        next_line = self.header[self.header.index(hl)+1]
                        end_pos_match = re.search(r'(([0-9][0-9]*\.?[0-9]*)|(\.[0-9]+))([Ee][+-]?[0-9]+)?',next_line)
                        end_pos = float(end_pos_match.groups()[0])
                        self.zscanRange_doubleSpinBox.setValue(np.abs(end_pos-starting_pos))
                    
                    if ftype == "Sample" and hl.find("Concentration") != -1:
                        concentration_match = re.search(r'(([0-9][0-9]*\.?[0-9]*\s*%)|(\.[0-9]()\s*%+))([Ee][+-]?[0-9]()\s*%+)?',hl) # Here make sure there is % symbol
                        self.concentr_percent = float(concentration_match.groups()[0].replace(' ','')[:-1]) # remove redundant space and % symbol and change to float
                        self.concentration_dataFittingTab_doubleSpinBox.setValue(self.concentr_percent)

    def solvent_autocomplete(self):
        self.solventDensity_lineEdit.setText(str(self.solvents[self.solventName_comboBox.currentText()]["density"])+' g/cm3')
        self.solventRefrIdx_lineEdit.setText(str(self.solvents[self.solventName_comboBox.currentText()]["index"]))
    
# THREAD CONTROLS
    def print_output(self, returned_value):
        print(returned_value)
    
    def thread_complete(self):
        print("THREAD COMPLETE!")
    
    def thread_it(self, func_to_execute):
        # Pass the function to execute
        worker = Worker(func_to_execute) # Any other args, kwargs are passed to the run function
        worker.signals.result.connect(self.print_output)
        worker.signals.finished.connect(self.thread_complete)

        if func_to_execute is self.mpositioner.movetostart:
            worker.signals.progress.connect(self.create_raw_log_line)

        # Execute
        self.threadpool.start(worker)

        return worker

# DIALOG BOXES
    def showdialog(self, msg_type:str, message:str):
        '''
        Message type (msg_type) is one of these: 'Error', 'Warning', 'Info'
        '''
        dialog_args = (msg_type, message)
        match msg_type:
            case "Error":
                button = QMessageBox.critical(self, *dialog_args)
            case "Warning":
                button = QMessageBox.warning(self, *dialog_args)
            case "Info":
                button = QMessageBox.information(self, *dialog_args)
        
        if button == QMessageBox.Ok:
           pass

# COLOR THEMES
    @QtCore.pyqtSlot()
    def changeSkinDark(self):
        # PALETTE
        dark_palette = QPalette()
        
        # Active
        dark_palette.setColor(QPalette.Window, QColor(35, 35, 40))
        dark_palette.setColor(QPalette.WindowText, QColor(200,200,200))
        dark_palette.setColor(QPalette.Base, QColor(60, 60, 65))
        dark_palette.setColor(QPalette.AlternateBase, QColor(35, 35, 40))
        #dark_palette.setColor(QPalette.ToolTipBase, QColor(255,255,255))
        #dark_palette.setColor(QPalette.ToolTipText, QColor(255,255,255))
        dark_palette.setColor(QPalette.Text, QColor(200,200,200))
        dark_palette.setColor(QPalette.Button, QColor(35, 35, 40))
        dark_palette.setColor(QPalette.ButtonText, QColor(200,200,200))
        #dark_palette.setColor(QPalette.BrightText, QColor(255,0,0))
        #dark_palette.setColor(QPalette.Link, QColor(42, 130, 218))
        dark_palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
        
        # Disabled
        dark_palette.setColor(QPalette.Disabled, QPalette.ButtonText, QColor(100,100,100))

        QtGui.QGuiApplication.setPalette(dark_palette)
        
        # Measurement Tab Charts
        self.rms_text.set_color("white")
        self.rms_text.set_bbox(dict(boxstyle='round', facecolor=(60/255,60/255,65/255), edgecolor=(200/255,200/255,200/255), alpha=1))

        for chart in self.charts.values():
            chart.fig.patch.set_facecolor((35/255,35/255,40/255,1))
            chart.axes.set_facecolor((35/255,35/255,40/255,1))
            
            for spine in chart.axes.spines.values():
                spine.set_color((200/255,200/255,200/255,1))
            
            chart.axes.set_xlabel(chart.axes.get_xlabel(), fontdict={'color': (200/255,200/255,200/255,1)})
            chart.axes.set_ylabel(chart.axes.get_ylabel(), fontdict={'color': (200/255,200/255,200/255,1)})
            chart.axes.tick_params(axis='both',which='both',colors=(200/255,200/255,200/255,1))
            chart.axes.grid(which="both",color=(60/255,60/255,65/255,1))
            chart.draw_idle()
        
        # Fitting Tab Charts
        for chart_types in self.fitting_charts.values():
            for chart in chart_types.values():
                chart.fig.patch.set_facecolor((35/255,35/255,40/255,1))
                chart.axes.set_facecolor((35/255,35/255,40/255,1))
                
                for spine in chart.axes.spines.values():
                    spine.set_color((200/255,200/255,200/255,1))
                
                chart.axes.set_title(chart.axes.get_title(), fontdict={'color': (200/255,200/255,200/255,1)})
                chart.axes.set_xlabel(chart.axes.get_xlabel(), fontdict={'color': (200/255,200/255,200/255,1)})
                chart.axes.set_ylabel(chart.axes.get_ylabel(), fontdict={'color': (200/255,200/255,200/255,1)})
                chart.axes.tick_params(axis='both',which='both',colors=(200/255,200/255,200/255,1))
                chart.axes.grid(which="both",color=(60/255,60/255,65/255,1))
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
        self.rms_text.set_bbox(dict(boxstyle='round', facecolor="white", edgecolor="black", alpha=1))

        for chart in self.charts.values():
            chart.fig.patch.set_facecolor((255/255,255/255,255/255,1))
            chart.axes.set_facecolor((255/255,255/255,255/255,1))

            for spine in chart.axes.spines.values():
                spine.set_color((0,0,0,1))

            chart.axes.set_xlabel(chart.axes.get_xlabel(), fontdict={'color': (0,0,0,1)})
            chart.axes.set_ylabel(chart.axes.get_ylabel(), fontdict={'color': (0,0,0,1)})
            chart.axes.tick_params(axis='both',which='both',colors=(0,0,0,1))
            chart.axes.grid(which="both",color='darkgrey')
            chart.draw_idle()
        
        # Fitting Tab Charts
        for chart_types in self.fitting_charts.values():
            for chart in chart_types.values():
                chart.fig.patch.set_facecolor((255/255,255/255,255/255,1))
                chart.axes.set_facecolor((255/255,255/255,255/255,1))
                
                for spine in chart.axes.spines.values():
                    spine.set_color((0,0,0,1))
                
                chart.axes.set_title(chart.axes.get_title(), fontdict={'color': (0,0,0,1)})
                chart.axes.set_xlabel(chart.axes.get_xlabel(), fontdict={'color': (0,0,0,1)})
                chart.axes.set_ylabel(chart.axes.get_ylabel(), fontdict={'color': (0,0,0,1)})
                chart.axes.tick_params(axis='both',which='both',colors=(0,0,0,1))
                chart.axes.grid(which="both",color='darkgrey')
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
            
## OTHER CLASSES
class Integration():
    '''Integrates the electric field according to procedure from Sheik-Bahae, given the initial parameters'''
    def __init__(self, beta, n2, DPhi0, positions, d0, aperture_radius, wavelength, beamwaist, n_components, integration_steps, stype="CA"):
        super(Integration, self).__init__()

        # Data range
        self.z = positions # evenly spaced

        # Apertures and distances
        self.d0 = d0 # distance from z=0 to the aperture plane
        self.ra = aperture_radius # in meters

        # Beam properties
        self.lda = wavelength # in meters
        self.w0 = beamwaist # in meters

        # Sample properties
        self.n2 = n2 # non-linear refractive index
        try:
            self.T = beta*self.lda/self.n2 # non-linear transmittance
        except ZeroDivisionError:
            self.T = 0
        self.DPhi0 = DPhi0 # on-axis phase shift in z=0
        
        # Integration parameters
        self.mm = n_components # number of electric field components (for Gaussian decomposition)
        self.ir = integration_steps
        
        self.stype = stype
        self.derive(self.DPhi0, self.w0, self.d0, self.ra, stype)
    
    def derive(self, DPhi0, w0, d0, ra, stype):
        '''1st method called. Called on __init__'''
        # Beam properties
        self.k = 2*np.pi/self.lda # wave vector in free space
        self.z0 = 0.5*self.k*w0**2 # diffraction length of the beam
        self.wa = w0*np.sqrt(1+d0**2/self.z0**2) # beam radius at the aperture plane

        # Aperture radius
        self.ra = ra

        # Sample properties
        self.Dphi0 = DPhi0/(1+self.z**2/self.z0**2)

        # Additional derived parameters
        self.wz = w0*np.sqrt(1+self.z**2/self.z0**2)

        self.Rz = self.z+self.z0**2/self.z
        self.d = d0-self.z
        self.g = 1+self.d/self.Rz

        match stype:
            case "CA":
                self.bigproduct()
            case "OA":
                self.calculate_Tz_for_OA(window.solventOA_absorptionModel_comboBox.currentText()) # CZYŻBY????????
        
    def calculate_Tz_for_OA(self, model):
        '''2nd method called. Called by derive() for stype = "OA".'''
        match window.fittingTabs.currentIndex():
            case 0:
                ftype = "Silica"
            case 1:
                ftype = "Solvent"
            case 2:
                ftype = "Sample"
        
        match ftype:
            case "Silica":
                return # Do not fit OA for silica
            case "Solvent":
                absorption_checkbox = window.solventOA_isAbsorption_checkBox.isChecked()
            case "Sample":
                absorption_checkbox = window.sampleOA_isAbsorption_checkBox.isChecked()
        
        if absorption_checkbox is False:
            self.Tz = 0
            # Psi_n = (n * beta_n * I0**n * L_eff)**(1/n) (n+1)-photon absorption
            # Psi1 = T*DPhi0/2/pi
            #Psi1 = self.T*self.DPhi0/(2*np.pi)
            if self.T != 0:
                #Psi1 = -(lambertw(-np.exp(-self.T)*self.T)-self.T)/self.T # Lambert W function calculates inverse of hyp2f1(1,1,2,-psi1) at self.z=0
                Psi1 = self.T
                #Psi2 = (2*window.laserI0*self.T*self.DPhi0/(2*np.pi))**(1/2)
            else:
                Psi1 = 0
            
            Psi2 = self.T*8 # 8 is a value out of the hat to obtain similar extreme point (valley or peak)
            
            match model:
                case "2PA":
                    psi1 = Psi1/(1+self.z**2/self.z0**2)
                    # self.Tz = hyp2f1(1/n,1/n,(n+1)/n,-psi**n)
                    self.Tz = hyp2f1(1,1,2,-psi1)
                    pass
                case "3PA":
                    psi2 = Psi2/(1+self.z**2/self.z0**2)
                    # self.Tz = hyp2f1(1/n,1/n,(n+1)/n,-psi**n)
                    self.Tz = hyp2f1(1/2,1/2,3/2,-(psi2)**2)
                
                case "2PA+3PA":
                    if Psi1 !=0:
                        psi1 = Psi1/(1+self.z**2/self.z0**2)
                        psi2 = Psi2/(1+self.z**2/self.z0**2)
                        f_x_psi1_psi2 = 1+psi1*(0.339*np.sin(0.498*psi2)-0.029)/(1+0.966*psi1*psi2**-0.718) # coupling function formula from OPTICS EXPRESS Vol. 13, No. 23 9231
                        self.Tz = hyp2f1(1,1,2,-psi1)*hyp2f1(1/2,1/2,3/2,-(psi2)**2)*f_x_psi1_psi2
                    
                    else:
                        self.Tz = np.ones_like(self.z)
                
                case "RSA":
                    self.Tz = np.ones_like(self.z)
                    
                case "SA":
                    self.Tz = np.ones_like(self.z)
                    
                case "2PA+SA":
                    self.Tz = np.ones_like(self.z)
                    
                case _:
                    self.Tz = np.ones_like(self.z)
            
            self.Tznorm = 2*self.Tz/(np.average(self.Tz[0:10])+np.average(self.Tz[len(self.Tz)-10:]))
        
        else:
            self.Tznorm = np.ones_like(self.z)

        return self.Tznorm

    def calculate_fm(self):
        '''3rd method called. Called by bigproduct()'''
        self.result = [(1j*self.Dphi0)**m/factorial(m)*self.product[m] for m in range(0,self.mm)]
        return self.result

    def bigproduct(self):
        '''2nd method called. Called by derive() for stype = "CA".'''
        # Big product operator (j from 1 to m)
        # if self.beta == 0:
        # self.product = np.ones(self.mm) # The result is already known. Uses numpy array for calculate_fm function to work properly

        #else: # THIS MUST BE CALCULATED TO TAKE INTO ACCOUNT TRANSMITTANCE THROUGH THE APERTURE
        self.product = []
        for m in range(0,self.mm):
            if m==0:
                self.product.append(1)
                continue # 0th value of m gives product=1 and in cumulative product it gives no contribution, so go to next iteration straight away.
            else:
                self.product.append(np.cumprod([1+1j*(j-1/2)/2/np.pi*self.T for j in range(1,m+1)])[-1]) # calculate cumulative product for given m
        #                                                                                                                          # and make 'product' contain all m products
        self.fm = self.calculate_fm()

        Tzo = self.open() # Returns final result for OA
        Tzc = self.closed() # Returns final result for CA
        
        self.Tznorm = Tzc/Tzo # This way, transmittance through aperture is taken into account
    
    def open(self):
        '''4th method called. Called by bigproduct()'''
        self.dr = 3*self.wa/self.ir # integration over radius step size
        self.open_sum = self.bigsum()
        return self.open_sum

    def closed(self):
        '''6th method called. Called by bigproduct()'''
        self.dr = self.ra/self.ir # integration over radius step size
        self.closed_sum = self.bigsum()
        return self.closed_sum
    
    def bigsum(self):
        '''5th and 7th method called. Called by open() and closed()'''
        # Big sum operator (m from 0 to "infinity")
        # integration over radius
        
        # Transmitted power
        self.Tz = 0

        for rr in range(self.ir):
            self.E = 0
            
            for m in range(0,self.mm):
                self.wm0 = self.wz/np.sqrt((2*m+1))
                self.dm = 0.5*self.k*self.wm0**2
                self.wm = self.wm0*np.sqrt(self.g**2+self.d**2/self.dm**2)
                self.tm = np.arctan(self.g/(self.d/self.dm))
                self.Rm = self.d/(1-self.g/(self.g**2+self.d**2/self.dm**2))

                self.E += self.fm[m]*np.exp(1j*self.tm)*self.wm0/self.wm/self.wz*np.exp((-1/self.wm**2+1j*np.pi/self.lda/self.Rm)*(rr*self.dr)**2)
        
            self.Tz += np.abs(self.E)**2*rr*self.dr # transmittance through aperture plane
        
        # T(z)=P_T/(S*P_i), where S=1-exp(-2*ra^2/wa^2)
        #S = 1-np.exp(-2*self.ra**2/self.wa**2)
        #I0 - instantaneous laser fluence - this is missing (we want to get it)
        #self.Tznorm = 3E8*8.854E-12*self.Tz/(S*self.w0**2/2*I0)
        
        # Normalize transmitted power before dividing CA/OA (the only way of normalization when we don't have I0)
        self.Tznorm = 2*self.Tz/(np.average(self.Tz[0:10])+np.average(self.Tz[len(self.Tz)-10:]))

        return self.Tznorm
class Fitting():
    def __init__(self, sample_type: Integration, amplitude, beamwaist, zero_level, centerpoint, nop, data):
        super(Fitting, self).__init__()
        self.sample_type = sample_type
        self.amplitude = amplitude
        self.beamwaist = beamwaist
        self.zero_level = zero_level
        self.centerpoint = centerpoint
        self.nop = nop
        self.ydata = data

    # The parametrized function to be plotted. It is also initial guess for automatic fitting.
    def manual(self, zero_level, centerpoint, amplitude, beamwaist, z_range, stype="CA") -> list:
        '''Returns integrated field for given input parameters
        
        ONLY FOR n2 FOR NOW!!!!!!!!!!!!!'''
        self.z_range = z_range # in meters
        self.sample_type.z = np.array([self.z_range*(zz - centerpoint)/self.nop-self.z_range/2 for zz in range(self.nop)]) # in meters
        window.get_general_parameters()
        if stype == "CA":
            self.sample_type.derive(amplitude,beamwaist,window.d0,window.ra,stype)
            cas = self.sample_type.closed_sum # Tznorm
            oas = self.sample_type.open_sum   # Tznorm
            result = (cas/oas)+(zero_level-1)
        elif stype == "OA":
            self.sample_type.derive(amplitude,beamwaist,window.d0,window.ra,stype)
            oas = self.sample_type.Tznorm
            result = oas+(zero_level-1)
        
        return result

    # The function to be minimized in automated fitting
    def fcn2min(self,params,*weights):
        pars = params.valuesdict().values()
        ynew = self.manual(*pars)
        return [w*((yn-yi)**2) for w,yn,yi in zip(weights,ynew,self.ydata)] # SSE

    # The actual processor for automatic fitting
    def automatic(self, z_range, ftype:str, stype:str, line_xydata):
        self.z_range = z_range

        if (ftype == "Solvent" and window.solventCA_customBeamwaist_checkBox.isChecked() is False):#\
            #or (ftype == "Sample" and window.sampleCA_customBeamwaist_checkBox.isChecked() is False):
            vary_beamwaist = False
        else:
            vary_beamwaist = True
        
        vary_centerpoint = True
        if ftype == "Solvent" and stype == "OA":
            if window.solventOA_customCenterPoint_checkBox.isChecked() is False:
                vary_centerpoint = False
            else:
                vary_centerpoint = True
        
        # Apply initial values
        self.params = Parameters()
        self.params.add('Zero', value=self.zero_level, min=0.75, max=1.25)
        self.params.add('Center', value=self.centerpoint, min=-50, max=50, vary=vary_centerpoint) # in number of datapoints
        
        if stype == "CA":
            self.params.add('DPhi0', value=self.amplitude, min=-2, max=2)
            self.params.add('Beamwaist', value=self.beamwaist, min=15E-6, max=150E-6, vary=vary_beamwaist)
            self.params.add('Zrange', value=self.z_range, vary=False)
        
        elif stype == "OA":
            self.params.add('T', value=self.amplitude, min=-2, max=2)
            self.params.add('Beamwaist', value=self.beamwaist, vary=False)
            self.params.add('Zrange', value=self.z_range, vary=False)
        
        if ftype == "Silica":
            #line = window.silica_figure.axes.get_lines()[0]
            #xs, ys = line.get_data()
            xs, ys = line_xydata
            weights = np.ones(np.shape(xs))
            if hasattr(window, 'silica_cursorPositions'):
                if len(window.silica_cursorPositions) == 2:
                    # CA ranges for weighting the fit
                    x1 = window.silica_cursorPositions[0][0]
                    x1_index = list(xs).index(min(xs, key=lambda x:abs(x-x1)))
                    # y1 = window.silica_cursorPositions[0][1]
                    #y1_index = list(ys).index(min(ys, key=lambda x:abs(x-y1)))
                    x2 = window.silica_cursorPositions[1][0]
                    x2_index = list(xs).index(min(xs, key=lambda x:abs(x-x2)))
                    # y2 = window.silica_cursorPositions[1][1]
                    #y2_index = list(ys).index(min(ys, key=lambda x:abs(x-y2)))
                    
                    x_sm, x_lg = sorted([x1_index, x2_index])
                    weights = [0 if (xi < x_sm or xi > x_lg) else 1 for xi in range(len(xs))]
            
            fitter = Minimizer(self.fcn2min,self.params,fcn_args=(weights))

        elif ftype == "Solvent":
            #line = window.solventCA_figure.axes.get_lines()[0]
            #xs, ys = line.get_data()
            xs, ys = line_xydata
            weights = np.ones(np.shape(xs))
            if stype == "CA":
                attribute_name = 'solventCA_cursorPositions'
            elif stype == "OA":
                attribute_name = 'solventOA_cursorPositions'
            
            if hasattr(window, attribute_name):
                if stype == "CA":
                    attribute = window.solventCA_cursorPositions
                elif stype == "OA":
                    attribute = window.solventOA_cursorPositions

                if len(attribute) == 2:
                    # CA ranges for weighting the fit
                    x1 = attribute[0][0]
                    x1_index = list(xs).index(min(xs, key=lambda x:abs(x-x1)))
                    # y1 = attribute[0][1]
                    # y1_index = list(ys).index(min(ys, key=lambda x:abs(x-y1)))
                    x2 = attribute[1][0]
                    x2_index = list(xs).index(min(xs, key=lambda x:abs(x-x2)))
                    # y2 = attribute[1][1]
                    # y2_index = list(ys).index(min(ys, key=lambda x:abs(x-y2)))
                    
                    x_sm, x_lg = sorted([x1_index, x2_index])
                    weights = [0 if (xi < x_sm or xi > x_lg) else 1 for xi in range(len(xs))]
            
            fitter = Minimizer(self.fcn2min,self.params,fcn_args=(weights))
        
        result = fitter.minimize(method='least_squares', max_nfev=1000)

        res_pars = result.params.valuesdict().values()
        result_line = self.manual(*res_pars, stype) # this will have to catch two lines (n2 and OA)

        result.params.pretty_print()

        return result, result_line

class MotorPositioner(QObject):
    def movetostart(self, progress_callback):
        start_pos = window.startPos_doubleSpinBox.value()
        end_pos = window.endPos_doubleSpinBox.value()
        if window.motor.position <= (start_pos+end_pos)/2:
            window.motor.move_to(start_pos)
            window.where_to_start = "start"
        else:
            window.motor.move_to(end_pos)
            window.where_to_start = "end"
        
        if window.motor.is_in_motion is True:
            print("Moving to starting position")
        
        while window.motor.is_in_motion is True:
            continue

        self.run(progress_callback)
    
    def movetocustompos(self, *args, **kwargs):
        text_val = window.custom_pos_dialog.new_pos.text()
        target = float(text_val.replace(",","."))

        if target > 100:
            return

        window.motor.move_to(target)
        while window.motor.is_in_motion is True:
            continue

        window.current_pos_chooser.setEnabled(True)

    def moveby(self, *args, **kwargs):
        move_step = (window.endPos_doubleSpinBox.value()-window.startPos_doubleSpinBox.value())/window.stepsScan_spinBox.value()
        if window.where_to_start == "end":
            window.motor.move_by(-move_step,blocking=True)
        else:
            window.motor.move_by(move_step,blocking=True)
        
        #while window.motor.is_in_motion is True:
        #    continue
    
    def movehome(self, *args, **kwargs):
        if window.motor.has_homing_been_completed is False:
            window.motor.move_home()
         
            time.sleep(0.2) # otherwise it may not move_home()
            if window.motor.is_in_motion is True:
                print("Homing now")
            while window.motor.has_homing_been_completed is False:
                continue # wait until homing is completed
            time.sleep(0.2) # wait a little more (so the motor.position gets exactly "0" position)
        
        return "Homing performed!"

    def run(self, progress_callback):
        window.data_acquisition_complete = False
        window.data_reversed = False # when backwards scan is performed, it later gets reversed (the data_reverse() method)

        time.sleep(0.2) # Sometimes the first datapoint is collected before the motor has settled

        nos = window.stepsScan_spinBox.value()

        for step in range(nos+1):
            if window.experiment_stopped is True:
                window.experiment_stopped = False
                window.running = False
                break
            self.step = step
            ### Create a task
            with nidaqmx.Task() as task:
                # Create MultiChannel "channel"
                task.ai_channels.add_ai_voltage_chan(window.detector_core_name+f"0:{window.number_of_channels_used}") # "Dev1/ai0:3"
                
                # Start Digital Edge
                task.triggers.start_trigger.dig_edge_src = "/Dev1/PFI0"
                task_trigger_src = task.triggers.start_trigger.dig_edge_src

                # Sample Clock
                rate0 = 1000 # 1 kHz (repetition rate of the laser)
                task.timing.cfg_samp_clk_timing(rate0,source=task_trigger_src,active_edge=Edge.FALLING, sample_mode=AcquisitionType.FINITE, samps_per_chan=window.samplesStep_spinBox.value())
                
                ### Data acquisition
                reader = nidaqmx.stream_readers.AnalogMultiChannelReader(task.in_stream)
                values_read = np.zeros((window.number_of_channels_used+1,window.samplesStep_spinBox.value()),dtype=np.float64)
                # Start Task
                task.start()
                
                time.sleep(0.2) # THIS IS THE TIME NEEDED TO COLLECT ALL SAMPLES

                # Acquire data
                window.data["positions"].append(window.motor.position)
                reader.read_many_sample(values_read, number_of_samples_per_channel=window.samplesStep_spinBox.value())
                
                # Take mean for each channel and append to "data" dictionary
                data_mean = np.mean(values_read,axis=1)

                for chan_no in range(window.number_of_channels_used):
                    window.data["absolute"][chan_no].append(data_mean[chan_no])
                # These loops need to be separated because "absolute" list at current step MUST BE defined at current step when accessed by "relative" list
                for chan_no in range(window.number_of_channels_used):
                    # Take last/current value (at index [step]) in "absolute" values list and divide by channel [1] (reference)
                    # And append to "relative" values list
                    window.data["relative"][chan_no].append(data_mean[chan_no]/window.data["absolute"][1][step])
                
                # 2) display data
                
                for type in window.charts.keys():
                    for chan_no in range(window.number_of_channels_used):
                        window.measurement_lines[type][chan_no].set_xdata(window.data["positions"]) # and set data of lines in "lines" dictionary to empty lists of positions and values
                        window.measurement_lines[type][chan_no].set_ydata(window.data[type][chan_no])

                y = window.data["absolute"][1]
                window.rms_value = np.abs(np.sqrt(np.mean([yi**2 for yi in y])) - y[0])/y[0]
                window.rms_text.set_text(f"RMS noise = {window.rms_value*100:.3f}%")
                
                window.measurement_plot_rescale(window.focusAt_comboBox.currentText())

                progress_callback.emit(step)

                # 3) move the motor
                if step == nos: # prevent useless additional step
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
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    #-db DATABASE -u USERNAME -p PASSWORD -size 20000
    parser.add_argument("-s",
                        "--settings",
                        dest = "settings_ini",
                        default = f"{os.path.join(os.path.dirname(__file__), 'settings.ini')}",
                        help="Path to settings file (with the filename.*ini)",
                        type=str)
    args = parser.parse_args()
    
    app = QtWidgets.QApplication(sys.argv)
    window = Window(args.settings_ini)
    app.setStyle("Fusion")
    window.default_palette = QtGui.QGuiApplication.palette()
    window.changeSkinDark() # Make sure the additional changes are applied
    app.exec_()
