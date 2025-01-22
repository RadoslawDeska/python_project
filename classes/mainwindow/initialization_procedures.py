import sys  # noqa: F401

from classes.custom_qwidgets.autocomplete import *
from classes.mainwindow.charts.chart_initialization import *
from classes.mainwindow.managers import *
from classes.mainwindow.motor_control import *
from classes.mainwindow.user_triggered_signals import *


def setup_triggers(self):
    window_menu_signals(self)
    measurement_tab_signals(self)
    saving_tab_signals(self)
    fitting_tab_signals(self)

def set_initial_states(self):
    # Measurement states
    self.initializing = False # variable for watching if initialize method has been called
    self.initialized = False
    self.clearing = False
    
    self.running = False # Variable for watching if run method has been called
    
    self.data_acquisition_complete = False
    
    # Fitting lines
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

def set_additional_variables(self):
    self.offset = 0  # margin control in measurement charts
    self.scanning_from = "start" # or "end" - whether the scan starts from start_pos or from end_pos
    self.current_run = 0
    
    self.last_solvent = None
    self.allsamples_db = {}  # dictionary consisting of name and properties of given specimen (silica/solvent/sample) (for database connections)
    self.user_samples = {}
    self.boxIndexChangedFlags = {}  # keeps reference to previous index values in QComboBox (its name stored as key, index as value)
    
    self.data = {"positions": [], "relative": {}, "absolute": {}}
    self.batch_data = {}  # batch data dictionary contains run number and data associated
    self.already_sorted = {}  # required for effective sorting of data for calculating the mean in data acquisition
    self.all_positions = []
    self.all_absolutes = {}
    self.all_relatives = {}
    
    self.previous_end_pos = self.endPos_doubleSpinBox.value()
    self.previous_start_pos = self.startPos_doubleSpinBox.value()
    self.previous_steps_scan = self.stepsScan_spinBox.value()

    self.header_correct = False
    self.rms_value = 0.0 # RMS of Reference signal

def modify_gui_to_initial_looks(self):
    _select_JSON_file(self, "Startup")  # load default solvents
    self.add_nonqt_widgets()
    self.solventName_dataSavingTab_comboBox.setCurrentIndex(-1)
    self.update_labels()
    self.set_logging_items((self.rawData_listWidget,self.fullExp_listWidget), self.numberOfScans_spinBox.value())
    self.toggle_absorption_model("solvent")
    self.toggle_saturation_model("solvent")
    
    initialize_measurement_charts(self)
    initialize_fitting_charts(self)

def initialize(self, *args, **kwargs):
    """Accessed by button click in GUI"""
    self.initializing = True
    self.initialize_pushButton.setEnabled(False)
    self.initLED_pushButton.setEnabled(True)
    self.toggle_user_input_lock(False)
    
    def configure_detectors(settings):
        device_name = settings.value("Hardware/nidaqmx_device_name")
        return device_name+"/ai", int(settings.value("Hardware/nidaqmx_no_of_channels"))
    
    self.detector_core_name, self.number_of_channels_used = configure_detectors(self.settings)
    configure_charts(self)
    
    # Initialize translation stage motor
    motor_detect(self)
    set_motor_home(self)  # this starts in separate thread

