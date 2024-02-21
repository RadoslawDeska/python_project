from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtCore import QSettings,QFile
from os import chmod, access, W_OK
from os.path import abspath, dirname, join, pardir

__all__ = [
    'DEFAULT_SETTINGS_LINES',
    'apply_settings',
    'save_settings',
    'save_as',
    'sort_settings',
    'write_settings_string',
    'ensure_default_settings_file'
    ]

settings_lines = [
                "[UI]",
                "ui_path=./window.ui",
                "defaults_location=./default_settings.ini",# expected location of default settings
                "[MeasurementTab]",
                "starting_position=35",
                "ending_position=75",
                "step_per_cycle=200",
                "samples_per_position=200",
                "[SavingTab]",
                "silica_thickness=4",
                "concentration=0",
                "wavelength=800",
                "main_directory=C:/z-scan/_wyniki",
                "[FittingTab]",
                "data_directory=C:/z-scan/_wyniki",
                "silica_thickness=4",
                "wavelength=800",
                "zrange=40",
                "aperture_diameter=1",
                "distance_from_focus_to_CA=260",
                "solvents_path=./solvents.json",# expected location of solvents file
                "previous_solvent=0",# index of solvent previously used
                "[Hardware]",
                "nidaqmx_device_name=/Dev1",
                "nidaqmx_no_of_channels=3",
                "nidaqmx_dig_edge_src=/Dev1/PFI0",
                "thorlabs_motor_id=40180184",
                "laser_repetition_rate=1000", # Hz
                "[Perks]",
                "end_beep_duration_ms=500", # ms
                "end_beep_tone_frequency=800" # Hz
                ]
DEFAULT_SETTINGS_LINES = '\n'.join(settings_lines)
            
def apply_settings(w, settings):
    # Measurement Tab
    w.startPos_doubleSpinBox.setValue(float(settings.value('MeasurementTab/starting_position')))
    w.endPos_doubleSpinBox.setValue(float(settings.value('MeasurementTab/ending_position')))
    w.stepsScan_spinBox.setValue(int(settings.value('MeasurementTab/step_per_cycle')))
    w.samplesStep_spinBox.setValue(int(settings.value('MeasurementTab/samples_per_position')))
    
    # Data saving Tab
    w.silicaThickness_dataSavingTab_doubleSpinBox.setValue(float(settings.value('SavingTab/silica_thickness')))
    w.concentration_dataSavingTab_doubleSpinBox.setValue(float(settings.value('SavingTab/concentration')))
    w.wavelength_dataSavingTab_doubleSpinBox.setValue(float(settings.value('SavingTab/wavelength')))
    w.mainDirectory_lineEdit.setText(settings.value('SavingTab/main_directory').replace("/","\\"))
    
    # Data fitting Tab
    w.dataDirectory_lineEdit.setText(settings.value('FittingTab/data_directory').replace("/","\\"))
    w.silicaThickness_dataFittingTab_doubleSpinBox.setValue(float(settings.value('FittingTab/silica_thickness')))
    w.wavelength_dataFittingTab_doubleSpinBox.setValue(float(settings.value('FittingTab/wavelength')))
    w.zscanRange_doubleSpinBox.setValue(float(settings.value('FittingTab/zrange')))
    w.apertureDiameter_doubleSpinBox.setValue(float(settings.value('FittingTab/aperture_diameter')))
    w.apertureToFocusDistance_doubleSpinBox.setValue(float(settings.value('FittingTab/distance_from_focus_to_CA')))
    w.solventName_comboBox.setCurrentIndex(int(settings.value('FittingTab/previous_solvent')))
    
    # w.dataDirectory_lineEdit.setText(settings.value('FittingTab/solvents_path').replace("/","\\"))

def save_settings(w, settings):
    # Measurement Tab
    settings.setValue('MeasurementTab/starting_position',w.startPos_doubleSpinBox.value())
    settings.setValue('MeasurementTab/ending_position', w.endPos_doubleSpinBox.value())
    settings.setValue('MeasurementTab/step_per_cycle', w.stepsScan_spinBox.value())
    settings.setValue('MeasurementTab/samples_per_position', w.samplesStep_spinBox.value())

    # Data saving Tab
    settings.setValue('SavingTab/silica_thickness', w.silicaThickness_dataSavingTab_doubleSpinBox.value())
    settings.setValue('SavingTab/concentration', w.concentration_dataSavingTab_doubleSpinBox.value())
    settings.setValue('SavingTab/wavelength', w.wavelength_dataSavingTab_doubleSpinBox.value())
    settings.setValue('SavingTab/main_directory', w.mainDirectory_lineEdit.text())

    # Data fitting Tab
    settings.setValue('FittingTab/data_directory', w.dataDirectory_lineEdit.text())
    settings.setValue('FittingTab/silica_thickness', w.silicaThickness_dataFittingTab_doubleSpinBox.value())
    settings.setValue('FittingTab/wavelength', w.wavelength_dataFittingTab_doubleSpinBox.value())
    settings.setValue('FittingTab/zrange', w.zscanRange_doubleSpinBox.value())
    settings.setValue('FittingTab/aperture_diameter', w.apertureDiameter_doubleSpinBox.value())
    settings.setValue('FittingTab/distance_from_focus_to_CA', w.apertureToFocusDistance_doubleSpinBox.value())
    settings.setValue('FittingTab/previous_solvent', w.solventName_comboBox.currentIndex())
    
def save_as(w, settings):
    path, _ = QFileDialog.getSaveFileName(w, "Save settings", join(abspath(dirname(__file__)),pardir), filter="INI file (*.ini)")
    
    save_settings(w,settings)
    result = sort_settings(settings)
    
    if path:
        write_settings_string(path,result)
    
        w.settings = QSettings(path, QSettings.IniFormat)

def sort_settings(settings:QSettings) -> str:
    """Sorts settings by group name in ascending order
    and sorts settings within the group by keyname in ascending order

    Args:
        settings (QSettings): 

    Returns:
        str: Sorted settings items for writing to a settings file
    """    
    result = ""
    groups = {}
    for key in settings.allKeys():
        value = settings.value(key)
        group, _, subkey = key.partition('/')
        if group not in groups:
            groups[group] = {}
        groups[group][subkey] = value

    for group, settings in sorted(groups.items()):
        result += f"[{group}]\n"
        for key, value in sorted(settings.items()):
            result += f"{key}={value}\n"
        result += "\n"
    
    return result

def write_settings_string(path, s:str):
    with open(path, mode="w", encoding="utf-8") as f:
        f.write(s)

def ensure_default_settings_file():
    '''Make sure default settings file exists and is not writeable by the user'''
    if not QFile(join(dirname(__file__),pardir, "default_settings.ini")).exists():
        # create default settings file
        with open(join(dirname(__file__),pardir, "default_settings.ini"), mode="w", encoding="utf-8") as fi:
            fi.write(DEFAULT_SETTINGS_LINES)
        
        temp_settings = QSettings(join(dirname(__file__), pardir, "default_settings.ini"), QSettings.IniFormat)
        write_settings_string(join(dirname(__file__), pardir, "default_settings.ini"),sort_settings(temp_settings))
    
    # prevent editing by setting permissions to read-only
    if access(join(dirname(__file__), pardir, "default_settings.ini"), W_OK):
        chmod(join(dirname(__file__), pardir, "default_settings.ini"), 0o444)