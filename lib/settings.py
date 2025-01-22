import os
from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtCore import QSettings,QFile,QFileInfo
from os import chmod, access, W_OK
from os.path import join

__all__ = [
    'DEFAULT_SETTINGS_LINES',
    'apply_settings',
    'load_settings',
    'save_settings',
    'save_as',
    'sort_settings',
    'write_settings_string',
    'ensure_default_settings_file'
    ]

settings_lines = [
                "[MeasurementTab]",
                "starting_position=35",
                "ending_position=75",
                "steps_per_scan=200",
                "samples_per_position=200",
                "number_of_scans=1",
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
                "solvents_path=./solvents.json",  # expected location of solvents file
                "previous_solvent=0",  # index of solvent previously used
                "[Hardware]",
                "nidaqmx_device_name=/Dev1",
                "nidaqmx_no_of_channels=3",
                "nidaqmx_dig_edge_src=/Dev1/PFI0",
                "thorlabs_motor_id=40180184",
                "laser_repetition_rate=1000",  # Hz
                "[Perks]",
                "end_beep_duration_ms=500",  # ms
                "end_beep_tone_frequency=800",  # Hz
                "end_beep_emit=True",  # whether the beep should be emitted or not
                "[UI]",
                "ui_path=./window.ui",
                "defaults_location=./default_settings.ini"  # expected location of default settings
                ]
DEFAULT_SETTINGS_LINES = '\n'.join(settings_lines)
            
def apply_settings(w, settings):
    # Measurement Tab
    w.startPos_doubleSpinBox.setValue(float(settings.value('MeasurementTab/starting_position')))
    w.endPos_doubleSpinBox.setValue(float(settings.value('MeasurementTab/ending_position')))
    w.stepsScan_spinBox.setValue(int(settings.value('MeasurementTab/steps_per_scan')))
    w.samplesStep_spinBox.setValue(int(settings.value('MeasurementTab/samples_per_position')))
    w.numberOfScans_spinBox.setValue(int(settings.value('MeasurementTab/number_of_scans')))
    
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
    w.solventName_comboBox.setCurrentIndex(min(0, int(settings.value('FittingTab/previous_solvent'))))  # Fallbacks to 0 if value from settings is invalid
    # w.dataDirectory_lineEdit.setText(settings.value('FittingTab/solvents_path').replace("/","\\"))

def load_settings(w=None,user=False):
    temp_path = join(os.getcwd(),'settings.ini')
    temp_settings = QSettings(temp_path, QSettings.IniFormat)
    # this holds the value of last path of the settings file
    
    if w is None:  # for the purpose of running the program with previous settings file
        # this part is supposed to be run ONLY before `Window` instance is created
        # print("Making sure default settings file exists")
        ensure_default_settings_file()
        
        # print("Checking temporary settings for last path...")
        last_path = temp_settings.value('UI/last_settings_location')
        if last_path is not None:
            print(f"Last remembered path is: {last_path}")
            info = QFileInfo(last_path)
            if os.path.exists(last_path):
                if bool(info.permissions() & QFile.WriteUser):
                    # print("The file exists. Continuing to start the program.")
                    return last_path
                else:
                    try:
                        chmod(last_path, 0o777)
                        return last_path
                    except PermissionError:
                        print(f"Cannot write to the file in: {last_path}. Reverting to {temp_path}")
        
        info = QFileInfo(temp_path)
        if os.path.exists(temp_path):
            if bool(info.permissions() & QFile.WriteUser):
                print(f"No path found/File not found. Found file: {temp_path}.",
                    "Continuing to start the program with the path of temporary settings file.",sep="\n")
                temp_settings.setValue('UI/last_settings_location',temp_path)
            else:
                try:
                    chmod(temp_path, 0o777)
                except PermissionError:
                    print(f"Cannot write to the file in: {temp_path}. Ensure you have write access to the program directory.")
                    return
        else:
            # print(f"No path found/File not found. Creating defaults in: {temp_path}.")
            write_settings_string(temp_path, DEFAULT_SETTINGS_LINES+f"\nlast_settings_location={temp_path}")
        
        # AT THIS POINT WHATEVER THE PROGRAM STARTED WITH, THERE IS AT LEAST THE DEFAULT SETTINGS IN THE (default_)settings.ini FILES.
        # APART FROM NO WRITE ACCESS TO THE PROGRAM DIRECTORY. THIS HAS TO BE SOLVED BY THE USER - THE PROGRAM WON'T START.
        return temp_path
    
    if user:  # for the purpose of loading the settings with user-selected file
        # print("User selects the settings file")
        path, _ = QFileDialog.getOpenFileName(w, "Load Settings", os.getcwd(), filter="INI file (*.ini)")
        
        if path:
            # print(f"User selected the settings file: {path}")
            temp_settings.setValue('UI/last_settings_location',path)
            temp_settings.sync()
            # print("Settings written to temp_settings.")
            w.settings = QSettings(path, QSettings.IniFormat)
            apply_settings(w,w.settings)
            # print(f"New settings applied from file: {path}")

def save_settings(w, settings):
    # Measurement Tab
    settings.setValue('MeasurementTab/starting_position',w.startPos_doubleSpinBox.value())
    settings.setValue('MeasurementTab/ending_position', w.endPos_doubleSpinBox.value())
    settings.setValue('MeasurementTab/steps_per_scan', w.stepsScan_spinBox.value())
    settings.setValue('MeasurementTab/samples_per_position', w.samplesStep_spinBox.value())
    settings.setValue('MeasurementTab/number_of_scans', w.numberOfScans_spinBox.value())

    # Data saving Tab
    settings.setValue('SavingTab/silica_thickness', w.silicaThickness_dataSavingTab_doubleSpinBox.value())
    # settings.setValue('SavingTab/concentration', w.concentration_dataSavingTab_doubleSpinBox.value())
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
    path, _ = QFileDialog.getSaveFileName(w, "Save settings", os.getcwd(), filter="INI file (*.ini)")
    
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

def ensure_default_settings_file(default=True):
    '''Make sure default settings file exists and is not writeable by the user'''
    if default:
        if not QFile("default_settings.ini").exists():
            # create default settings file
            print("Defaults don't exist. Creating default settings file.")
            write_settings_string("default_settings.ini",DEFAULT_SETTINGS_LINES+"\nlast_settings_location=./default_settings.ini")
            
            temp_settings = QSettings("default_settings.ini", QSettings.IniFormat)
            write_settings_string("default_settings.ini",sort_settings(temp_settings))
        
        # prevent editing by setting permissions to read-only
        if access("default_settings.ini", W_OK):
            chmod("default_settings.ini", 0o444)
        
        print("Default setting file has been found or successfully created.")