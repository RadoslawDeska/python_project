from PyQt5.QtCore import QSettings
from PyQt5.QtWidgets import QComboBox

from classes.custom_qwidgets.autocomplete import *
from classes.custom_qwidgets.autocomplete import MyComboBox, MyComboBox_NonEditable
from classes.mainwindow.managers import SampleManager as SM
from classes.mainwindow.measurement_processing import update_fullLog
from classes.mainwindow.colorthemes import improperValueStyle
from lib.settings import apply_settings, load_settings, save_as, save_settings

__all__ = ["window_menu_signals", "measurement_tab_signals", "saving_tab_signals", "fitting_tab_signals"]


def window_menu_signals(self):
    def menu():
        # File
        self.actionExit.triggered.connect(self.close)
        # Samples
        self.actionLoadSamples.triggered.connect(lambda: _select_JSON_file(self, "actionLoadSamples"))
        self.actionSaveSamples.triggered.connect(lambda: _save_JSON_file(self, "actionSaveSamples"))
        # Solvents
        self.actionLoadSolvents.triggered.connect(lambda: _select_JSON_file(self, "actionLoadSolvents"))
        self.actionSaveSolvents.triggered.connect(lambda: _save_JSON_file(self, "actionSaveSolvents"))
        # Settings
        self.actionLoadSettings.triggered.connect(lambda: load_settings(self, user=True))  # noqa: F405
        self.actionSaveSettings.triggered.connect(lambda: save_settings(self, self.settings))  # noqa: F405
        self.actionSaveAsSettings.triggered.connect(lambda: save_as(self, self.settings))  # noqa: F405
        self.actionRestoreDefaultSettings.triggered.connect(
            lambda: apply_settings(self, QSettings(self.settings.value("UI/defaults_location"), QSettings.IniFormat))
        )  # noqa: F405
        # View
        self.actionLight.triggered.connect(self.changeSkinLight)
        self.actionDark.triggered.connect(self.changeSkinDark)

    def topbar():
        self.stopRed_pushButton.clicked.connect(self.stop_experiment)

    menu()
    topbar()


def measurement_tab_signals(self):
    def measurement_control(self):
        self.initialize_pushButton.clicked.connect(self.initialize)
        self.run_pushButton.clicked.connect(self.set_motor_to_start)
        self.clear_pushButton.clicked.connect(self.clear_measurement)

    def measurement_configuration(self):
        self.extremes_radioButton.toggled.connect(self.toggle_measurement_range_selector)
        self.centerRange_radioButton.toggled.connect(self.toggle_measurement_range_selector)

        self.endPos_doubleSpinBox.editingFinished.connect(lambda: self.rescale_measurement_plots(who_called="end"))
        self.startPos_doubleSpinBox.editingFinished.connect(lambda: self.rescale_measurement_plots(who_called="start"))
        self.focalPoint_doubleSpinBox.editingFinished.connect(lambda: self.rescale_measurement_plots(who_called="focal"))
        self.zscanRange_measurementTab_doubleSpinBox.editingFinished.connect(lambda: self.rescale_measurement_plots(who_called="range"))
        self.stepsScan_spinBox.editingFinished.connect(lambda: self.rescale_measurement_plots())
        self.numberOfScans_spinBox.editingFinished.connect(
            lambda: self.set_logging_items((self.rawData_listWidget, self.fullExp_listWidget), self.numberOfScans_spinBox.value())
        )
        self.numberOfScans_spinBox.editingFinished.connect(self.update_labels)

    def charts_display_control(self):
        self.showScans_comboBox.currentIndexChanged.connect(lambda: self.plot_updated_data(self.batch_data))
        self.focusAt_comboBox.currentIndexChanged.connect(lambda: self.rescale_measurement_plots(self.focusAt_comboBox.currentText()))

    measurement_control(self)
    measurement_configuration(self)
    charts_display_control(self)


def saving_tab_signals(self):
    def measurement_parameters():
        self.apply_pushButton.clicked.connect(lambda: update_indicator(self))
        
        self.cuvetteType_comboBox.currentTextChanged.connect(lambda current_cuvette: toggle_codeOfSample(self, current_cuvette))
        # codeOfSample POPUP
        self.codeOfSample_comboBox.popupShow.connect(lambda: on_popup_show(self.codeOfSample_comboBox))
        self.codeOfSample_comboBox.popupHide.connect(lambda: self.codeOfSample_comboBox.setEditable(True))
        self.codeOfSample_comboBox.view().focusOutSignal.connect(lambda box: on_focus_out(self, box))
        # codeOfSample ITEM SELECTION/DELETION
        self.codeOfSample_comboBox.textActivated.connect(lambda: on_focus_out(self, self.codeOfSample_comboBox))
        self.codeOfSample_comboBox.lineEdit().editingFinished.connect(lambda: on_focus_out(self, self.codeOfSample_comboBox))
        self.codeOfSample_comboBox.focusOutSignal.connect(lambda: on_focus_out(self, self.codeOfSample_comboBox))
        self.codeOfSample_comboBox.currentIndexChanged.connect(lambda: on_focus_out(self, self.codeOfSample_comboBox))
        self.codeOfSample_comboBox.itemDeleted.connect(lambda sample_code: on_item_delete(self, sample_code))
        
        self.concentration_dataSavingTab_doubleSpinBox.editingFinished.connect(
            lambda: _populate_combo_with_params(self,
                                                self.codeOfSample_comboBox,
                                                caller="concentration"))

        self.solventName_dataSavingTab_comboBox.popupShow.connect(lambda: on_popup_show(self.solventName_dataSavingTab_comboBox))

        # Synchronize solvent list in Data Fitting tab
        # self.solventName_dataSavingTab_comboBox.currentIndexChanged.connect(
        #     lambda: sync_combos(self, self.solventName_dataSavingTab_comboBox, self.solventName_comboBox)
        # )
        # self.solventName_dataSavingTab_comboBox.lineEdit().editingFinished.connect(
        #     lambda: sync_combos(self, self.solventName_dataSavingTab_comboBox, self.solventName_comboBox)
        # )

        self.update_pushButton.clicked.connect(lambda: update_fullLog(self, update_button_clicked=True))

    def filenames():
        self.chooseDirectory_pushButton.clicked.connect(lambda: self.choose_dir(caller="DataSaving"))

    def measurement_saving():
        self.saveData_pushButton.clicked.connect(self.data_save)
        self.sendToFit_pushButton.clicked.connect(
            lambda: self.data_loader(caller="Current Measurement", ftype=self.cuvetteType_comboBox.currentText())
        )

    def logging_tabs():
        self.rawData_listWidget.itemClicked.connect(self.switch_widget_items)
        self.fullExp_listWidget.itemClicked.connect(self.switch_widget_items)

    measurement_parameters()
    filenames()
    measurement_saving()
    logging_tabs()


def fitting_tab_signals(self):
    def files():
        self.customDirectory_pushButton.clicked.connect(lambda: self.choose_dir(caller="DataFitting"))
        self.customSilicaFile_pushButton.clicked.connect(lambda: self.data_loader(caller="Load From File", ftype="silica"))
        self.customSolventFile_pushButton.clicked.connect(lambda: self.data_loader(caller="Load From File", ftype="solvent"))
        self.customSampleFile_pushButton.clicked.connect(lambda: self.data_loader(caller="Load From File", ftype="sample"))

    def general_parameters():
        self.customApertureDiameter_checkBox.stateChanged.connect(lambda: self.enable_custom("ApertureDiameter"))
        self.customApertureToFocusDistance_checkBox.stateChanged.connect(lambda: self.enable_custom("ApertureDistance"))
        self.customSilicaThickness_checkBox.stateChanged.connect(lambda: self.enable_custom("SilicaThickness"))
        self.customWavelength_checkBox.stateChanged.connect(lambda: self.enable_custom("Wavelength"))

        self.zscanRange_doubleSpinBox.editingFinished.connect(self.set_new_positions)
        self.customZscanRange_checkBox.stateChanged.connect(lambda: self.enable_custom("ZscanRange"))

    def solvent_properties():
        self.solventName_comboBox.popupShow.connect(lambda: on_popup_show(self.solventName_comboBox))
        self.solventName_comboBox.popupHide.connect(lambda: self.solventName_comboBox.setEditable(True))
        self.solventName_comboBox.currentIndexChanged.connect(lambda: _load_item(self, self.solventName_comboBox))

        self.solventDensity_doubleSpinBox.editingFinished.connect(lambda: _populate_combo_with_params(self, self.solventName_comboBox))
        self.solventRefrIdx_doubleSpinBox.editingFinished.connect(lambda: _populate_combo_with_params(self, self.solventName_comboBox))

        self.solventDensity_doubleSpinBox.valueChanged.connect(lambda: improperValueStyle(self, self.solventDensity_doubleSpinBox))
        self.solventRefrIdx_doubleSpinBox.valueChanged.connect(lambda: improperValueStyle(self, self.solventRefrIdx_doubleSpinBox))

    def sample_properties():
        self.customConcentration_checkBox.stateChanged.connect(lambda: self.enable_custom("Concentration"))

    def silicaCA_tab():
        self.silica_fit_pushButton.clicked.connect(lambda: self.fit_automatically(ftype="silica", stype="CA"))
        self.silica_fixROI_checkBox.stateChanged.connect(lambda: self.enable_cursors(ftype="silica", stype="CA"))

        self.slider_fit_manually_connect(self.silica_RayleighLength_slider, "Connect")
        self.slider_fit_manually_connect(self.silica_centerPoint_slider, "Connect")
        self.slider_fit_manually_connect(self.silica_zeroLevel_slider, "Connect")
        self.slider_fit_manually_connect(self.silica_DPhi0_slider, "Connect")
        # self.silica_RayleighLength_slider.valueChanged.connect(lambda: self.fit_manually(ftype="silica", stype="CA"))
        # self.silica_centerPoint_slider.valueChanged.connect(lambda: self.fit_manually(ftype="silica", stype="CA"))
        # self.silica_zeroLevel_slider.valueChanged.connect(lambda: self.fit_manually(ftype="silica", stype="CA"))
        # self.silica_DPhi0_slider.valueChanged.connect(lambda: self.fit_manually(ftype="silica", stype="CA"))
        self.silica_filterSize_slider.valueChanged.connect(lambda: self.reduce_noise_in_data(self.silica_data_set, ftype="silica", stype="CA"))

    def solventCA_tab():
        self.solventCA_RayleighLength_slider.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="CA"))
        self.solventCA_centerPoint_slider.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="CA"))
        self.solventCA_zeroLevel_slider.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="CA"))
        self.solventCA_DPhi0_slider.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="CA"))
        self.solventCA_filterSize_slider.valueChanged.connect(lambda: self.reduce_noise_in_data(self.solvent_data_set, ftype="solvent", stype="CA"))

        self.solventCA_fit_pushButton.clicked.connect(lambda: self.fit_automatically(ftype="solvent", stype="CA"))
        self.solventCA_fixROI_checkBox.stateChanged.connect(lambda: self.enable_cursors(ftype="solvent", stype="CA"))
        self.solventCA_customBeamwaist_checkBox.stateChanged.connect(lambda: self.enable_custom("SolventBeamwaist"))

    def solventOA_tab():
        self.solventOA_centerPoint_slider.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="OA"))
        self.solventOA_zeroLevel_slider.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="OA"))
        self.solventOA_T_slider.valueChanged.connect(lambda: self.fit_manually(ftype="solvent", stype="OA"))
        self.solventOA_filterSize_slider.valueChanged.connect(lambda: self.reduce_noise_in_data(self.solvent_data_set, ftype="solvent", stype="OA"))

        self.solventOA_fit_pushButton.clicked.connect(lambda: self.fit_automatically(ftype="solvent", stype="OA"))
        self.solventOA_fixROI_checkBox.stateChanged.connect(lambda: self.enable_cursors(ftype="solvent", stype="OA"))
        self.solventOA_isAbsorption_checkBox.stateChanged.connect(lambda: self.toggle_absorption_model(ftype="solvent"))
        self.solventOA_absorptionModel_comboBox.currentIndexChanged.connect(lambda: self.toggle_saturation_model(ftype="solvent"))
        self.solventOA_customCenterPoint_checkBox.stateChanged.connect(lambda: self.enable_custom("SolventCenterPoint"))

    files()
    general_parameters()
    solvent_properties()
    sample_properties()
    silicaCA_tab()
    solventCA_tab()
    solventOA_tab()


# CUSTOM SLOTS
def check_sample_type(self):
    '''This function ensures that a sample code is matched with its corresponding cuvette type
    (silica, solvent, or sample). If the pair hasn't been previously saved, it records them;
    otherwise, it verifies that the current selections align with the saved pair and requests
    to change the cuvette type selection.'''
    cuvette = self.cuvetteType_comboBox.currentText()
    sample_code = self.codeOfSample_comboBox.currentText()
    
    print(f"{sample_code=}, {cuvette=}")
    # Prevent empty keys
    if not sample_code:
        return
    
    # Save non-existing pair
    elif sample_code not in SM().allsamples_db:
        print("Saving parameters")
        SM().register_sample(sample_code, cuvette, **SM().read_parameters(self))
    
    # Non-matching selections found, toggle cuvette type
    elif cuvette != SM().allsamples_db[sample_code][0]:
        print("Non-matching selections found, toggle cuvette type")
        toggle_cuvetteType(self, SM().allsamples_db[sample_code][0])
    
    else:
        print("all is well")  # debug print

def update_indicator(self):
    sample_code = self.codeOfSample_comboBox.currentText()
    cuvette = self.cuvetteType_comboBox.currentText().lower()
    
    if sample_code in SM().allsamples_db:
        params = SM().allsamples_db[sample_code][1]
        match cuvette:
            case "silica":
                self.thickness_indicator.setText(str(params["silica_thicknessMM"]))
                self.density_indicator.setText("-")
                self.index_indicator.setText("-")
                self.concentration_indicator.setText("-")
                self.solvent_indicator.setText("-")
            case "solvent":
                self.density_indicator.setText(str(params["solvent_density"]))
                self.index_indicator.setText(str(params["solvent_index"]))
                self.thickness_indicator.setText("-")
                self.concentration_indicator.setText("-")
                self.solvent_indicator.setText("-")
            case "sample":
                self.concentration_indicator.setText(str(params["sample_concentration"]))
                self.solvent_indicator.setText(str(params["sample_solvent"]))
                self.thickness_indicator.setText("-")
                self.density_indicator.setText("-")
                self.index_indicator.setText("-")
   
def toggle_codeOfSample(self, current_cuvette):
    """
    Adjusts the state of various UI elements based on the type of cuvette selected.
    Parameters:
    current_cuvette (str): The type of cuvette selected. Can be "silica", "solvent", or "sample".
    Behavior:
    - For "silica":
        - Sets the codeOfSample_comboBox to "silica" and disables it.
        - Makes silicaThickness_dataSavingTab_doubleSpinBox editable.
        - Sets concentration_dataSavingTab_doubleSpinBox to 0 and makes it read-only.
        - Disables solventName_dataSavingTab_comboBox and clears its selection.
    - For "solvent":
        - Sets the codeOfSample_comboBox to the last selected solvent (if any) and enables it.
        - Makes silicaThickness_dataSavingTab_doubleSpinBox read-only.
        - Sets concentration_dataSavingTab_doubleSpinBox to 0 and makes it read-only.
        - Disables solventName_dataSavingTab_comboBox and clears its selection.
    - For "sample":
        - Clears the selection of codeOfSample_comboBox and enables it.
        - Makes silicaThickness_dataSavingTab_doubleSpinBox read-only.
        - Makes concentration_dataSavingTab_doubleSpinBox editable.
        - Sets solventName_dataSavingTab_comboBox to the last selected solvent and enables it.
    - For any other value:
        - No action is taken.
    Additionally, the method calls changeReadOnlyStyle() to update the UI styles based on the read-only states.
    """
    
    match current_cuvette:
        case "silica":
            print("silica selected")
            self.codeOfSample_comboBox.setCurrentText("silica")
            self.codeOfSample_comboBox.setEnabled(False)
            self.silicaThickness_dataSavingTab_doubleSpinBox.setReadOnly(False)
            self.concentration_dataSavingTab_doubleSpinBox.setValue(0)
            self.concentration_dataSavingTab_doubleSpinBox.setReadOnly(True)
            self.solventName_dataSavingTab_comboBox.setEnabled(False)
            self.solventName_dataSavingTab_comboBox.setCurrentIndex(-1)
        case "solvent":
            print("solvent selected")
            self.codeOfSample_comboBox.setCurrentText(self.last_solvent if self.last_solvent else "")
            self.codeOfSample_comboBox.setEnabled(True)
            self.silicaThickness_dataSavingTab_doubleSpinBox.setReadOnly(True)
            self.concentration_dataSavingTab_doubleSpinBox.setValue(0)
            self.concentration_dataSavingTab_doubleSpinBox.setReadOnly(True)
            self.solventName_dataSavingTab_comboBox.setEnabled(False)
            self.solventName_dataSavingTab_comboBox.setCurrentIndex(-1)
        case "sample":
            print("1) sample selected")
            self.codeOfSample_comboBox.blockSignals(True)
            self.codeOfSample_comboBox.setCurrentIndex(-1)
            self.codeOfSample_comboBox.blockSignals(False)
            print("2) code of sample removed")
            self.codeOfSample_comboBox.setEnabled(True)
            print("3) enabled editing code of sample")
            self.silicaThickness_dataSavingTab_doubleSpinBox.setReadOnly(True)
            print("4) disabled editing silica thickness")
            self.concentration_dataSavingTab_doubleSpinBox.setReadOnly(False)
            print("5) enabled editing concentration")
            if self.last_solvent:
                self.solventName_dataSavingTab_comboBox.setCurrentText(self.last_solvent)
                print(f"6) solvent combo set to {self.last_solvent}")
            else:
                print(f"6) solvent combo NOT SET, because {self.last_solvent=}")
            self.solventName_dataSavingTab_comboBox.setEnabled(True)
            print("7) solvent combo enabled")

    self.changeReadOnlyStyle()

def toggle_cuvetteType(self, target_cuvette):
    '''Executed on demand by check_sample_type()'''
    match target_cuvette:
        case "silica":
            self.silicaThickness_dataSavingTab_doubleSpinBox.setReadOnly(False)
            self.concentration_dataSavingTab_doubleSpinBox.setValue(0)
            self.solventName_dataSavingTab_comboBox.setCurrentIndex(-1)
            self.codeOfSample_comboBox.setEnabled(False)
            self.concentration_dataSavingTab_doubleSpinBox.setReadOnly(True)
            self.cuvetteType_comboBox.blockSignals(True)  # This prevents self.toggle_codeOfSample() - otherwise sample code index is reset to -1
            self.cuvetteType_comboBox.setCurrentText("silica")  # This otherwise would trigger self.toggle_codeOfSample(current_cuvette="silica").
            self.cuvetteType_comboBox.blockSignals(False)  # reset the signal management
        
        case "solvent":
            response = self.showdialog("Question", "Did you mean SOLVENT cuvette type?")
            # self.codeOfSample_comboBox.setCurrentText(current_code)
            if response:
                self.last_solvent = self.codeOfSample_comboBox.currentText()
                self.silicaThickness_dataSavingTab_doubleSpinBox.setReadOnly(True)
                self.cuvetteType_comboBox.blockSignals(True)  # This prevents self.toggle_codeOfSample() - otherwise sample code index is reset to -1
                self.cuvetteType_comboBox.setCurrentText("solvent")  # This otherwise would trigger self.toggle_codeOfSample(current_cuvette="solvent").
                # The order matters because changing cuvetteType uses last_solvent,
                # so last_solvent has to be already updated to newly entered value.
                self.codeOfSample_comboBox.setEnabled(True)  # update manually
                self.concentration_dataSavingTab_doubleSpinBox.setValue(0)
                self.concentration_dataSavingTab_doubleSpinBox.setReadOnly(True)  # update manually
                self.solventName_dataSavingTab_comboBox.setEnabled(False)  # update manually

                self.cuvetteType_comboBox.blockSignals(False)  # reset the signal management
            
            else:
                SM().register_sample(self.codeOfSample_comboBox.currentText(), self.cuvetteType_comboBox.currentText(), **SM().read_parameters(self))
        
        case "sample":
            response = self.showdialog("Question", "Did you mean SAMPLE cuvette type?")
            if response:
                current_code = self.codeOfSample_comboBox.currentText()
                sample_params = SM().allsamples_db[current_code][0]
                self.silicaThickness_dataSavingTab_doubleSpinBox.setReadOnly(True)
                self.concentration_dataSavingTab_doubleSpinBox.setValue(sample_params["sample_concentration"])
                self.solventName_dataSavingTab_comboBox.setCurrentText(sample_params["sample_solvent"])
            else:
                SM().register_sample(self.codeOfSample_comboBox.currentText(), self.cuvetteType_comboBox.currentText(), **SM().read_parameters(self))

    self.changeReadOnlyStyle()

    _load_item(self, self.codeOfSample_comboBox)

def indexChangedFlag(self, box: QComboBox, box_name: str):
    """Check if indexChanged signal is triggered for QComboBox with variable name: `box_name`."""
    if box_name in self.boxIndexChangedFlags:
        if self.boxIndexChangedFlags[box_name] != box.currentIndex():
            # remember the current value
            self.boxIndexChangedFlags[box_name] = box.currentIndex()
            # inform that index was changed
            return True
        else:
            # inform that index wasn't changed
            return False
    else:  # remember the current value
        self.boxIndexChangedFlags[box_name] = box.currentIndex()
        # inform that index was changed, because even if it didn't,
        # the functionality should lead to bringing up parameters related to current QComboBox values
        return True

# Define the on_popup_show function
def on_popup_show(box: MyComboBox):
    all_items = box.allItems()
    
    print("Popup shown. Items:",", ".join(item for item in all_items))

    # Fix disappearing elements of solventName elements
    if isinstance(box, MyComboBox_NonEditable):
        box.setEditable(False)
        box.setStyleSheet("")
        return
    
    # Ensure that the combobox doesn't lock itself
    if all_items:
        box.setEditable(False)
    else:
        box.setEditable(True)  # in non-editable state when there is no item in the list


def on_popup_close(box: MyComboBox_NonEditable):
    box.setEditable(False)
    
# Define focus out actions
def on_focus_out(self, box):
    if not box:
        return
    print("Focus is out... ",end="")
    if box == self.codeOfSample_comboBox:
        print("of self.codeOfSample_comboBox")
        if box.focus_out_processed:
            # Already processed
            return
        box.focus_out_processed = True  # reset the flag, so that the FocusOut event is not re-processed
    
        if len(box.allItems()) == 1 or box.currentText().lower() == "silica":
            toggle_cuvetteType(self, "silica")  # This is the last element and it is of type "silica"
            # return
    
    elif box == self.codeOfSample_comboBox.view():
        print("of self.codeOfSample_comboBox.view()")
        if len(self.codeOfSample_comboBox.allItems()) == 1 or self.codeOfSample_comboBox.currentText().lower() == "silica":
            toggle_cuvetteType(self, "silica")  # This is the last element and it is of type "silica"
            # return
    else:
        print()
    
    check_sample_type(self)
    print("="*50)

# Define item deletetion action
def on_item_delete(self, sample_code):
    index = self.solventName_dataSavingTab_comboBox.findText(sample_code)
    self.solventName_dataSavingTab_comboBox.removeItem(index)
    SM().unregister_sample(sample_code)

# Get values of Data Saving tab solventName QComboBox and insert into Data Fitting tab solventName QComboBox
# def sync_combos(self, source, target):
#     # Get current items in source combo
#     source_items = [source.itemText(i) for i in range(source.count())]

#     # Get current items in target combo
#     target_items = [target.itemText(i) for i in range(target.count())]

#     # Add items to target
#     for item in source_items:
#         if item not in target_items:
#             try:
#                 insert_index = target.findInsertIndex(item.lower())
#                 target.insertItem(insert_index, item)  # Assume `None` is okay for user data
#             except Exception as ex:
#                 print(f"Error inserting item: {ex}")  # Handle or log the error appropriately

#             index = target.findText(item)
#             target.setItemData(index, get_data(self, source))