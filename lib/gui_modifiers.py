import sys

from PyQt5.QtCore import QSize
from PyQt5.QtWidgets import QApplication, QLayout, QListWidgetItem, QWidget

from lib.mgmotor import MG17Motor

try:
    from lib.window_gui import Ui_MainWindow  # type: ignore
except ModuleNotFoundError:
    # Attempt to compile the UI file
    from lib.compiling import compile_GUI
    from lib.exceptions import GUI_Error

    if not compile_GUI("./window.ui"):  # Replace './window.ui' with your actual UI file path
        # If compilation failed, show the error dialog
        app = QApplication(sys.argv)
        error_dialog = GUI_Error()
        error_dialog.show()
        sys.exit(app.exec_())
    else:
        # Compilation succeeded, attempt to import again
        from lib.window_gui import Ui_MainWindow  # type: ignore


class GUI_Modifiers(Ui_MainWindow):
    
    @staticmethod
    def switch_LED(led, switch):
        # Switch on/off if switch == True/False
        def switch_in_loop():
            for led_, switch_ in zip(led, switch):
                led_.setEnabled(switch_)

        if isinstance(led, list):
            if isinstance(switch, list):
                if len(led) == len(switch):
                    switch_in_loop()
            else:
                switch = [switch] * len(led)
                switch_in_loop()
        else:
            if isinstance(switch, list):
                led.setEnabled(switch[0])
            else:
                led.setEnabled(switch)

    @staticmethod
    def cover_widget(widget):
        if isinstance(widget, QLayout):
            for child in widget.parent().findChildren(QWidget):
                try:
                    csizepolicy = widget.sizePolicy()
                    csizepolicy.setRetainSizeWhenHidden(True)
                    child.setSizePolicy(csizepolicy)
                    child.hide()
                except AttributeError:
                    pass
        else:
            wsizepolicy = widget.sizePolicy()
            wsizepolicy.setRetainSizeWhenHidden(True)
            widget.setSizePolicy(wsizepolicy)
            widget.hide()

    @staticmethod
    def set_logging_items(widgets: list, final_number):
        size_hint = QSize(100, 30)  # Set the desired size width x height
        
        for widget in widgets:
            widget.clear()  # Clear existing items to avoid duplication
            
            # Insert additional items based on final_number
            if final_number > 1:
                # Create and insert the first item
                mean_item = QListWidgetItem("Mean values")
                mean_item.setSizeHint(size_hint)  # Set size hint for the Mean values item
                widget.insertItem(0, mean_item)  # Insert the Mean values item at the beginning
            
                for i in range(1, final_number + 1):
                    run_item = QListWidgetItem(f"Run {i}")  # Create new list items
                    run_item.setSizeHint(size_hint)  # Set size hint for each run item
                    widget.insertItem(i, run_item)  # Insert each run item accordingly
                
                # for i in range(widget.count()):
                #     widget.item(i).setData(Qt.UserRole, {"raw": "No data."})
            else:
                # Create and insert the first item
                mean_item = QListWidgetItem("Run 1")
                mean_item.setSizeHint(size_hint)  # Set size hint for the Mean values item
                widget.insertItem(0, mean_item)  # Insert the Mean values item at the beginning
                # widget.item(0).setData(Qt.UserRole, {"raw": "No data."})

            if final_number > 1:
                widget.show()
            else:
                widget.hide()
            
            widget.setCurrentRow(0)

    def toggle_user_input_lock(self, unlock):
        """Toggles the lock on widgets editable by user important for measurement"""
        elements = [
            self.extremes_radioButton,
            self.centerRange_radioButton,
            self.startPos_doubleSpinBox,
            self.focalPoint_doubleSpinBox,
            self.endPos_doubleSpinBox,
            self.zscanRange_measurementTab_doubleSpinBox,
            self.stepsScan_spinBox,
            self.samplesStep_spinBox,
            self.numberOfScans_spinBox,
        ]

        for element in elements:
            element.setEnabled(unlock)

        if unlock:
            self.toggle_measurement_range_selector()

    def toggle_measurement_buttons_lock(self, unlock):
        """Toggles the lock on measurement control buttons"""
        buttons = [self.run_pushButton, self.clear_pushButton, self.stopRed_pushButton]

        for button in buttons:
            button.setEnabled(unlock)
            if button is self.stopRed_pushButton:
                button.setEnabled(not unlock)
    
    def toggle_datasaving_buttons_lock(self, unlock, add_buttons=[]):
        """Toggles the lock on data saving buttons"""
        buttons = [self.update_pushButton]
        for button in add_buttons:
            buttons.append(button)

        for button in buttons:
            button.setEnabled(unlock)

    def toggle_measurement_range_selector(self):
        """Activate and deactivate range selectors based on checked radio button selector"""
        if self.extremes_radioButton.isChecked() is True:
            # enable first pair of settings
            self.startPos_doubleSpinBox.setEnabled(True)
            self.endPos_doubleSpinBox.setEnabled(True)
            # disable second pair of settings
            self.focalPoint_doubleSpinBox.setEnabled(False)
            self.zscanRange_measurementTab_doubleSpinBox.setEnabled(False)

        elif self.centerRange_radioButton.isChecked() is True:
            # disable first pair of settings
            self.startPos_doubleSpinBox.setEnabled(False)
            self.endPos_doubleSpinBox.setEnabled(False)
            # enable second pair of settings
            self.focalPoint_doubleSpinBox.setEnabled(True)
            self.zscanRange_measurementTab_doubleSpinBox.setEnabled(True)

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
                    models = ["SA", "2PA+SA"]  # , "RSA"]
                    if self.solventOA_absorptionModel_comboBox.currentText() in models:
                        self.solventOA_saturationModel_label.setVisible(True)
                        self.solventOA_saturationModel_comboBox.setVisible(True)
                    else:
                        self.solventOA_saturationModel_label.setVisible(False)
                        self.solventOA_saturationModel_comboBox.setVisible(False)
            case "Sample":
                if self.sampleOA_isAbsorption_checkBox.isChecked() is False:
                    models = ["SA", "2PA+SA"]  # , "RSA"]
                    if self.sampleOA_absorptionModel_comboBox.currentText() in models:
                        self.sampleOA_saturationModel_label.setVisible(True)
                        self.sampleOA_saturationModel_comboBox.setVisible(True)
                    else:
                        self.sampleOA_saturationModel_label.setVisible(False)
                        self.sampleOA_saturationModel_comboBox.setVisible(False)

    def add_nonqt_widgets(self):
        self.ocx = MG17Motor()
        ocx_layout = self.mg17motor_control_vlayout
        ocx_layout.addWidget(self.ocx)

    def update_labels(self):
        self.currentRun_label.setText(
            '<html><head/><body><p align="right">'
            f'<span style=" font-size:16pt; font-weight:600;"> {self.current_run + 1}'
            '</span></p></body></html>'
        )
        self.numberOfScans_label.setText(
            '<html><head/><body><p align="left">'
            f'<span style=" font-size:16pt; font-weight:600;">{self.numberOfScans_spinBox.value()} </span>'
            '</p></body></html>'
        )