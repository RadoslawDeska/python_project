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