from PyQt5.QtCore import pyqtSlot

from classes.mainwindow.measurement_plotting import plot_updated_data
from classes.mainwindow.measurement_processing import calculate_mean_run, create_string_list_to_display, update_fullLog, write_raw_data_log
from lib.gui_modifiers import GUI_Modifiers
from lib.worker import Worker

__all__ = ["print_output", "thread_complete", "thread_it"]


def print_output(self, returned_value):
    # print(returned_value)
    pass


def thread_complete(self):
    # print("THREAD COMPLETE!")
    pass


def thread_it(self, func_to_execute, *args, **kwargs):
    def check_function(self, func_to_execute):
        # Check directly in self
        if hasattr(self, func_to_execute.__name__):
            return getattr(self, func_to_execute.__name__) == func_to_execute

        # Now check in the nested mpositioner attribute
        if hasattr(self, "mpositioner"):
            if hasattr(self.mpositioner, func_to_execute.__name__):
                return getattr(self.mpositioner, func_to_execute.__name__) == func_to_execute

        # If none of those checks pass, return False or handle as needed
        return False

    # Pass the function to execute
    self.worker = Worker(func_to_execute, *args, **kwargs)  # Any other args, kwargs are passed to the run function
    self.worker.signals.result.connect(self.print_output)
    self.worker.signals.finished.connect(self.thread_complete)

    if check_function(self, func_to_execute):
        match func_to_execute:
            case self.mpositioner.movetostart:
                self.worker.signals.job_starting.connect(self.acquisition_starting)
                self.worker.signals.data_ready.connect(self.process_current_datapoint)
                self.worker.signals.acquisition_complete.connect(self.acquisition_complete)
                self.worker.signals.finished.connect(self.experiment_thread_finished)
                self.resume_scanning.connect(self.worker.resume)
            case self.mpositioner.movehome:
                label_base = f"Translation stage {self.motor_id} connected and "
                self.worker.signals.job_starting.connect(lambda: self.motor_label.setText(label_base + "homing..."))
                self.worker.signals.job_starting.connect(self.homing_start)
                self.worker.signals.finished.connect(lambda: self.motor_label.setText(label_base + "homed."))
                self.worker.signals.finished.connect(self.homing_complete)

    # Execute
    self.threadpool.start(self.worker)

    return self.worker


@pyqtSlot()
def homing_start(self):
    self.toggle_user_input_lock(False)


@pyqtSlot()
def homing_complete(self):
    self.toggle_measurement_buttons_lock(True)
    self.toggle_user_input_lock(True)
    GUI_Modifiers.switch_LED(self.initLED_pushButton, False)

    self.initializing = False
    self.initialized = True


@pyqtSlot()
def acquisition_starting(self):
    self.toggle_user_input_lock(False)
    self.toggle_measurement_buttons_lock(False)
    self.toggle_datasaving_buttons_lock(False, add_buttons=[self.saveData_pushButton])
    leds = [self.runLED_pushButton, self.waitLED_pushButton]
    GUI_Modifiers.switch_LED(leds, True)
    self.experiment_finished = False


@pyqtSlot(bool, int)
def acquisition_complete(self, *args):
    """When data acquisition single scan is finished"""
    self.data_acquisition_complete = args[0]
    run_no = args[1]  # run_no is in range from 0 to (number of scans - 1, including)

    # When a scan is finished
    if self.data_acquisition_complete:
        # When the last scan is finished
        if run_no == self.numberOfScans_spinBox.value() - 1:
            print("="*100)
            print("Finishing...")
            self.experiment_finished = True
            self.running = False

            # Calculate final mean values
            self.batch_data["Mean values"] = calculate_mean_run(self, self.batch_data, run_no, final=True)

            # Make sure that if "Mean values" is selected, it is updated in both `textBrowser`s
            string_list_to_display = create_string_list_to_display(self)
            write_raw_data_log(self, string_list_to_display)
            self.current_run += 1
            update_fullLog(self, update_button_clicked=False)
            self.current_run -= 1

            # Update the plot (important if "Mean" is selected)
            plot_updated_data(self, self.batch_data)

        self.update_labels()

    # While a scan is running
    # This condition is managed by data_ready signal in self.process_current_datapoint method
    else:
        pass


@pyqtSlot()
def experiment_thread_finished(self):
    self.experiment_stopped.clear()
    self.toggle_user_input_lock(True)
    self.toggle_measurement_buttons_lock(True)
    self.toggle_datasaving_buttons_lock(True)  # unlock only the "Update" button (forces user to click it to continue on saving the data)
    leds = [self.initLED_pushButton, self.runLED_pushButton, self.clearLED_pushButton, self.waitLED_pushButton]
    GUI_Modifiers.switch_LED(leds, False)
