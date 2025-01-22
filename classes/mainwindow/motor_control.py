# MOTOR NAVIGATION
import logging
import traceback

from classes.motor import motor as mot
from lib.settings import save_settings


def motor_detect(self, *args, **kwargs):
    self.motor_label.setText("Initializing translation stage...")
    self.motor_id = int(self.settings.value("Hardware/thorlabs_motor_id"))

    # only now (not on startup) import the package, because it takes a few seconds to load it
    from packages import thorlabs_apt as apt  # importing it from custom location allows to place APT.dll

    # in the package directory to read it (it is more user-friendly)
    try:
        self.motor = apt.Motor(self.motor_id)
        self.ocx.configure(self.motor_id)
        self.motor_label.setText(f"Translation stage {self.motor_id} connected.")

        if self.motor.get_stage_axis_info()[2] == 1:  # STAGE_UNITS_MM
            self.startPos_doubleSpinBox.setMinimum(self.motor.get_stage_axis_info()[0])
            self.startPos_doubleSpinBox.setMaximum(self.motor.get_stage_axis_info()[1])
            self.endPos_doubleSpinBox.setMinimum(self.motor.get_stage_axis_info()[0])
            self.endPos_doubleSpinBox.setMaximum(self.motor.get_stage_axis_info()[1])

    except Exception:
        logging.error(traceback.format_exc())
        print("Try closing other programs that may be accessing the device.")
        print("Try re-plugging the USB device and then run this program.")
        print(
            f"If you are simulating the motion controller, make sure Thorlabs APT Server has been configured with the HWSerialNumber {self.motor_id}"
        )
        self.showdialog("Error", "Motion controller not found! For more information, see the console output.")

    self.mpositioner = mot.MotorPositioner(self, self.motor)
    self.motor.backlash_distance = 0


def set_motor_home(self):
    self.thread_it(self.mpositioner.movehome)


def set_motor_to_start(self):
    # save settings so that if there is any crash afterwards, the settings are preserved
    save_settings(self, self.settings)
    self.experiment_stopped.clear()  # Reset stop event
    self.running = True
    self.current_run = 0
    start_pos = self.startPos_doubleSpinBox.value()
    end_pos = self.endPos_doubleSpinBox.value()

    if self.motor.position <= (start_pos + end_pos) / 2:
        self.motor.move_to(start_pos)
        self.scanning_from = "start"
    else:
        self.motor.move_to(end_pos)
        self.scanning_from = "end"

    self.mpositioner.resumed = False
    self.data = {"positions": [],
                 "relative": {k: [] for k in range(self.number_of_channels_used)},
                 "absolute": {k: [] for k in range(self.number_of_channels_used)}}
    self.thread_it(
        self.mpositioner.movetostart,
        self.experiment_stopped,
        detector_core_name=self.detector_core_name,
        data=self.data,
        batch_data=self.batch_data,
        number_of_scans=self.numberOfScans_spinBox.value(),
        scanning_from=self.scanning_from,
        number_of_steps=int(self.settings.value("MeasurementTab/steps_per_scan")),
        mpositioner=self.mpositioner,
        stepsize=(self.endPos_doubleSpinBox.value() - self.startPos_doubleSpinBox.value())
        / int(self.settings.value("MeasurementTab/steps_per_scan")),
    )


def stop_experiment(self):
    if hasattr(self, "motor"):
        self.experiment_stopped.set()  # Set the stop event to signal the worker to stop
        self.motor.stop_profiled()

        self.initializing = False
        self.running = False
        self.clearing = False

        self.clear_measurement()
