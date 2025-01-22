import time
import winsound

# import nidaqmx
# import numpy as np
# from nidaqmx import stream_readers  # noqa: F401
# from nidaqmx.constants import AcquisitionType, Edge
from classes.DAQ.manager import perform_acquisition
from PyQt5.QtCore import QObject

from lib.exceptions import ExperimentStopped


class MotorPositioner(QObject):
    def __init__(self, parent_window, motor):
        self.parent_window = parent_window
        self.motor = motor
        self.stop_flag = None
        self.resumed = False

    def movetostart(self, *args, **kwargs):
        self.stop_flag = args[0]

        if self.motor.is_in_motion is True:
            print("Moving to starting position")

        while self.motor.is_in_motion is True:
            continue

        self.run(**kwargs)

    def moveby(self, move_step):
        self.motor.move_by(move_step, blocking=True)

    def movehome(self, *args, **kwargs):
        if self.motor.has_homing_been_completed is False:
            self.motor.move_home(blocking=True)
            time.sleep(0.2)  # otherwise it may not move_home()
            if self.motor.is_in_motion is True:
                print("Homing now")
            while self.motor.has_homing_been_completed is False:
                kwargs["motion_finished_callback"].emit(0)  # wait until homing is completed
            time.sleep(0.2)  # wait a little more (so the motor.position gets exactly "0" position)

        kwargs["motion_finished_callback"].emit(1)

    def run(self, **kwargs):
        """This function collects all the data from measurement. The function runs in a separate thread so that
        freezes of GUI such as dialog boxes doesn't affect the measurement process."""
        self.parent_window.data_acquisition_complete = False
        samples_per_pos = int(self.parent_window.settings.value("MeasurementTab/samples_per_position"))
        number_of_steps = int(self.parent_window.settings.value("MeasurementTab/steps_per_scan"))

        laser_rep_rate = float(self.parent_window.settings.value("Hardware/laser_repetition_rate"))
        dig_edge_src = self.parent_window.settings.value("Hardware/nidaqmx_dig_edge_src")
        number_of_channels = self.parent_window.number_of_channels_used

        slp_time = samples_per_pos / laser_rep_rate
        time.sleep(slp_time)  # Sometimes the first datapoint is collected before the motor has settled

        kwargs["batch_data"]["merged"] = {"positions": []}
        for sub_key in ["relative", "absolute"]:
            kwargs["batch_data"]["merged"][sub_key] = {
                i: [] for i in range(number_of_channels)
            }  # this is before the batch loop to merge on the fly

        try:
            for run_no in range(kwargs["number_of_scans"]):
                if self.stop_flag.is_set():
                    raise ExperimentStopped

                self.parent_window.current_run = run_no
                for step in range(number_of_steps + 1):
                    if self.stop_flag.is_set():
                        raise ExperimentStopped

                    ### Create a task - OLD VERSION
                    # with nidaqmx.Task() as task:
                    #     # Create MultiChannel "channel"
                    #     task.ai_channels.add_ai_voltage_chan(self.parent_window.detector_core_name+f"0:{number_of_channels}") # "Dev1/ai0:3"

                    #     # Start Digital Edge
                    #     task.triggers.start_trigger.dig_edge_src = dig_edge_src
                    #     task_trigger_src = task.triggers.start_trigger.dig_edge_src

                    #     # Sample Clock
                    #     rate0 = laser_rep_rate
                    #     task.timing.cfg_samp_clk_timing(
                    #         rate0,
                    #         source=task_trigger_src,
                    #         active_edge=Edge.FALLING,
                    #         sample_mode=AcquisitionType.FINITE,
                    #         samps_per_chan=samples_per_pos
                    #         )

                    #     sampling_period = samples_per_pos/rate0 # THIS IS THE TIME NEEDED TO COLLECT ALL SAMPLES

                    #     ### Data acquisition
                    #     reader = nidaqmx.stream_readers.AnalogMultiChannelReader(task.in_stream)
                    #     values_read = np.zeros((number_of_channels+1,samples_per_pos),dtype=np.float64)
                    #     # Start Task
                    #     task.start()

                    #     time.sleep(sampling_period)

                    #     # Acquire data
                    #     self.parent_window.data["positions"].append(self.motor.position)
                    #     reader.read_many_sample(values_read, number_of_samples_per_channel=samples_per_pos)

                    #     # Take mean for each channel and append to "data" dictionary
                    #     data_mean = np.mean(values_read,axis=1)

                    #     for chan_no in range(number_of_channels):
                    #         self.parent_window.data["absolute"][chan_no].append(data_mean[chan_no])
                    #     # These loops need to be separated because "absolute" list at current step MUST BE defined at current step when accessed by "relative" list
                    #     for chan_no in range(number_of_channels):
                    #         # Take last/current value (at index [step]) in "absolute" values list and divide by channel [1] (reference)
                    #         # And append to "relative" values list
                    #         self.parent_window.data["relative"][chan_no].append(data_mean[chan_no]/self.parent_window.data["absolute"][1][step])

                    #     # Merge in batch_data (for display of all data on the fly over multiple runs)
                    #     self.parent_window.batch_data['merged']['positions'].append(self.parent_window.data["positions"][-1])
                    #     for sub_key in ['relative', 'absolute']:
                    #         for chan_no in range(number_of_channels):
                    #             self.parent_window.batch_data['merged'][sub_key][chan_no].append(self.parent_window.data[sub_key][chan_no][-1])

                    #     # Save current data to batch_data dictionary
                    #     old = copy.deepcopy(self.parent_window.data)
                    #     self.parent_window.batch_data[run_no] = old

                    #     # 2) display data in GUI thread
                    #     kwargs['progress_callback'].emit(step)
                    #     kwargs['data_ready_callback'].emit(run_no, step, old, scanning_from, step==number_of_steps)

                    #     # 3) move the motor
                    #     if step < number_of_steps:
                    #         self.parent_window.mpositioner.moveby(stepsize)

                    #     task.stop()

                    # kwargs['acquisition_complete'].emit(step==number_of_steps, run_no)

                    ## NEW VERSION (to be checked)
                    self.resumed = False
                    perform_acquisition(self, run_no, step, number_of_channels, samples_per_pos, laser_rep_rate, dig_edge_src, kwargs)
                    while not self.resumed:  # wait for GUI ready
                        time.sleep(0.1)
                    # Move the motor if more steps are available
                    if step < kwargs["number_of_steps"]:
                        kwargs["mpositioner"].moveby(kwargs["stepsize"])

                if kwargs["number_of_scans"] > 1:
                    # CHANGE DIRECTION OF THE SCAN
                    if kwargs["scanning_from"] == "end":
                        kwargs["scanning_from"] = "start"
                    else:
                        kwargs["scanning_from"] = "end"

                    # Clear current data if the loop is still running
                    if run_no < kwargs["number_of_scans"] - 1:
                        kwargs["data"]["positions"] = []  # reset positions to empty list
                        for typ, _ in self.parent_window.charts.items():  # e.g.: take the tuple ("relative", "self.rel_chart")
                            for chan_no in range(number_of_channels):
                                kwargs["data"][typ][chan_no] = []

                kwargs["acquisition_complete"].emit(step == kwargs["number_of_steps"], run_no)  # run_no is in range from 0 to (number of scans - 1, including)

            # EMIT SOUND AT FINISH
            duration = int(self.parent_window.settings.value("Perks/end_beep_duration_ms"))
            freq = int(self.parent_window.settings.value("Perks/end_beep_tone_frequency"))
            if self.parent_window.settings.value("Perks/end_beep_emit") == "True":
                for _ in range(3):
                    winsound.Beep(freq, duration)
                    time.sleep(0.05)

            return

        except ExperimentStopped:
            print("Experiment stopped with red button.")
            return
