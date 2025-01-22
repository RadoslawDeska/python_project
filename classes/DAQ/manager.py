import copy
import logging
import random
import time
import traceback

import nidaqmx
import nidaqmx.stream_readers
import numpy as np


class DaqTaskManager:
    def __init__(self, detector_core_name, number_of_channels, laser_rep_rate, samples_per_pos, dig_edge_src):
        self.detector_core_name = detector_core_name
        self.number_of_channels = number_of_channels
        self.laser_rep_rate = laser_rep_rate
        self.samples_per_pos = samples_per_pos
        self.dig_edge_src = dig_edge_src
        self.task = nidaqmx.Task()

    def __enter__(self):
        self.task.ai_channels.add_ai_voltage_chan(f"{self.detector_core_name}0:{self.number_of_channels - 1}")

        # Setup triggers
        self.task.triggers.start_trigger.dig_edge_src = self.dig_edge_src

        # Configure sample clock
        self.task.timing.cfg_samp_clk_timing(
            self.laser_rep_rate,
            source=self.dig_edge_src,
            active_edge=nidaqmx.constants.Edge.FALLING,
            sample_mode=nidaqmx.constants.AcquisitionType.FINITE,
            samps_per_chan=self.samples_per_pos,
        )

        return self

    def read_data(self):
        values_read = np.zeros((self.number_of_channels, self.samples_per_pos), dtype=np.float64)
        reader = nidaqmx.stream_readers.AnalogMultiChannelReader(self.task.in_stream)
        self.task.start()

        # Wait for data acquisition to complete
        time.sleep(self.samples_per_pos / self.laser_rep_rate)
        reader.read_many_sample(values_read, number_of_samples_per_channel=self.samples_per_pos)
        self.task.stop()

        return values_read

    def __exit__(self, exc_type, exc_value, traceback):
        self.task.close()


def perform_acquisition(self, run_no, step, number_of_channels, samples_per_pos, laser_rep_rate, dig_edge_src, kwargs: dict):
    keys = [
        "detector_core_name",
        "data",
        "batch_data",
        "progress_callback",
        "data_ready_callback",
        "scanning_from",
        "number_of_steps",
        "mpositioner",
        "acquisition_complete",
        "stepsize",
    ]

    try:
        if len(kwargs["data"]["positions"]) == 0:
            kwargs["data"]["positions"] = [self.motor.position]
        else:
            kwargs["data"]["positions"].append(self.motor.position)

        for key in keys:
            if key not in kwargs:
                print("Key missing:", key)
                raise Exception

    except Exception:
        print("Missing keys in kwargs of perform_acquisition method.")
        logging.error(traceback.format_exc())

    stepsize = kwargs["stepsize"]
    scanning_from = kwargs["scanning_from"]
    if (scanning_from == "end" and stepsize > 0) or (scanning_from == "start" and stepsize < 0):
        stepsize *= -1

    with DaqTaskManager(
        detector_core_name=kwargs["detector_core_name"],
        number_of_channels=number_of_channels,
        laser_rep_rate=laser_rep_rate,
        samples_per_pos=samples_per_pos,
        dig_edge_src=dig_edge_src,
    ) as task_manager:
        values_read = task_manager.read_data()

        # Data Processing
        data_mean = np.mean(values_read, axis=1)
        
        # Try this to simulate the job
        data_mean[0] = (1/((self.motor.position-58)**2+10)-1/((self.motor.position-52)**2+10))+0.21
        data_mean[0] = random.gauss(data_mean[0], 0.003053)
        data_mean[1] = 0.198
        data_mean[1] = random.gauss(data_mean[1], 0.000682)
        data_mean[2] = 0.20
        data_mean[2] = random.gauss(data_mean[2], 0.001087)

        try:
            for chan_no in range(number_of_channels):
                kwargs["data"]["absolute"][chan_no].append(data_mean[chan_no])
        
            for chan_no in range(number_of_channels):
                try:
                    kwargs["data"]["relative"][chan_no].append(data_mean[chan_no] / kwargs["data"]["absolute"][1][step])
                except IndexError:
                    raise
                except ZeroDivisionError:
                    kwargs["data"]["relative"][chan_no].append(None)
        
            kwargs["batch_data"]["merged"]["positions"].append(kwargs["data"]["positions"][-1])
            
            for sub_key in ["relative", "absolute"]:
                for chan_no in range(number_of_channels):
                    kwargs["batch_data"]["merged"][sub_key][chan_no].append(kwargs["data"][sub_key][chan_no][-1])
        
        except IndexError:  # known to happen on experiment interruption
            logging.error(traceback.format_exc())
        
        old = copy.deepcopy(kwargs["data"])
        kwargs["batch_data"][f"Run {run_no + 1}"] = old

        # Emit signals for GUI
        kwargs["stepsize"] = stepsize
        kwargs["data_ready_callback"].emit(run_no, step, kwargs)  # run_no is in range from 0 to (number of scans - 1, including)
