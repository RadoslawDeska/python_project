import logging
import time

import numpy as np


def clear_measurement(self):
    """clears all in the first two tabs (Measurement and Data Saving)"""
    self.clearing = True
    self.data_acquisition_complete = False
    self.current_run = 0

    self.data["positions"] = []  # reset positions to empty list
    self.batch_data = {}
    
    # SORTING VARIABLES
    self.already_sorted = {}
    self.all_positions = []
    self.all_absolutes = {k: [] for k in range(self.number_of_channels_used)}
    self.all_relatives = {k: [] for k in range(self.number_of_channels_used)}

    self.rms_value = 0.0
    self.rms_text.txt.set_text(f"RMS noise = {self.rms_value*100:.3f}%")

    for signal_type, chart in self.charts.items():  # e.g.: take the tuple ("relative", "self.rel_chart")
        for chan_no in range(self.number_of_channels_used):
            self.data[signal_type][chan_no] = []  # then fill "relative" dictionary in "self.data" dictionary empty list of values per channel

            # UPDATE LINES INSTEAD OF DELETING AND REINSTANTIATING
            self.measurement_lines[signal_type][chan_no]["current"].set_xdata(
                self.data["positions"]
            )  # and set data of lines in "lines" dictionary to empty lists of positions and values
            self.measurement_lines[signal_type][chan_no]["current"].set_ydata(self.data[signal_type][chan_no])
            self.measurement_lines[signal_type][chan_no]["mean"].set_xdata(
                self.data["positions"]
            )  # and set data of lines in "lines" dictionary to empty lists of positions and values
            self.measurement_lines[signal_type][chan_no]["mean"].set_ydata(self.data[signal_type][chan_no])

        chart.axes.relim()
        chart.axes.autoscale_view()
        chart.draw_idle()

    self.update_labels()
    self.set_logging_items((self.rawData_listWidget, self.fullExp_listWidget), self.numberOfScans_spinBox.value())
    self.rawData_textBrowser.clear()
    self.fullLogData_textBrowser.clear()

    self.toggle_datasaving_buttons_lock(False, add_buttons=[self.saveData_pushButton])

    self.clearing = False


def plot_updated_data(self, batch_data):
    if not hasattr(self, 'number_of_channels_used'):
        return
    
    y = []  # Initialize y to ensure it is defined
    showscans = "Current" if self.current_run == 0 else self.showScans_comboBox.currentText()

    try:
        if showscans == "All":
            for k in self.charts.keys():
                for chan_no in range(self.number_of_channels_used):
                    self.measurement_lines[k][chan_no]["current"].set_xdata(batch_data["merged"]["positions"])
                    self.measurement_lines[k][chan_no]["current"].set_ydata(batch_data["merged"][k][chan_no])

                    self.measurement_lines[k][chan_no]["mean"].set_xdata([])
                    self.measurement_lines[k][chan_no]["mean"].set_ydata([])

            y = batch_data["merged"]["absolute"][1]

        elif showscans == "Mean":  # Active only after the first run (self.current_run > 1)
            for k in self.charts.keys():
                for chan_no in range(self.number_of_channels_used):
                    cur_pos = self.data["positions"]
                    cur_data = self.data[k][chan_no]
                    mn_current = min(len(cur_data), len(cur_pos))
                    
                    mean_pos = batch_data["Mean values"]["positions"]
                    mean_data = batch_data["Mean values"][k][chan_no]
                    mn_mean = min(len(mean_data), len(mean_pos))

                    if self.scanning_from == "start":
                        self.measurement_lines[k][chan_no]["current"].set_xdata(cur_pos[:mn_current])
                        self.measurement_lines[k][chan_no]["current"].set_ydata(cur_data[:mn_current])
                        
                        self.measurement_lines[k][chan_no]["mean"].set_xdata(mean_pos[:mn_mean])
                        self.measurement_lines[k][chan_no]["mean"].set_ydata(mean_data[:mn_mean])
                    else:
                        self.measurement_lines[k][chan_no]["current"].set_xdata(cur_pos[:mn_current])
                        self.measurement_lines[k][chan_no]["current"].set_ydata(cur_data[:mn_current])
                        
                        self.measurement_lines[k][chan_no]["mean"].set_xdata(mean_pos[-mn_mean:][::-1])
                        self.measurement_lines[k][chan_no]["mean"].set_ydata(mean_data[-mn_mean:][::-1])

                    if k == "absolute" and chan_no == 1:
                        y = batch_data["Mean values"][k][chan_no]

        elif showscans == "Current":
            for k in self.charts.keys():
                for chan_no in range(self.number_of_channels_used):
                    self.measurement_lines[k][chan_no]["current"].set_xdata(self.data["positions"])
                    self.measurement_lines[k][chan_no]["current"].set_ydata(self.data[k][chan_no])

                    self.measurement_lines[k][chan_no]["mean"].set_xdata([])
                    self.measurement_lines[k][chan_no]["mean"].set_ydata([])

            y = self.data["absolute"][1]
        
        for chart in self.charts.values():
            chart.draw_idle()  # Only redraw the canvas once after all updates
        
    except (KeyError, ValueError) as e:
        logging.error(f"Error in processing data: {str(e)}")
        # return  # Return here to avoid undefined behavior later
    else:
        if len(y) > 0:
            try:
                self.rms_value = np.linalg.norm([yi - np.mean(y) for yi in y]) / np.sqrt(len(y))  # RMS calculation
                self.rms_text.txt.set_text(f"RMS noise = {self.rms_value*100:.3f}%")
            except IndexError as e:
                logging.error(f"Error calculating RMS: {str(e)}")
            except Exception as e:
                logging.error(f"Unexpected error: {str(e)}")
        tic = time.perf_counter_ns()
        rescale_measurement_plots(self, self.focusAt_comboBox.currentText())
        toc = time.perf_counter_ns()
        print(f"rescale_measurement_plots() executed in {(toc-tic)/1E6} ms" )


def rescale_measurement_plots(self, who_called=""):
    """Rescales plots in the 'Measurement' Tab based on
    values given in the 'Measurement control' panel\n
    If user changes the focus of chart at specific line, rescaling fits the charts
    to display the line in its min-max y-range."""

    s = self.startPos_doubleSpinBox.value()
    e = self.endPos_doubleSpinBox.value()
    f = self.focalPoint_doubleSpinBox.value()
    z = self.zscanRange_measurementTab_doubleSpinBox.value()

    smin = self.startPos_doubleSpinBox.minimum()
    emax = self.endPos_doubleSpinBox.maximum()

    def rescale_horizontally(who_called):
        nonlocal s, e, f, z  # Declaring these variables as nonlocal so we can modify them

        # RESCALE HORIZONTALLY (upon start/end and focal/range change values)
        match who_called:
            case "start":
                if s >= e:
                    self.startPos_doubleSpinBox.setValue(self.previous_start_pos)
            case "end":
                if e <= s:
                    self.endPos_doubleSpinBox.setValue(self.previous_end_pos)
            # The 'focal' and 'range' cases ensure that start and end positions are in range of motion
            case "focal":
                if f + z / 2 > emax:
                    self.zscanRange_measurementTab_doubleSpinBox.setValue(2 * (emax - f))
                    z = self.zscanRange_measurementTab_doubleSpinBox.value()
                elif f - z / 2 < smin:
                    self.zscanRange_measurementTab_doubleSpinBox.setValue(2 * f)
                    z = self.zscanRange_measurementTab_doubleSpinBox.value()
            case "range":
                if f + z / 2 > emax:
                    self.focalPoint_doubleSpinBox.setValue(emax - z / 2)
                    f = self.focalPoint_doubleSpinBox.value()
                elif f - z / 2 < smin:
                    self.focalPoint_doubleSpinBox.setValue(z / 2)
                    f = self.focalPoint_doubleSpinBox.value()

        # RESCALE HORIZONTALLY
        if self.extremes_radioButton.isChecked() is True:
            self.focalPoint_doubleSpinBox.setValue((e + s) / 2)
            self.zscanRange_measurementTab_doubleSpinBox.setValue(np.abs(e - s))
        elif self.centerRange_radioButton.isChecked() is True:
            self.startPos_doubleSpinBox.setValue(f - z / 2)
            self.endPos_doubleSpinBox.setValue(f + z / 2)

        # refresh variables content after changes
        s = self.startPos_doubleSpinBox.value()
        e = self.endPos_doubleSpinBox.value()
        f = self.focalPoint_doubleSpinBox.value()
        z = self.zscanRange_measurementTab_doubleSpinBox.value()

        self.offset = (e - s) / self.stepsScan_spinBox.value()

        try:
            for chart in self.charts.values():
                chart.axes.set_xlim(s - self.offset, e + self.offset)
        except ValueError:
            pass

    def rescale_vertically():
        try:
            if self.initialized is False:  # prevent the problem of accessing data for rescale when there is no data created yet.
                return

            focus_at = self.focusAt_comboBox.currentText()
            show_scans = self.showScans_comboBox.currentText()

            for signal_type, chart in self.charts.items():
                match focus_at:
                    case "All":
                        set_new_limits(self, -1, chart=chart)
                    case "Closed":
                        set_new_limits(self, 0, show_scans=show_scans, signal_type=signal_type, chart=chart)
                    case "Reference":
                        set_new_limits(self, 1, show_scans=show_scans, signal_type=signal_type, chart=chart)
                    case "Open":
                        set_new_limits(self, 2, show_scans=show_scans, signal_type=signal_type, chart=chart)

        except ValueError:
            pass

    if who_called in ["start", "end", "focal", "range"]:
        rescale_horizontally(who_called)
    else:  # if the who_called value is any of the values in the list of "Focus At" dropdown OR "ANYTHING ELSE"
        rescale_vertically()


def set_new_limits(self, channel_number, **kwargs):
    try:
        chart = kwargs["chart"]

        if channel_number == -1:
            chart.axes.relim()
            chart.axes.autoscale(axis="y")
            return
        else:
            show_scans = kwargs["show_scans"]
            signal_type = kwargs["signal_type"]

            match show_scans:
                case "All" | "Current":
                    data = self.measurement_lines[signal_type][channel_number]["current"].get_ydata()
                    mn = min(data)
                    mx = max(data)
                case "Mean":
                    data = self.measurement_lines[signal_type][channel_number]["current"].get_ydata()
                    mn_current = min(data)
                    mx_current = max(data)
                    data = self.measurement_lines[signal_type][channel_number]["mean"].get_ydata()
                    mn_mean = min(data)
                    mx_mean = max(data)

                    mn = min(mn_current, mn_mean)
                    mx = max(mx_current, mx_mean)

            if mn == mx:  # prevent UserWarning: Attempting to set identical low and high ylims
                mn *= 0.95
                mx *= 1.05
            chart.axes.set_ylim(mn, mx)

    except ValueError:
        return
