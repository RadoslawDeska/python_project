"""This file contains methods responsible for:
- displaying the data incoming from currently running measurement
- saving measurement data after measurement is over
"""

import logging
import os
import time
import traceback
import tracemalloc

from datetime import datetime

import numpy as np

from classes.mainwindow.measurement_plotting import plot_updated_data
from lib.functions import get_base_directory, is_not_ragged

__all__ = ["process_current_datapoint", "update_fullLog", "data_save"]

tracemalloc.start()


def switch_widget_items(self, item, write_full_log=True):
    parent_widget = item.listWidget()
    selection_row = parent_widget.row(item)

    # Select the same item in the list of the second data tab
    if parent_widget is self.rawData_listWidget:
        self.fullExp_listWidget.setCurrentItem(self.fullExp_listWidget.item(selection_row))
    elif parent_widget is self.fullExp_listWidget:
        self.rawData_listWidget.setCurrentItem(self.rawData_listWidget.item(selection_row))

    string_for_selected_item = create_string_list_to_display(self)
    write_raw_data_log(self, string_for_selected_item)

    if write_full_log:
        update_fullLog(self, update_button_clicked=False)


def sort_batch_data(self, batch_data, current_run_no):
    """
    Sorting batch data only for completed runs
    """

    for run_no in range(current_run_no + 1):  # make sure that the function runs if current_run_no == 0
        key = f"Run {run_no + 1}"
        if key in batch_data.keys():
            positions = batch_data[key]["positions"]
            absolute = batch_data[key]["absolute"]
            relative = batch_data[key]["relative"]

            if run_no == current_run_no:  # currently acquired data
                if len(positions) != self.stepsScan_spinBox.value() + 1:  # don't sort if incomplete
                    print(f'Data not complete. Not sorting yet. "Run {run_no + 1}".')
                    continue

            if (
                key in self.already_sorted and len(self.already_sorted[key]["positions"]) == self.stepsScan_spinBox.value() + 1
            ):  # don't sort if data is complete and sorted
                print(f'Data complete and sorted for "{key}".')
                continue

            print(f'Data complete and unsorted for "{key}". Sorting...')
            # Transpose all channels absolute and relative values into lists of tuples
            absolute_values = [tuple(item) for item in zip(*absolute.values())]
            relative_values = [tuple(item) for item in zip(*relative.values())]

            # Create final tuples by zipping positions with absolute and relative values
            combined = [(pos, *abs_val, *rel_val) for pos, abs_val, rel_val in zip(positions, absolute_values, relative_values)]

            sorted_combined = sorted(combined, key=lambda x: x[0])
            positions_sorted, absolute_sorted, relative_sorted = zip(*[(item[0], item[1:4], item[4:]) for item in sorted_combined])

            self.already_sorted[key] = {
                "positions": list(positions_sorted),
                "absolute": {i: list(entry) for i, entry in enumerate(zip(*absolute_sorted))},
                "relative": {i: list(entry) for i, entry in enumerate(zip(*relative_sorted))},
            }


def calculate_mean_run(self, batch_data, current_run_no, final=False):
    """
    Calculate the mean values of channel data across multiple runs for given positions.

    This method computes the arithmetic mean of channel data for all runs up to the
    provided `current_run_no`. It ensures that the first run is included in the calculation.
    The positions are sorted in ascending order, and mean values are computed for both
    relative and absolute data across the specified channels.

    Parameters:
    -----------
    batch_data : List[Dict]
        A list of dictionaries where each dictionary contains:
          - "positions": A list of positions (e.g., numerical coordinates)
          - "relative": A list of lists containing relative channel data for each channel
          - "absolute": A list of lists containing absolute channel data for each channel

    current_run_no : int
        The index of the current run. This determines which runs are included in the mean calculation.
        Note that runs are indexed starting from 0.

    Returns:
    --------
    mean_run_data : Dict
        A dictionary containing the mean run data structured as follows:
        {
            "positions": List,
            "relative": {
                channel_no: List,  # Mean values for the relative channel data
            },
            "absolute": {
                channel_no: List,  # Mean values for the absolute channel data
            }
        }
        - "positions": A sorted list of positions from the first run.
        - "relative": A dictionary mapping each channel index to its mean relative values.
        - "absolute": A dictionary mapping each channel index to its mean absolute values.

    Notes:
    -------
    - This method assumes that the structure of `batch_data` is consistent with the description above.
    - The mean is calculated only for fully available datasets in previous runs. If a run is incomplete
      (missing data), it will not be considered in the mean computation.
    - If `current_run_no` is 0, the method will only include the data from the first run.
    """

    mean_run_data = {"positions": [], "relative": {}, "absolute": {}}

    sort_batch_data(self, batch_data, current_run_no)

    # error_ = False
    for key in self.already_sorted:
        for ch in range(self.number_of_channels_used):
            if (
                is_not_ragged(
                    [self.already_sorted[key]["positions"], self.already_sorted[key]["relative"][ch], self.already_sorted[key]["absolute"][ch]]
                )
                is False
            ):
                print("number of datapoints - positions:", len(self.already_sorted[key]["positions"]))
                print(f"number of datapoints - relative channel {ch}: {len(self.already_sorted[key]['relative'][ch])}")
                print(f"number of datapoints - absolute channel {ch}: {len(self.already_sorted[key]['absolute'][ch])}")
                print(f"{key} {ch} - unequal lengths of sorted data.")

                # error_ = True

    if current_run_no == 0 and final is False:  # don't calculate mean of the first data
        return batch_data[f"Run {current_run_no + 1}"]
    elif final is True:
        # print("Calculating final mean")
        pass
    else:  # continue execution
        pass

    if not self.all_positions:  # preallocate the data for current run
        self.all_positions = []
        self.all_absolutes = {k: [] for k in range(self.number_of_channels_used)}
        self.all_relatives = {k: [] for k in range(self.number_of_channels_used)}

    # Collect positions and data for runs less than or equal to the current run number
    for run_no in range(current_run_no + 1):
        if f"Run {run_no + 1}" in self.already_sorted:
            self.all_positions.append(self.already_sorted[f"Run {run_no + 1}"]["positions"])
            for k in self.all_absolutes:
                self.all_absolutes[k].append(self.already_sorted[f"Run {run_no + 1}"]["absolute"][k])
                self.all_relatives[k].append(self.already_sorted[f"Run {run_no + 1}"]["relative"][k])

    # Calculate means based on the collected data
    try:
        mean_run_data["positions"] = np.mean(self.all_positions, axis=0).tolist()
    except ValueError:
        print("Data is missing. Trying again.")
        calculate_mean_run(self, batch_data, current_run_no, final=False)
        print("Data is missing, filling with NaN's.")
        max_length = max(len(sublist) for sublist in self.all_positions)
        nan_padded_positions = np.array([sublist + [np.nan] * (max_length - len(sublist)) for sublist in self.all_positions])
        mean_run_data["positions"] = np.nanmean(nan_padded_positions, axis=0).tolist()

    for k in range(self.number_of_channels_used):
        mean_run_data["absolute"][k] = np.mean(self.all_absolutes[k], axis=0).tolist()
        mean_run_data["relative"][k] = np.mean(self.all_relatives[k], axis=0).tolist()

    return mean_run_data


def process_current_datapoint(self, *args):
    """
    Processes a current data point as part of a scanning process within a batch data context.

    This method handles the processing of the current run, updating corresponding data structures
    and visualizations, saving results, and managing the state of the scanning operation.

    Parameters:
        *args: Unpacked arguments that are expected to be provided in the following order:
            - current_run_no (int): The index of the current run being processed.
            - current_step (int): The current step in the scanning process.
            - current_kwargs (dict): A dictionary containing keys:
                - "batch_data" (dict): The current batch data being processed, which includes data relevant to the runs.
                - "number_of_steps" (int): Total number of steps for the current run, used to determine if the run is finished.
                - "scanning_from" (str): The initial scanning point or identifier.

    Returns:
        None

    Exceptions:
        This method logs any exceptions that occur during execution using Python's logging module
        and stack trace capturing with traceback.format_exc().

    Functionality:
        - Determines if the current run has finished based on the provided step and number of steps.
        - Updates current batch data to include the mean of the run using the `calculate_mean_run` utility.
        - Plots the updated data using `plot_updated_data`.
        - Constructs a string representation of the current data point using `create_string_list_from_data`.
        - Writes the constructed string to a temporary file with `write_temp_file`.
        - Saves the string to appropriate list widget item data, updating both the current run and
          the mean run data.
        - Logs the raw data using the `write_raw_data_log`.
        - Empties the temporary lines if the current run is finished.
        - Deep copies the current batch data to the instance's `batch_data` attribute.

    Notes:
        - The method triggers the emission of a scanning resume signal at the end, allowing any
          connected observers to react to the completion of the processing step.
        - Ensure that the list widget (self.rawData_listWidget) is properly initialized and contains
          the appropriate number of items before attempting to access its child items.

    Example:
        # Example usage of this method within a scanning process context
        self.process_current_datapoint(run_no, step_no, {
            "batch_data": batch_data,
            "number_of_steps": total_steps,
            "scanning_from": scanning_from_identifier
        })
    """
    try:
        current_run_no, current_step, current_kwargs = args  # current_run_no is in range from 0 to (number of scans - 1, including)
        current_batch_data = current_kwargs["batch_data"]
        current_scanning_from = current_kwargs["scanning_from"]
        self.scanning_from = current_scanning_from

        self.batch_data[f"Run {current_run_no + 1}"] = current_batch_data[f"Run {current_run_no + 1}"]

        tic = time.perf_counter_ns()
        # Get mean values from previous runs or calculate the mean
        if current_run_no > 0:
            # calculate the mean only when data for current run is complete
            if current_step == self.stepsScan_spinBox.value():
                self.batch_data["Mean values"] = calculate_mean_run(self, current_batch_data, current_run_no)
            else:  # don't change the currently stored value of mean values
                print("Not changing currently stored mean values")
        else:  # take the first run as the values for mean
            print("Taking first run data as mean values")
            self.batch_data["Mean values"] = self.batch_data["Run 1"]

        toc = time.perf_counter_ns()
        self.meanruntook.setText(f"Calculating mean took {(toc-tic)/1E6} ms.")

        # Lazily update the plot (only when in view, or at the experiment end)
        tic = time.perf_counter_ns()
        if self.mainTabs.currentIndex() == 0 or self.experiment_finished:
            plot_updated_data(self, self.batch_data)
            toc = time.perf_counter_ns()
            print(f"Plotting took {(toc-tic)/1E6} ms.")
        else:
            toc = time.perf_counter_ns()
            print(f"Plotting took {(toc-tic)/1E6} ms.")

        tic = time.perf_counter_ns()
        string_list_to_write = create_string_list_for_temp_file(self, current_run_no)
        toc = time.perf_counter_ns()
        print(f"Creating string for temp file took {(toc-tic)/1E6} ms.")

        tic = time.perf_counter_ns()
        write_temp_file(self, current_run_no, current_step, string_list_to_write)
        toc = time.perf_counter_ns()
        print(f"Writing temp file took {(toc-tic)/1E6} ms.")

        tic = time.perf_counter_ns()

        # Lazily update the entry in raw data log (only when in view, or at the experiment end)
        if (self.mainTabs.currentIndex() == 1 and self.datalogTabs.currentIndex() == 0) or self.experiment_finished:
            tic1 = time.perf_counter_ns()
            string_list_to_display = create_string_list_to_display(self)
            toc1 = time.perf_counter_ns()
            print(f"Creating string list to display took {(toc1-tic1)/1E6} ms.")

            tic2 = time.perf_counter_ns()
            write_raw_data_log(self, string_list_to_display)
            toc2 = time.perf_counter_ns()
            print(f"Writing raw data log took {(toc2-tic2)/1E6} ms.")
        else:
            toc = time.perf_counter_ns()
            print(f"Creating string list to display and writing raw data log took {(toc-tic)/1E6} ms.")

    except Exception:
        logging.error(traceback.format_exc())

    self.resume_scanning.emit()

    print("\n" + "=" * 100 + "\n" * 2)


def create_line(self, data: dict, current_step: int = 0) -> str:
    if data["positions"][-1] < data["positions"][0]:  # scanning_from = "end"
        step_to_write = np.abs(current_step - self.stepsScan_spinBox.value())
    else:
        step_to_write = current_step

    line = f"{step_to_write+1:4d}{' '}"

    # Format the current_step data
    for chan_no in data["absolute"].keys():
        line += f"{2*' '}{data['absolute'][chan_no][current_step]:12.8f}{9*' '}"
    # Add the zero column to fit the formatting of older files
    line += f"{2*' '}{0:12.8f}"
    line += "\n"

    return line


def create_string_list_for_temp_file(self, current_run):
    try:
        # Get the current run data
        itemdata = self.batch_data[f"Run {current_run + 1}"]

        # Convert the item data to list of strings
        lines = []
        for step in range(len(itemdata["positions"])):
            line = create_line(self, itemdata, step)
            lines.append(line)
        # Reverse the list if the scanning is from "end"
        if itemdata["positions"][-1] < itemdata["positions"][0]:
            lines = list(reversed(lines))
        # Otherwise the scan is from "start" or there is only one datapoint. Anyway, this means no need to reverse the data in these cases.

    except KeyError:
        # No data yet
        lines = []
        
    except IndexError:
        # End of program most likely
        lines = []

    return lines


def create_string_list_to_display(self):
    # Check which item is selected
    selected_itemtext = self.rawData_listWidget.currentItem().text()

    try:
        # Get the data related to the selected item
        itemdata = self.batch_data[selected_itemtext]

        # Convert the item data to list of strings
        lines = []
        # print(f"Selected item: {selected_itemtext}")
        # print("Steps available in create_line:", itemdata['absolute'][0])
        for step in range(len(itemdata["positions"])):
            line = create_line(self, itemdata, step)
            lines.append(line)
        # Reverse the list if the scanning is from "end"
        if itemdata["positions"][-1] < itemdata["positions"][0]:
            lines = list(reversed(lines))
        # Otherwise the scan is from "start" or there is only one datapoint. Anyway, this means no need to reverse the data in these cases.

    except KeyError:
        # No data yet
        lines = []

    return lines


def write_temp_file(self, current_run_no, current_step, string_list_to_write):
    """
    Creates and writes a temporary text file containing specified string data.

    The method generates a temporary file with a name based on the current date and time,
    as well as the provided run number. This file will only be created after the first
    execution step of a scan (when `current_step == 0`). The file is saved in a
    designated temporary folder within the base directory of the project.

    Parameters:
    ----------
    current_run_no : int
        The current run number of the scan. This is used to help format the temporary
        filename and ensure uniqueness across multiple runs.

    current_step : int
        The current step number in a scan process. If this is the first step (0), a new
        temporary file will be created.

    string_list_to_write : list of str
        A list of strings to be written to the temporary text file. Each string in the
        list will be written as a separate line in the file.

    Raises:
    ------
    OSError:
        If there is an issue with writing the file or creating the directory.

    Example:
    --------
    # Assuming this method is part of a class and instance is created:
    instance.write_temp_file(current_run_no=0, current_step=0, string_list_to_write=["First line", "Second line"])

    This will create a temporary file named with the current date and time, such as:
    "2023_10_18__14_30_45__run_1.txt" containing:
    First line
    Second line

    Note:
    -----
    The temporary file is stored in a "temporary" directory relative to the base
    directory of the project.
    """
    now = datetime.now()
    if current_step == 0:  # after the first current_step of a scan
        # Create the temp file
        self.cur_date = now.strftime("%Y_%m_%d")
        self.cur_time = now.strftime("%H_%M_%S")
        self.temp_fname = f"{self.cur_date}__{self.cur_time}__run_{current_run_no+1}.txt"

    base_folder_name = "python_project"
    current_path = os.path.abspath(os.path.dirname(__file__))
    current_path = get_base_directory(base_folder_name, current_path)
    fpath = os.path.join(current_path, "./temporary", self.temp_fname)

    with open(fpath, "w") as f:
        for line in string_list_to_write:
            f.write(line)


def write_raw_data_log(self, string_list_to_write: list):
    """
    Writes a list of strings to a text browser widget, updating the display
    only if the new content is different from the existing content.

    Parameters:
        string_list_to_write (list): A list of strings to be written
                                       to the text browser. Each string
                                       represents a line of text.

    Behavior:
        - The function retrieves the current text displayed in the text browser
          and compares it with the new content provided in `string_list_to_write`.
        - If the current content matches the new content exactly, the function
          exits early without making any changes.
        - If there is a difference, the function will:
            1. Retrieve the current vertical scrollbar position of the text browser.
            2. Clear the existing content of the text browser.
            3. Append each line from `string_list_to_write` to the text browser,
               ensuring that any trailing newline characters are removed.
            4. Restore the scrollbar position to its previous state after
               updating the content.

    Exception Handling:
        - The function includes a try-except block to catch any exceptions that
          may arise during the execution of the content update. If an exception
          occurs, it is logged with the traceback for debugging purposes.

    Returns:
        None
    """
    # Get the current content of the text browser in the incoming data format
    old_content = self.rawData_textBrowser.toPlainText().splitlines()

    # Compare old content with new incoming content
    if old_content == string_list_to_write:
        return  # Exit the function if they are the same

    try:
        # Get the text browser current scrollbar position
        scrollbar = self.rawData_textBrowser.verticalScrollBar()
        current_position = scrollbar.value()

        # Clear the text browser before writing
        self.rawData_textBrowser.clear()

        # Write new data in the browser
        self.rawData_textBrowser.setPlainText("".join(string_list_to_write))

        # Restore the scrollbar position
        scrollbar.setValue(current_position)

    except Exception:
        logging.error(traceback.format_exc())


def get_experiment_description_header(self, mean_values=False):
    # Full description header
    exp_descr = self.experimentDescription_plainTextEdit.toPlainText()
    header_1 = (
        "Z-scan Measurement\n"
        f"Code: {self.codeOfSample_comboBox.currentText()}\n"
        f"Silica thickness: {self.silicaThickness_dataSavingTab_doubleSpinBox.text()}\n"
        f"Concentration: {self.concentration_dataSavingTab_doubleSpinBox.text()}\n"
        f"Wavelength: {self.wavelength_dataSavingTab_doubleSpinBox.text()}\n"
    )
    if exp_descr != "":
        header_2 = "\nExperiment description:\n" f"{exp_descr}\n"
    else:
        header_2 = ""

    header_3 = f"\nThe data is the arithmetic mean of {self.numberOfScans_spinBox.value()} scans performed in one go.\n" if mean_values else "\n"

    header_4 = (
        "--------------------------------------------------------------------------------------------\n\n"
        f"Starting pos: {self.startPos_doubleSpinBox.value()}\n"
        f"Ending pos: {self.endPos_doubleSpinBox.value()}\n"
        "CH1:   Closed aperture\n"
        "CH2:   Reference\n"
        "CH3:   Open aperture\n"
        "CH4:   Empty channel\n\n"
        "--------------------------------------------------------------------------------------------\n\n"
        "SNo.   CA channel [V]        Ref channel [V]         OA channel [V]        Empty channel [V]\n\n"
        "--------------------------------------------------------------------------------------------\n"
    )

    return header_1 + header_2 + header_3 + header_4


def update_fullLog(self, update_button_clicked=True):
    """Runs at the end of experiment and on click of UPDATE button"""
    if update_button_clicked:
        self.toggle_datasaving_buttons_lock(True, add_buttons=[self.saveData_pushButton])  # unlock the data saving button

    if self.numberOfScans_spinBox.value() > 1:
        mean_values = True
    else:
        mean_values = False
    header = get_experiment_description_header(self, mean_values)

    try:
        assert self.rawData_listWidget.currentRow() == self.fullExp_listWidget.currentRow()
    except AssertionError:
        # Switch raw log item
        switch_widget_items(self, self.fullExp_listWidget.currentItem(), write_full_log=False)
    finally:
        # Get the content of fullLogData browser
        old_content = self.fullLogData_textBrowser.toPlainText().splitlines()

        # Get the new content for fullLogData browser based on currently selected item
        new_content = header.splitlines() + ["\n"] + self.rawData_textBrowser.toPlainText().splitlines()

    # If the new content is different than the previous one, display it
    if new_content != old_content:
        try:
            self.fullLogData_textBrowser.clear()  # Clear the text browser
            for line in new_content:
                self.fullLogData_textBrowser.append(line.rstrip("\n"))  # Add new lines
        except Exception:
            logging.error(traceback.format_exc())

    self.fullLogData_textBrowser.verticalScrollBar().setValue(0)

    if update_button_clicked:
        self.datalogTabs.setCurrentIndex(1)

    # Save mean values in data_set
    if self.data_acquisition_complete is True:
        self.data_set = (
            [i for i in range(self.stepsScan_spinBox.value() + 1)],
            self.batch_data["Mean values"]["absolute"][0],
            self.batch_data["Mean values"]["absolute"][1],
            self.batch_data["Mean values"]["absolute"][2],
        )  # , self.batch_data["Mean values"]["absolute"][3] not using the last column with zeros

    if update_button_clicked:  # and self.current_run == self.numberOfScans_spinBox.value():
        # FILENAMES
        now = datetime.now()
        self.cur_date = now.strftime("%Y_%m_%d")
        self.cur_time = now.strftime("%H_%M")
        sample_type = self.codeOfSample_comboBox.currentText()

        concentration = self.concentration_dataSavingTab_doubleSpinBox.value()
        conc_hyphen = str(f"{concentration:.2f}").replace(".", "-")

        wavelength = self.wavelength_dataSavingTab_doubleSpinBox.value()
        wavel_hyphen = str(f"{wavelength:.1f}").replace(".", "-")

        self.rawLogFilename_lineEdit.setText(f"{self.cur_date}__{self.cur_time}__{sample_type}_{conc_hyphen}_{wavel_hyphen}.txt")
        self.fullLogFilename_lineEdit.setText(f"{self.cur_date}__{self.cur_time}__{sample_type}_{conc_hyphen}_{wavel_hyphen}_2.txt")

        self.files = (self.rawLogFilename_lineEdit.text(), self.fullLogFilename_lineEdit.text())


def data_save(self):
    ## DATABASE CONNECTION
    # save to database before, so that if any error occurs in writing to file, the data can be retrieved from db.
    
    # sample_code = self.codeOfSample_comboBox.currentText()
    # sample_type = self.cuvetteType_comboBox.currentText()
    
    # match sample_type:
    #     case "Silica":
    #         params = {'silica_thicknessMM': self.silicaThickness_dataSavingTab_doubleSpinBox.value()}
    #     case "Solvent":
    #         params = {'solvent_density': 0,
    #                   'solvent_index': 0}
    #     case "Sample":
    #         params = {'sample_solvent': self.solventName_dataSavingTab_comboBox.currentText(),
    #                   'sample_concentration': self.concentration_dataSavingTab_doubleSpinBox.value()}
        # no other cases are available in GUI
        
    # def write_sample_to_db(): ...
    
    # write_sample_to_db(sample_code=sample_code,
    #                    sample_type=sample_type, **params)
    
    self.accurate_path = os.path.join(self.mainDirectory_lineEdit.text(), self.cur_date)
    
    try:
        os.makedirs(self.accurate_path)  # creates all folders that do not exist on the way to the last folder
    except FileExistsError:
        print(f"Directory already exists at: {self.accurate_path}")
    except Exception:
        logging.error(traceback.format_exc())

    try:
        for file in self.files:
            with open(os.path.join(self.accurate_path, file), "w") as f:
                if self.files.index(file) == 0:
                    f.write(self.rawData_textBrowser.toPlainText())
                else:
                    f.write(self.fullLogData_textBrowser.toPlainText())

        self.showdialog("Info", "Files successfully written!")

    except PermissionError:
        self.showdialog("Error", "Permission denied!\nCannot write the file in this directory.")
    except FileNotFoundError:
        self.showdialog("Error", "File not found.")
    except Exception:
        self.showdialog("Error", "Unknown error. Contact the software developer.")

    else:  # no error
        try:
            base_folder_name = "python_project"
            current_path = os.path.abspath(os.path.dirname(__file__))
            current_path = get_base_directory(base_folder_name, current_path)
            fpath = os.path.join(current_path, "./temporary")

            for filename in os.listdir(fpath):
                file_path = os.path.join(fpath, filename)
                if os.path.isfile(file_path):  # check if it's a file (and not a directory)
                    os.remove(file_path)
        
        except FileNotFoundError:
            logging.error(traceback.format_exc())
        except Exception:
            logging.error(traceback.format_exc())
        
        else:  # only after successful process (otherwise the user won't be able to save)
            self.saveData_pushButton.setEnabled(False)
            self.sendToFit_pushButton.setEnabled(True)
