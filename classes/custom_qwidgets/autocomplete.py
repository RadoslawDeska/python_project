__author__ = "Radosław Deska"

import json
import logging
import os
import traceback

from PyQt5.QtGui import QFocusEvent, QKeyEvent, QRegularExpressionValidator
from PyQt5.QtWidgets import QApplication, QComboBox, QDialog, QFileDialog, QListView, QStyleOptionComboBox, QStyle, QStyledItemDelegate
from PyQt5.QtCore import Qt, QEvent, pyqtSignal, QRegularExpression

from classes.mainwindow.managers import FocusManager, SampleManager as SM
from classes.custom_qwidgets.solvent_inputbox import SolventProperties
from lib.functions import get_base_directory

__all__ = ["get_data", "_select_JSON_file", "_save_JSON_file", "_populate_combo_with_params", "_load_item"]

class CustomListView(QListView):
    elementsReceived = pyqtSignal(object)
    focusOutSignal = pyqtSignal(object)
    
    def __init__(self, parent=None):
        super(CustomListView, self).__init__(parent)
        self.elementsReceived.connect(self.handle_elements)
        self.installEventFilter(self)
    
    def handle_elements(self, widgets_to_check):
        self.widgets_to_check = widgets_to_check

    def eventFilter(self, obj, event):
        if event.type() == QEvent.FocusIn:
            if obj in self.widgets_to_check:
                FocusManager()._previously_focused = obj
                print(f"LISTVIEW: FocusInEvent {self}")
        elif event.type() == QEvent.FocusOut:
            self.focusOutSignal.emit(self)
        
        return super(QListView, self).eventFilter(obj, event)
    
class MyComboBox(QComboBox):
    """The code was created with the help of Microsoft® Copilot AI Assistant (OpenAI)
    and DeepAI.org chat based on deep discussions about specific solutions
    """

    popupShow = pyqtSignal()
    popupHide = pyqtSignal()
    focusOutSignal = pyqtSignal()
    itemDeleted = pyqtSignal(object)
    elementsReceived = pyqtSignal(object)

    protected_items = ["Silica"]
    
    def __init__(self, parent=None):
        super(MyComboBox, self).__init__(parent)
        self.elementsReceived.connect(self.handle_elements)
        self._view = CustomListView(self)
        self.setView(self._view)
        self.setEditable(True)
        
        # GEOMETRY
        delegate = CustomDelegate(self)
        self.setItemDelegate(delegate)
        
        # EVENTS
        self.view().installEventFilter(self)
        self.focus_out_processed = False
        self.parent_ = parent
        
        # INPUT VALIDATION
        self.forbidden_symbols = ["\\", "/", ":", "*", "?", '"', "<", ">", "|"]
        self.replacement_symbols = ["\u29f5", "\u2215", "\ua789", "\u204e", "\uff1f", "\u201d", "\uff1c", "\uff1e", "\u23d0"]
        pattern = "^[^" + "".join(self.forbidden_symbols) + "]*$"
        regex = QRegularExpression(pattern)
        validator = QRegularExpressionValidator(regex, self)
        
        self.lineEdit().setValidator(validator)

    def handle_elements(self, widgets_to_check):
        self.widgets_to_check = widgets_to_check
        self._view.elementsReceived.emit(self.widgets_to_check)
    
    def eventFilter(self, obj, event):
        if event.type() == QEvent.KeyPress:
            if event.key() == Qt.Key_Delete:
                index = self.view().currentIndex().row()
                current_text = self.itemText(index)
                if current_text in self.protected_items:
                    print("Protected item. Cannot delete.")
                    return True # Ignore the delete event
                
                if index >= 0:
                    delete_item = self.itemText(index)
                    self.removeItem(index)
                    print(f"Item {delete_item} deleted.")
                    self.itemDeleted.emit(delete_item)
                
                return True
        
        elif event.type() == QEvent.Show and obj == self.view():
            self.popupShow.emit()
        elif event.type() == QEvent.Hide and obj == self.view():
            self.popupHide.emit()
        
        elif event.type() == QEvent.ShortcutOverride:
            if event.key() == Qt.Key_Tab:
                self.hidePopup()
                self.setCurrentIndex(self.view().currentIndex().row())
                self.focusNextPrevChild(True)
                self.focusOutSignal.emit()
        
        elif event.type() == QEvent.FocusOut and obj in {self, self.lineEdit()}:
            print("Event filter checks the focus_out_processed to: True")
            self.focus_out_processed = True
        elif event.type() == QEvent.FocusIn:
            print("Objects in self.widgets_to_check:")
            for o in self.widgets_to_check:
                print(o)
            print("")
            print(f"Currently FocusIn event is triggered for {obj}")
            if obj in self.widgets_to_check:
                FocusManager()._previously_focused = obj
                print(f"Focused in on {self.widgets_to_check[self.widgets_to_check.index(obj)]}")
            
            if obj in {self, self.lineEdit()}:
                print("Event filter sets the focus_out_processed to: False")
                self.focus_out_processed = False
        
        return super(MyComboBox, self).eventFilter(obj, event)

    def allItems(self) -> list:
        return [self.itemText(i) for i in range(self.count())]

    def allData(self) -> dict:
        return {self.itemText(i): self.itemData(i) for i in range(self.count())}

    def sort(self, col_num=0, order=Qt.AscendingOrder):
        """Remove this functionality completely"""
        self.model().sort(col_num, order)
        return
    
    def showPopup(self) -> None:
        # print("showPopup called") # Debug print
        # Calculate the available screen geometry
        screen_geometry = QApplication.desktop().availableGeometry(self)
        combo_geometry = self.geometry()
        global_pos = self.mapToGlobal(combo_geometry.bottomLeft())

        # Get the geometry of the dropdown list
        option = QStyleOptionComboBox()
        self.initStyleOption(option)
        popup_rect = self.style().subControlRect(QStyle.CC_ComboBox, option, QStyle.SC_ComboBoxListBoxPopup, self)
        
        # Calculate the position to ensure it shows below
        popup_rect.moveTo(global_pos)

        # Ensure the dropdown is entirely visible within the screen
        if popup_rect.bottom() > screen_geometry.bottom():
            popup_rect.moveBottom(screen_geometry.bottom())

        # Calculate the height of the dropdown based on the number of items
        item_height = self.view().sizeHintForRow(0)
        num_items = min(self.count(), self.maxVisibleItems())
        popup_height = item_height * num_items # Set the minimum height for the view to ensure it displays correctly
        self.view().setMinimumHeight(popup_height)
        
        self.view().setGeometry(popup_rect)
        super().showPopup()

    def hidePopup(self) -> None:
        self.popupHide.emit()
        super(MyComboBox, self).hidePopup()

    def setViewPosition(self, below):
        if below:
            # This is the default behavior; the dropdown will show below
            pass
        else:
            # Calculate the position to ensure it shows below
            option = QStyleOptionComboBox()
            self.initStyleOption(option)
            rect = self.style().subControlRect(QStyle.CC_ComboBox, option, QStyle.SC_ComboBoxListBoxPopup, self)
            rect.moveTo(self.mapToGlobal(self.rect().bottomLeft()))
            self.view().setGeometry(rect)
    
    def focusInEvent(self, e: QFocusEvent) -> None:
        print(f"COMBOBOX: FocusInEvent {self}")
        self.focus_out_processed = False
        super(MyComboBox, self).focusInEvent(e)
    
    def focusOutEvent(self, e: QFocusEvent | None) -> None:
        if e is not None:
            current_text = self.currentText()
            allitems_lower = [item.lower() for item in self.allItems()]

            # Check for valid current text
            if current_text.strip() and current_text.lower() not in allitems_lower:
                try:
                    insert_index = self._findInsertIndex(current_text.lower())
                    self.insertItem(insert_index, current_text)  # Assume `None` is okay for user data

                except Exception as ex:
                    print(f"Error inserting item: {ex}")  # Handle or log the error appropriately
        
        super(MyComboBox, self).focusOutEvent(e)

    def _findInsertIndex(self, new_text) -> int:
        """
        Finds the index at which new_text should be inserted to maintain
        sorted order (case-insensitive) in the list of current items.

        :param new_text: The text to be inserted.
        :return: The index at which new_text should be inserted.
        """
        for i in range(self.count()):
            if self.itemText(i).lower() > new_text.lower():
                return i

        return self.count()

    # def findText(self)
    
    def keyPressEvent(self, e: QKeyEvent | None) -> None:
        if e.key() in [Qt.Key_Enter, Qt.Key_Return]:
            # Change focus to the next GUI element
            # This will trigger focusOutEvent
            # and will add an item if it doesn't exist
            # in the items list
            self.focusNextPrevChild(True)
        
        # Get the current input text
        current_text = self.lineEdit().text()

        # Get the character that is pressed
        char = e.text()

        # Check if the character is invalid
        if char in self.forbidden_symbols:
            # Replace the invalid character with the replacement character
            new_text = current_text + self.replacement_symbols[self.forbidden_symbols.index(char)]
            self.lineEdit().setText(new_text)
            # Move the cursor to the end of the input
            self.lineEdit().setCursorPosition(len(new_text))
            return

        super(MyComboBox, self).keyPressEvent(e)
        
    def setCurrentText(self, text: str) -> None:
        if text:
            allitems_lower = [item.lower() for item in self.allItems()]
            if text.strip() and text.lower() not in allitems_lower:
                try:
                    insert_index = self._findInsertIndex(text.lower())
                    self.insertItem(insert_index, text)
                except Exception as ex:
                    print(f"Error inserting item: {ex}")
        
        super(MyComboBox, self).setCurrentText(text)
    
    def setCurrentTextWithData(self, text: str, data):
        # Add the item if it doesn't exist
        self.setCurrentText(text)
        
        # Find the index of the newly added or existing item
        index = self.findText(text, flags=Qt.MatchFixedString)
        if index != -1:
            # Add or update the associated data
            self.setItemData(index, data)
        else:
            print(f"Error: Item '{text}' not found in combo box.")

class MyComboBox_NonEditable(MyComboBox):
    def __init__(self, parent=None):
        super(MyComboBox_NonEditable, self).__init__(parent)
        
        self.setEditable(False)
    
class CustomDelegate(QStyledItemDelegate):
    def sizeHint(self, option, index):
        size = super().sizeHint(option, index)
        # Adjust the height if necessary
        size.setHeight(20)  # Example height, adjust as needed
        return size

def get_data(self, box):
    match box:
        case self.codeOfSample_comboBox:
            name = self.codeOfSample_comboBox.currentText()
            match self.cuvetteType_comboBox.currentText():
                case "silica":
                    thickness = self.silicaThickness_dataSavingTab_doubleSpinBox.value()
                    data = (name, {"thickness": thickness})
                
                case "solvent":
                    if name in self.solventName_dataSavingTab_comboBox.allItems():
                        data = self.solventName_dataSavingTab_comboBox.allData()[name]
                    else:
                        dialog = SolventProperties()
                        if dialog.exec_() == QDialog.Accepted:
                            density, index = dialog.get_values()
                            if not density:
                                density = 0
                            if not index:
                                index = 0
                        else:
                            density, index = 0, 0
                        
                        data = (name, {"density": float(density), "index": float(index)})
                    
                case "sample":
                    concentration = self.concentration_dataSavingTab_doubleSpinBox.value()
                    solvent = self.solventName_dataSavingTab_comboBox.currentText()
                    data = (name, {"concentration": concentration, "solvent": solvent})
        
        case self.solventName_comboBox:
            name = self.solventName_comboBox.currentText()
            density = self.solventDensity_doubleSpinBox.value()
            rind = self.solventRefrIdx_doubleSpinBox.value()
            data = (name, {"density": density, "index": rind})
        
        case self.solventName_dataSavingTab_comboBox:
            name = self.solventName_dataSavingTab_comboBox.currentText()
            data = (name, {"density": 0, "index": 0})
    
    return data


def _select_JSON_file(self, caller="Startup"):
    """This method is used for selecting the JSON file containing samples/solvents data. It also takes care of non-existing file in settings
    and of dialog box closing without selecting a file. If all is well, the selected JSON file is passed on for data extraction.

    Args:
        `caller` (str, optional): This tells what triggered this method (either the program startup ("Startup") or menu action ("ActionLoadSolvents"
        or "ActionLoadSamples"). Defaults to "Startup".
    """
    base_folder_name = "python_project"
    current_path = os.path.abspath(os.path.dirname(__file__))
            
    match caller:
        ## Startup action
        case "Startup":
            # Load solvents file from settings
            try:
                with open(self.settings.value("FittingTab/solvents_path"), mode="r", encoding="utf-8") as json_file:
                    new_data = json.load(json_file)

            # or from user-selected file, if not found
            except FileNotFoundError:
                self.showdialog("Warning", "solvents.json not found in the default location. Select the file.")
                path = os.path.abspath(os.path.dirname(__file__))  # this is where solvents.json is expected to be
                file, _ = QFileDialog.getOpenFileName(self, "Open File", path, filter="JSON file (*.json)")

                if file:  # if dialog was not cancelled
                    with open(file, mode="r", encoding="utf-8") as json_file:
                        new_data = json.load(json_file)
                else:
                    new_data = None

            # or set new_data as None if other errors occur
            except Exception:
                logging.error(traceback.format_exc())
                print("Setting new_data as None")
                new_data = None

            # if no error
            else:
                box = self.solventName_comboBox

        ## Menu actions
        case "actionLoadSolvents":
            try:
                current_path = get_base_directory(base_folder_name, current_path)
                    
                file, _ = QFileDialog.getOpenFileName(self, "Open Solvents File", current_path, filter="JSON file (*.json)")
                if file:  # if dialog was not cancelled
                    with open(file, mode="r", encoding="utf-8") as json_file:
                        new_data = json.load(json_file)
                else:
                    new_data = None

            # or set new_data as None if errors occur
            except Exception:
                logging.error(traceback.format_exc())
                print("Setting new_data as None")
                new_data = None

            # # if no error
            # else:
            #     box = self.solventName_comboBox

        case "actionLoadSamples":
            try:
                current_path = get_base_directory(base_folder_name, current_path)
                
                file, _ = QFileDialog.getOpenFileName(self, "Open Samples File", current_path, filter="JSON file (*.json)")
                if file:  # if dialog was not cancelled
                    with open(file, mode="r", encoding="utf-8") as json_file:
                        new_data = json.load(json_file)
                else:
                    new_data = None

            # or set new_data as None if errors occur
            except Exception:
                logging.error(traceback.format_exc())
                print("Setting new_data as None")
                new_data = None

            # if no error
            else:
                box = self.codeOfSample_comboBox

    if new_data:
        if caller == "actionLoadSolvents" or caller == "Startup":
            _populate_combo_with_params(self, self.solventName_comboBox, caller, from_file=True, new_data=new_data)
            _populate_combo_with_params(self, self.solventName_dataSavingTab_comboBox, caller, from_file=True, new_data=new_data)
        else:
            _populate_combo_with_params(self, box, caller, from_file=True, new_data=new_data)


def _save_JSON_file(self, caller=None):
    base_folder_name = "python_project"
    current_path = os.path.abspath(os.path.dirname(__file__))
    
    match caller:
        case "actionSaveSamples":
            combo_data = self.codeOfSample_comboBox.allData()
            saving = "Samples"
        case "actionSaveSolvents":
            combo_data = self.solventName_comboBox.allData()
            saving = "Solvents"
        case _:
            combo_data = None

    if not combo_data:
        self.showdialog(msg_type="Warning", message="No items found.")
        return

    to_save = {}
    for _, lst in combo_data.items():
        k, v = lst
        to_save[k] = v

    try:
        current_path = get_base_directory(base_folder_name, current_path)
                
        file, _ = QFileDialog.getSaveFileName(self, f"Save {saving} File", current_path, filter="JSON file (*.json)")
        if file:  # if dialog was not cancelled
            with open(file, mode="w", encoding="utf-8") as json_file:
                json.dump(to_save, json_file, indent=4)  # Use 'indent' for pretty formatting

    except Exception:
        logging.error(traceback.format_exc())


def _populate_combo_with_params(self, box: MyComboBox, caller=None, from_file: bool = False, new_data=None):
    """Checks for duplicates between existing `box` items and `new_data` and populates `box` with selected items
    or starts anew. Assigns itemData from `new_data` to each item of `box` after the check.

    The function is called either by `_select_JSON_file` or upon editingFinished signal from proper box's parameters
    carrying elements.
    
    ## ISSUES:
    this most likely will result in double call of dialog box for duplicate entries (from_file=True) because the method is called twice:
    for two comboboxes with solvents (one after another in select_JSON_file method)
    """
    try:
        box.currentIndexChanged.disconnect()
    except Exception:
        pass
    
    print("[0]: _populate_combo_with_params method has: box.currentText =", box.currentText())
    
    match box:
        case self.codeOfSample_comboBox:
            data_type = "samples"
        case self.solventName_comboBox:
            data_type = "solvents"

    if from_file:
        old_data = box.allData()
        # check if new data exist in already existing items list
        for item in new_data.keys():
            if item in box.allData().keys():
                fresh, keep_new = self.showdialog(
                    msg_type="Combobox",
                    message=f"Some {data_type} already exist in the old list. What do you wish to do?",
                    buttons=("fresh", "new_dup", "old_dup"),
                )
                if fresh:
                    combodata = new_data  # replaces the old list
                    break  # DUPLICATE FOUND; no more is important
                if keep_new:
                    combodata = dict(old_data.values()) | new_data  # keeps new duplicates
                else:
                    combodata = new_data | dict(old_data.values())  # keeps old duplicates
                break  # DUPLICATE FOUND; no more is important
        ######## DON'T INDENT THE ELSE BELOW (for-else loop) ########
        else:  # NO DUPLICATES FOUND; doesn't execute if break occured
            if caller != "Startup" and old_data:
                fresh, add_ = self.showdialog(msg_type="Combobox", message="What do you wish to do?", buttons=("fresh", "append"))
            else:
                fresh = True
            if fresh:
                combodata = new_data  # replaces the old list
            elif not fresh and add_:
                # print('here')
                combodata = old_data | new_data

        # Fill the names/codes from json
        box.clear()
        for i, (k, v) in enumerate(zip(combodata.keys(), combodata.values())):
            box.insertItem(i, k, v)
        # box.sort(0, Qt.AscendingOrder)  # DON'T SORT, IT BREAKS CASE-INSENSITIVE INSERTING USING BOX._findInsertIndex
        box.setCurrentIndex(0)

        # Fill the rest of data for solvents/samples
        for index, item in enumerate(box.allData().items()):
            box.setItemData(index, item)

        # print(box.allData())

    else:  # from_file = False; if user edits the parameters (editingFinished signal)
        print("getting data")
        cb_index = box.currentIndex()
        print(f"{cb_index=}")
        res = get_data(self, box)  # opens dialog box for solvent parameters
        box.setItemData(cb_index, res)
        box.setCurrentText(res[0])  # prevent item insertion being reverted to the previous value upon the first insertion
        print("[after setItemData]: _populate_combo_with_params method has: box.currentText =", box.currentText())
        if box.itemData(cb_index):
            print('data written')
            print(box.itemData(cb_index))
        else:
            print('no data!!')
        # box.sort(0, Qt.AscendingOrder)  # DON'T SORT, IT BREAKS CASE-INSENSITIVE INSERTING USING BOX._findInsertIndex

    print("[4]: _populate_combo_with_params method has: box.currentText =", box.currentText())
    
    box.currentIndexChanged.connect(lambda: _load_item(self, box))
    _load_item(self, box)
    print('End of populate_combo_with_params') 


def _load_item(self, box: QComboBox) -> None:
    """Loads itemData to elements displaying the parameters.\\        
        The function is called by programmatic changes in the `box` itself (`_populate_combo_with_params` method)
        and through emission of `currentIndexChanged` signal on the box.
        
        Args:
            `box` (QtWidgets.QComboBox) : should be the combo box containing the names/codes to which further
            the data is related to.
        """
    curindex = box.currentIndex()
    current_code = self.codeOfSample_comboBox.currentText()
    
    if curindex == -1:  # No selection
        return
    
    match self.cuvetteType_comboBox.currentText():
        case "silica":
            pass
        case "solvent":
            if self.codeOfSample_comboBox.currentText() == self.solventName_dataSavingTab_comboBox.currentText():  # ????
                return
        case "sample":
            pass
        
    # Item selected    
    match box:
        case self.solventName_comboBox:  # Data fitting tab
            if not self.solventName_comboBox.currentText():
                return
            self.solventDensity_doubleSpinBox.setValue(box.itemData(curindex)[1]["density"])
            self.solventRefrIdx_doubleSpinBox.setValue(box.itemData(curindex)[1]["index"])
        
        case self.codeOfSample_comboBox:  # Data saving tab
            if not current_code or current_code not in SM().allsamples_db:
                return
                        
            # data = get_data(self, box)
            cuvette = self.cuvetteType_comboBox.currentText()
            match cuvette:
                case "silica":
                    thickness_box = self.silicaThickness_dataSavingTab_doubleSpinBox
                    thickness_box.setValue(float(SM().allsamples_db[current_code][1]['silica_thicknessMM']))
                
                case "solvent":
                    solvent_box = self.solventName_dataSavingTab_comboBox
                    # solvent_box.setCurrentText(data[0])
                    # solvent_box.setItemData(solvent_box.currentIndex(), data)
                
                case "sample":
                    conc_box = self.concentration_dataSavingTab_doubleSpinBox
                    solvent_box = self.solventName_dataSavingTab_comboBox
                    
                    conc_box.setValue(float(SM().allsamples_db[current_code][1]['sample_concentration']))
                    solvent_box.setCurrentText(SM().allsamples_db[current_code][1]['sample_solvent'])
                    
                