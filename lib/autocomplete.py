"""The code was created with the help of Microsoft® Copilot AI Assistant (OpenAI)
and DeepAI.org chat based on deep discussions about specific solutions
"""

__author__ = "Radosław Deska"

import json
import logging
import os
import traceback

from PyQt5.QtGui import QFocusEvent, QKeyEvent, QRegularExpressionValidator
from PyQt5.QtWidgets import QComboBox, QFileDialog
from PyQt5.QtCore import Qt, QEvent, pyqtSignal, QRegularExpression

class MyComboBox(QComboBox):
    popupShow = pyqtSignal()
    popupHide = pyqtSignal()
    
    def __init__(self, parent=None):
        super(MyComboBox, self).__init__(parent)
        self.setEditable(True)
        self.view().installEventFilter(self)
        
        self.forbidden_symbols = ["\\","/",":","*","?","\"","<",">","|"]
        self.replacement_symbols = ["\u29F5","\u2215","\uA789","\u204E","\uFF1F","\u201D","\uFF1C","\uFF1E","\u23D0"]
        pattern = "^[^" + "".join(self.forbidden_symbols) + "]*$"
        regex = QRegularExpression(pattern)
        validator = QRegularExpressionValidator(regex, self)
        self.lineEdit().setValidator(validator)

    def eventFilter(self, obj, event):
        if event.type() == QEvent.KeyPress:
            if event.key() == Qt.Key_Delete:
                index = self.view().currentIndex().row()
                if index >= 0:
                    self.removeItem(index)
                return True
        return super(MyComboBox, self).eventFilter(obj, event)
    
    def allItems(self) -> list:
        return [self.itemText(i) for i in range(self.count())]
    
    def allData(self) -> dict:
        return {self.itemText(i): self.itemData(i) for i in range(self.count())}
    
    def sort(self, col_num=0, order=Qt.AscendingOrder):
        self.model().sort(col_num, order)
    
    def showPopup(self) -> None:
        self.popupShow.emit()
        super(MyComboBox, self).showPopup()
    
    def hidePopup(self) -> None:
        self.popupHide.emit()
        super(MyComboBox, self).hidePopup()
    
    def focusOutEvent(self, e: QFocusEvent | None) -> None:
        current_text = self.currentText()
        allitems_lower = [it.lower() for it in self.allItems()]
        
        if current_text.strip().lower() and current_text.lower() not in allitems_lower:
            insert_index = self.findInsertIndex(current_text.lower())
            self.insertItem(insert_index, current_text, None)
            # self.sort(0, Qt.AscendingOrder)

        super(MyComboBox, self).focusOutEvent(e)
    
    def findInsertIndex(self, new_text):
        for i in range(self.count()):
            if self.itemText(i).lower() > new_text.lower():
                return i

        return self.count()
    
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

class ComboBoxParametrization():
    '''The class is inherited by the parent Window class, so it uses the parent's elements without need of passing them as arguments.'''
    def get_data(self, box):
        match box:
            case self.codeOfSample_comboBox:
                name = self.codeOfSample_comboBox.currentText()
                concentration = self.concentration_dataSavingTab_doubleSpinBox.value()
                data = (name, {"concentration": concentration})
            case self.solventName_comboBox:
                name = self.solventName_comboBox.currentText()
                density = self.solventDensity_doubleSpinBox.value()
                rind = self.solventRefrIdx_doubleSpinBox.value()
                data = (name, {"density": density, "index": rind})
        return data
    
    def _select_JSON_file(self, caller="Startup"):
        '''This method is used for selecting the JSON file containing samples/solvents data. It also takes care of non-existing file in settings
        and of dialog box closing without selecting a file. If all is well, the selected JSON file is passed on for data extraction.
        
        Args:
            `caller` (str, optional): This tells what triggered this method (either the program startup ("Startup") or menu action ("ActionLoadSolvents"
            or "ActionLoadSamples"). Defaults to "Startup".
        '''

        match caller:
            ## Startup action
            case "Startup":
                # Load solvents file from settings
                try:
                    with open(self.settings.value('FittingTab/solvents_path'), mode="r", encoding="utf-8") as json_file:
                        new_data = json.load(json_file)
                
                # or from user-selected file, if not found
                except FileNotFoundError:
                    self.showdialog('Warning','solvents.json not found in the default location. Select the file.')
                    path = os.path.abspath(os.path.dirname(__file__)) # this is where solvents.json is expected to be
                    file, _ = QFileDialog.getOpenFileName(self, "Open File", path,filter="JSON file (*.json)")
                    
                    if file: # if dialog was not cancelled
                        with open(file, mode="r", encoding="utf-8") as json_file:
                            new_data = json.load(json_file)
                    else:
                        new_data = None
                
                # or set new_data as None if other errors occur
                except Exception:
                    logging.error(traceback.format_exc())
                    print('Setting new_data as None')
                    new_data = None
                
                # if no error    
                else:
                    box = self.solventName_comboBox

            ## Menu actions
            case "actionLoadSolvents":
                try:
                    path = os.path.abspath(os.path.join(os.path.dirname(__file__),'..'))
                    file, _ = QFileDialog.getOpenFileName(self, "Open Solvents File", path, filter="JSON file (*.json)")
                    if file: # if dialog was not cancelled
                        with open(file, mode="r", encoding="utf-8") as json_file:
                            new_data = json.load(json_file)
                    else:
                        new_data = None
                
                # or set new_data as None if errors occur
                except Exception:
                    logging.error(traceback.format_exc())
                    print('Setting new_data as None')
                    new_data = None
                    
                # if no error    
                else:
                    box = self.solventName_comboBox
            
            case "actionLoadSamples":
                try:
                    path = os.path.abspath(os.path.join(os.path.dirname(__file__),'..'))
                    file, _ = QFileDialog.getOpenFileName(self, "Open Samples File", path, filter="JSON file (*.json)")
                    if file: # if dialog was not cancelled
                        with open(file, mode="r", encoding="utf-8") as json_file:
                            new_data = json.load(json_file)
                    else:
                        new_data = None
                
                # or set new_data as None if errors occur
                except Exception:
                    logging.error(traceback.format_exc())
                    print('Setting new_data as None')
                    new_data = None
                
                # if no error    
                else:
                    box = self.codeOfSample_comboBox
        
        if new_data:
            self._populate_combo_with_params(box, caller, from_file=True, new_data=new_data)
        else:
            print("No data.")
            pass
    
    def _save_JSON_file(self, caller=None):
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
            self.showdialog(msg_type='Warning', 
                            message='No items found.')
            return
        
        to_save = {}
        for _, lst in combo_data.items():
            k, v = lst
            to_save[k] = v
        
        try:
            path = os.path.abspath(os.path.join(os.path.dirname(__file__),'..'))
            file, _ = QFileDialog.getSaveFileName(self, f"Save {saving} File", path, filter="JSON file (*.json)")
            if file: # if dialog was not cancelled
                with open(file, mode="w", encoding="utf-8") as json_file:
                    json.dump(to_save, json_file, indent=4)  # Use 'indent' for pretty formatting
        
        except Exception:
            logging.error(traceback.format_exc())

    def _populate_combo_with_params(self, box: MyComboBox, caller=None, from_file: bool = False, new_data=None):
        '''Checks for duplicates between existing `box` items and `new_data` and populates `box` with selected items 
        or starts anew. Assigns itemData from `new_data` to each item of `box` after the check.
        
        The function is called either by `_select_JSON_file` or upon editingFinished signal from proper box's parameters
        carrying elements.'''
        if not isinstance(box, MyComboBox):
            print("Provided combo box is not instance of MyComboBox.")
            return

        try:
            box.currentIndexChanged.disconnect()
        except Exception:
            pass
        
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
                    fresh, keep_new = self.showdialog(msg_type='Combobox', 
                                                    message=f'Some {data_type} already exist in the old list. What do you wish to do?',
                                                    buttons=('fresh','new_dup','old_dup'))
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
                    fresh, add_ = self.showdialog(msg_type='Combobox', 
                                                        message='What do you wish to do?',
                                                        buttons=('fresh','append'))
                else:
                    fresh = True
                if fresh:
                    combodata = new_data  # replaces the old list
                elif not fresh and add_:
                    # print('here')
                    combodata = old_data | new_data
            
            # Fill the names/codes from json
            box.clear()
            for i,(k,v) in enumerate(zip(combodata.keys(),combodata.values())):
                box.insertItem(i,k,v)
            # box.sort(0, Qt.AscendingOrder)  # DON'T SORT, IT BREAKS CASE-INSENSITIVE INSERTING USING BOX.findInsertIndex
            box.setCurrentIndex(0)
            
            # Fill the rest of data for solvents/samples
            for index, item in enumerate(box.allData().items()):
                box.setItemData(index, item)
            
            # print(box.allData())
        
        else:  # from_file = False; if user edits the parameters (editingFinished signal)
            cb_index = box.currentIndex()
            box.setItemData(cb_index, self.get_data(box))
            # box.sort(0, Qt.AscendingOrder)  # DON'T SORT, IT BREAKS CASE-INSENSITIVE INSERTING USING BOX.findInsertIndex
        
        box.currentIndexChanged.connect(lambda: self._load_item(box))
        self._load_item(box)
    
    def _load_item(self, box: QComboBox):
        '''Loads itemData to elements carrying the parameters.\\        
        The function is called by programmatic changes in the `box` itself (`_populate_combo_with_params` method)
        and through emission of `currentIndexChanged` signal on the box.
        
        Args:
            `box` (QtWidgets.QComboBox) : should be the combo box containing the names/codes to which further
            the data is related to.
        '''
        
        curindex = box.currentIndex()
        if curindex >= 0: # if list is not emptied
            if box.itemData(curindex):
                itemdata = box.itemData(curindex)[1]

                match box:
                    case self.codeOfSample_comboBox:
                        self.concentration_dataSavingTab_doubleSpinBox.setValue(itemdata['concentration'])
                    case self.solventName_comboBox:
                        self.solventDensity_doubleSpinBox.setValue(itemdata['density'])
                        self.solventRefrIdx_doubleSpinBox.setValue(itemdata['index'])
            else:
                self._populate_combo_with_params(box)
        else:
            pass