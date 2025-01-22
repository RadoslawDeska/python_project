from classes.custom_qwidgets.solvent_inputbox import SolventProperties
from PyQt5.QtWidgets import QDialog
from PyQt5.QtCore import pyqtSignal, QObject, Qt


class SampleManager():
    """
    SampleManager is a class that manages the registration, unregistration, and parameter extraction of samples.

    Attributes:
        allsamples_db (dict): A dictionary to store all registered samples.
        allowed_sample_types (set): A set of allowed sample types.

    Methods:
        __init__():
            Initializes the SampleManager instance.
        _ask_for_solvent():
            Opens a dialog to ask for solvent properties and returns the values.
        _write_solvent(parent_window, data):
            Updates the solvent selection in the parent window's data saving tab.
        register_sample(sample_name, sample_type, **kwargs):
            Registers a sample with the given name and type, along with additional parameters.
        unregister_sample(sample_name):
            Unregisters a sample from the allsamples_db.
        read_parameters(parent_window):
            Extracts parameters from the parent window's UI elements and returns them as a dictionary.
    """
    allsamples_db = {}
    allowed_sample_types = {"silica", "solvent", "sample"}
    
    def __init__(self):
        pass
    
    
    def _ask_for_solvent(self):
        """
        Opens a dialog to ask for solvent properties and returns the values.
        This method creates an instance of the SolventProperties dialog and 
        executes it. If the dialog is accepted, it retrieves the density and 
        index values from the dialog. If the values are not provided, they 
        default to 0. If the dialog is not accepted, both values default to 0.
        Returns:
            dict: A dictionary containing the solvent properties with keys:
                - "density" (float): The density of the solvent.
                - "index" (float): The index of the solvent.
        """
        dialog = SolventProperties()
        if dialog.exec_() == QDialog.Accepted:
            density, index = dialog.get_values()
            if not density:
                density = 0
            if not index:
                index = 0
        else:
            density, index = 0, 0
            
        return {"density": float(density), "index": float(index)}
    
    
    def _write_solvent(self, parent_window, data):
        """
        Updates the solvent selection in the parent window's data saving tab.

        Args:
            parent_window: The parent window object containing the solventName_dataSavingTab_comboBox.
            data (list): A list where the first element is the solvent name and the rest is additional data.

        Returns:
            None
        """
        solvents = parent_window.solventName_dataSavingTab_comboBox
        solvents.setCurrentTextWithData(data[0], data)
    
    
    def register_sample(self, sample_name, sample_type, **kwargs):
        """
        Registers a sample with the given name and type, along with additional parameters.
        ## Parameters:
        sample_name (str): The name of the sample to be registered.
        sample_type (str): The type of the sample. Must be one of the allowed sample types.
        **kwargs: Additional parameters required for the specific sample type.
        ### Sample Types and Required Parameters:
        - "silica": Requires "silica_thicknessMM".
        - "solvent": Requires "solvent_density" and "solvent_index".
        - "sample": Requires "sample_concentration" and "sample_solvent".
        ## Returns:
        None
        ## Prints:
        - A message indicating if the sample type is forbidden.
        - A message listing any missing parameters.
        - A message listing any parameters with incorrect types.
        - A success message if the sample is registered successfully.
        """
        if sample_type.lower() not in self.allowed_sample_types:
            print(f'Forbidden sample type: {sample_type}. Allowed values:', ', '.join(self.allowed_sample_types))
            return
    
        parameters = []
        missing = []
        wrong_type = []
        match sample_type.lower():
            case "silica":
                if "silica_thicknessMM" in kwargs:
                    parameters.append("silica_thicknessMM")
                else:
                    missing.append("silica_thicknessMM")
            case "solvent":
                if "solvent_density" in kwargs:
                    parameters.append("solvent_density")
                else:
                    missing.append("solvent_density")
                if "solvent_index" in kwargs:
                    parameters.append("solvent_index")
                else:
                    missing.append("solvent_index")
            case "sample":
                if "sample_concentration" in kwargs:
                    parameters.append("sample_concentration")
                else:
                    missing.append("sample_concentration")
                if "sample_solvent" in kwargs:
                    parameters.append("sample_solvent")
                else:
                    missing.append("sample_solvent")
        
        if missing:
            print('Missing parameters:', ', '.join(missing))
            return
        
        for param in parameters:
            if param != "sample_solvent":
                try:
                    float(kwargs[param])
                except (TypeError, ValueError):
                    wrong_type.append(param)
        
        if wrong_type:
            print("TypeError occurred for:", ', '.join(wrong_type))
            return
        
        self.allsamples_db[sample_name] = (
            sample_type.lower(), {param: kwargs[param] for param in kwargs} )
        
        print(f'Sample "{sample_name}" registered successfully as {sample_type.lower()}.')
    
    
    def unregister_sample(self, sample_name):
        """
        Unregisters a sample from the allsamples_db.

        This method removes the sample with the given name from the allsamples_db
        dictionary. If the sample does not exist, no action is taken.

        Args:
            sample_name (str): The name of the sample to unregister.

        Returns:
            None
        """
        self.allsamples_db.pop(sample_name, None)
        print(f'Sample "{sample_name}" unregistered successfully.')
    
    
    def read_parameters(self, parent_window):
        """
        Extracts parameters from the parent window's UI elements and returns them as a dictionary.
        Args:
            parent_window (QMainWindow): The main window containing the UI elements.
        Returns:
            dict: A dictionary containing the following keys:
            - "silica_thicknessMM" (float): The thickness of the silica in millimeters.
            - "solvent_density" (float): The density of the solvent. Defaults to 0 if not a solvent.
            - "solvent_index" (float): The index of the solvent. Defaults to 0 if not a solvent.
            - "sample_concentration" (float): The concentration of the sample.
            - "sample_solvent" (str): The name of the solvent.
        """
        cuvette = parent_window.cuvetteType_comboBox.currentText()
        name = parent_window.codeOfSample_comboBox.currentText()
        thickness = parent_window.silicaThickness_dataSavingTab_doubleSpinBox.value()
        concentration = parent_window.concentration_dataSavingTab_doubleSpinBox.value()
        solvent = parent_window.solventName_dataSavingTab_comboBox.currentText()
        
        if cuvette == "solvent":
            solvents = parent_window.solventName_dataSavingTab_comboBox
            item_number = solvents.findText(name, flags=Qt.MatchFixedString)
            if  item_number != -1:
                print(f"Item name entered is: {name}")
                print(f"Item name found is: {solvents.itemText(item_number)}")
                sv_density, sv_index = solvents.allData()[solvents.itemText(item_number)][1].values()
            else:
                user_data_dict = self._ask_for_solvent()
                self._write_solvent(parent_window, (name, user_data_dict))
                
                sv_density, sv_index = user_data_dict.values()
        else:
            sv_density, sv_index = 0, 0
        
        params = {
            "silica_thicknessMM": thickness,
            "solvent_density": sv_density,
            "solvent_index": sv_index,
            "sample_concentration": concentration,
            "sample_solvent": solvent
        }
        
        return params
    
    
class FocusManager(QObject):
    focusedChanged = pyqtSignal(QObject)
    _previously_focused = None
    
    def __init__(self):
        super().__init__()
        
    @property
    def previously_focused(self):
        return self._previously_focused
    
    @previously_focused.setter 
    def previously_focused(self, value):
        self._previously_focused = value
        self.focusedChanged.emit(value)
    
    