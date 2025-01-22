import sys
from PyQt5.QtWidgets import QApplication, QDialog, QVBoxLayout, QLabel, QLineEdit, QDialogButtonBox
from PyQt5.QtGui import QValidator
from PyQt5.QtCore import QLocale


class NoValuesBelowOneValidator(QValidator):
    def validate(self, input_str, pos):
        if not input_str:
            return QValidator.Intermediate, input_str, pos

        try:
            value = float(input_str)
            if value < 1 or ',' in input_str:
                return QValidator.Invalid, input_str, pos
            return QValidator.Acceptable, input_str, pos
        except ValueError:
            return QValidator.Invalid, input_str, pos

    def fixup(self, input_str):
        return input_str.replace(',', '')


class AboveZeroValidator(QValidator):
    def validate(self, input_str, pos):
        if not input_str:
            return QValidator.Intermediate, input_str, pos

        try:
            value = float(input_str)
            if value <= 0 or ',' in input_str:
                return QValidator.Invalid, input_str, pos
            return QValidator.Acceptable, input_str, pos
        except ValueError:
            return QValidator.Invalid, input_str, pos

    def fixup(self, input_str):
        return input_str.replace(',', '')


class SolventProperties(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.setWindowTitle("Set properties of the solvent")

        # Create layout
        layout = QVBoxLayout(self)
        
        # First input field
        self.label1 = QLabel()
        self.label1.setText('Solvent density (g/cm<sup>3</sup>)')
        self.lineEdit1 = QLineEdit(self)
        
        doubleValidator = AboveZeroValidator()
        doubleValidator.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        
        self.lineEdit1.setValidator(doubleValidator)
        
        layout.addWidget(self.label1)
        layout.addWidget(self.lineEdit1)

        # Second input field
        self.label2 = QLabel("Solvent refractive index:")
        self.lineEdit2 = QLineEdit(self)
        
        doubleValidator = NoValuesBelowOneValidator()
        doubleValidator.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        
        self.lineEdit2.setValidator(doubleValidator)
        
        layout.addWidget(self.label2)
        layout.addWidget(self.lineEdit2)

        # Dialog buttons
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self)
        layout.addWidget(self.buttonBox)

        # Connect signals
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)

    def get_values(self):
        return self.lineEdit1.text(), self.lineEdit2.text()

# Main application
if __name__ == "__main__":
    app = QApplication(sys.argv)
    dialog = SolventProperties()

    if dialog.exec_() == QDialog.Accepted:
        value1, value2 = dialog.get_values()
        print(f"First Value: {value1}")
        print(f"Second Value: {value2}")

    sys.exit(app.exec_())