import sys
import os # Add the project root directory to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    
from classes.mainwindow.managers import FocusManager, SampleManager
from classes.custom_qwidgets.autocomplete import MyComboBox, MyComboBox_NonEditable
from classes.custom_qwidgets.spinboxes import MyDoubleSpinBox

import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QLabel, QFormLayout


class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()

        central_widget = QWidget(self)
        self.setCentralWidget(central_widget)

        layout = QFormLayout(central_widget)

        self.cuvetteType_comboBox = MyComboBox()
        self.codeOfSample_comboBox = MyComboBox()
        self.solventName_dataSavingTab_comboBox = MyComboBox_NonEditable()
        self.silicaThickness_dataSavingTab_doubleSpinBox = MyDoubleSpinBox()
        self.concentration_dataSavingTab_doubleSpinBox = MyDoubleSpinBox()

        layout.addRow(QLabel("Cuvette type:"), self.cuvetteType_comboBox)
        layout.addRow(QLabel("Sample code:"), self.codeOfSample_comboBox)
        layout.addRow(QLabel("Silica thickness:"), self.silicaThickness_dataSavingTab_doubleSpinBox)
        layout.addRow(QLabel("Concentration:"), self.concentration_dataSavingTab_doubleSpinBox)

        self.setWindowTitle("Custom Widgets Application")
        self.setGeometry(100, 100, 400, 200)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
