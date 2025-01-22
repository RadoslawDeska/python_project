from PyQt5.QtWidgets import QApplication, QMainWindow, QMessageBox
import sys


class ExperimentStopped(Exception):
    pass

class GUI_Error(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("QErrorMessage Example")
        self.setGeometry(100, 100, 300, 200)
        self.center_on_screen()
        self.show_error_message()

    def show_error_message(self):
        error_message = QMessageBox(self)
        error_message.setIcon(QMessageBox.Critical)
        error_message.setWindowTitle("Error")
        error_message.setText("GUI file missing. Please reinstall the program.")
        error_message.setStandardButtons(QMessageBox.Ok)
        if error_message.exec_() == QMessageBox.Ok:
            sys.exit()
    
    def center_on_screen(self):
        # Get the screen geometry
        screen_geometry = QApplication.desktop().screenGeometry()
        
        # Calculate the center point
        self_center_x = (screen_geometry.width() - self.width()) // 2
        self_center_y = (screen_geometry.height() - self.height()) // 2

        # Move the widget to the center
        self.move(self_center_x, self_center_y)
