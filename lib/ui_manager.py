"""UI Manager - Handles window initialization and styling.

Manages UI setup, styling, and layout configuration for the main window.
"""

from PyQt5 import QtGui
from PyQt5.QtWidgets import QMainWindow


class UIManager:
    """Manages UI initialization and styling for the Z-scan application."""
    
    def __init__(self, window: QMainWindow):
        """Initialize UI manager.
        
        Args:
            window: The main PyQt5 window
        """
        self.window = window
    
    def setup_initial_visibility(self):
        """Set initial visibility for UI elements."""
        self.window.solventOA_absorptionModel_label.setVisible(False)
        self.window.solventOA_absorptionModel_comboBox.setVisible(False)
        self.window.solventOA_fixROI_checkBox.setVisible(False)
        
        self.window.solventOA_saturationModel_label.setVisible(False)
        self.window.solventOA_saturationModel_comboBox.setVisible(False)
    
    def setup_stylesheet(self):
        """Apply custom stylesheet to the main window."""
        stylesheet = """
            QTabBar::tab:hover {background: rgba(100, 150, 255, 0.3)} 
            QTabBar::tab:!selected {margin-top: 3px;}
            QTabBar::tab:selected {height: 15px; border: 1px solid rgba(42, 130, 218, 1); border-top-left-radius: 3px; border-top-right-radius: 3px;
                background: rgba(42, 130, 218, 1); color: white}
            """
        self.window.setStyleSheet(stylesheet)
    
    def store_default_palette(self):
        """Store the default application palette."""
        self.window.default_palette = QtGui.QGuiApplication.palette()
        return self.window.default_palette
