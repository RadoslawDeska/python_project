from PyQt5.QtWidgets import QDoubleSpinBox
from PyQt5.QtCore import pyqtSignal, QEvent

from classes.mainwindow.managers import FocusManager

class MyDoubleSpinBox(QDoubleSpinBox):
    elementsReceived = pyqtSignal(object)
                                  
    def __init__(self, parent=None):
        super(MyDoubleSpinBox, self).__init__(parent)
        self.elementsReceived.connect(self.handle_elements)
        
        self.installEventFilter(self)
    
    def handle_elements(self, widgets_to_check):
        self.widgets_to_check = widgets_to_check

    def eventFilter(self, obj, event):
        if event.type() == QEvent.FocusIn:
            if obj in self.widgets_to_check:
                FocusManager()._previously_focused = obj
                print(f"SPINBOX: FocusInEvent {self}")
        
        return super(QDoubleSpinBox, self).eventFilter(obj, event)