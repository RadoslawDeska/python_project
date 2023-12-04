from PyQt5.QAxContainer import QAxWidget
from PyQt5.QtCore import QVariant

class MG17Motor(QAxWidget):
    def __init__(self):
        super(MG17Motor, self).__init__()
        self.setControl('{3CE35BF3-1E13-4D2C-8C0B-DEF6314420B3}') # Thorlabs MGMotor Control
    
    def configure(self, hw_serial_num):
        self.dynamicCall('SetHWSerialNum(int)', QVariant(hw_serial_num))
        self.dynamicCall('StartCtrl()')