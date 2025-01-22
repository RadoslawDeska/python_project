from PyQt5.QtWidgets import QMessageBox

def showdialog(self, msg_type: str, message: str, **kwargs):
    '''
    Message type (msg_type) is one of these: 'Error', 'Warning', 'Info', 'Combobox'
    '''
    dialog_args = (msg_type, message)
    match msg_type:
        case "Error":
            button = QMessageBox.critical(self, *dialog_args)
        case "Warning":
            button = QMessageBox.warning(self, *dialog_args)
        case "Info":
            button = QMessageBox.information(self, *dialog_args)
        case "Question":
            button = QMessageBox.question(self, *dialog_args)                    
        case "Combobox":  # used for loading solvent/sample list from file
            box = QMessageBox()
            box.setIcon(QMessageBox.Question)
            box.setWindowTitle('Combobox items loading')
            box.setText(message)
            
            button_options = {
                'fresh': ('Start fresh list', QMessageBox.YesRole),
                'new_dup': ('Keep new duplicates', QMessageBox.YesRole),
                'old_dup': ('Keep older duplicates', QMessageBox.YesRole),
                'append': ('Append to old list', QMessageBox.YesRole),
            }
            
            buttons_added = {}
            if 'buttons' in kwargs:
                for button in kwargs['buttons']:
                    buttons_added[button] = box.addButton(*button_options[button])
                for button in button_options:
                    if button not in kwargs["buttons"]:  # all buttons that were not passed to the method
                        buttons_added[button] = None
                        
            buttons_added['cancel'] = box.addButton('Don\'t load anything', QMessageBox.RejectRole)
            box.setDefaultButton(buttons_added['cancel'])
            
            box.exec_()
            response = box.clickedButton()
            
            match response:
                case response if response is buttons_added['fresh']:
                    return (True,False)
                case response if response is buttons_added['new_dup']:
                    return (False,True)
                case response if response is buttons_added['append']:
                    return (False,True)
                case response if response is buttons_added['old_dup']:
                    return (False,False)
                case response if response is buttons_added['cancel']:
                    return None
                case _:
                    return None
            
    if button == QMessageBox.Yes:
        return button
