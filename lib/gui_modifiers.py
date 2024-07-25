from PyQt5.QtWidgets import QLayout, QWidget

def cover_widget(widget):
    if isinstance(widget, QLayout):
        for child in widget.parent().findChildren(QWidget):
            try:
                csizepolicy = widget.sizePolicy()
                csizepolicy.setRetainSizeWhenHidden(True)
                child.setSizePolicy(csizepolicy)
                child.hide()
            except AttributeError:
                pass
    else:
        wsizepolicy = widget.sizePolicy()
        wsizepolicy.setRetainSizeWhenHidden(True)
        widget.setSizePolicy(wsizepolicy)
        widget.hide()