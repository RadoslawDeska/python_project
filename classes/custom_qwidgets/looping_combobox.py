from PyQt5.QtWidgets import QComboBox
from PyQt5.QtCore import Qt

class LoopingComboBox(QComboBox):
    def __init__(self, *args):
        super().__init__(*args)

    def keyPressEvent(self, event):
        # If the user presses the up arrow
        if event.key() == Qt.Key_Up:
            if self.currentIndex() == 0:
                self.setCurrentIndex(self.count() - 1)
                return  # Avoid calling the base implementation
        # If the user presses the down arrow
        elif event.key() == Qt.Key_Down:
            if self.currentIndex() == self.count() - 1:
                self.setCurrentIndex(0)
                return  # Avoid calling the base implementation

        # Call the base implementation for other key events
        super().keyPressEvent(event)

    def wheelEvent(self, event):
        # Determine the scroll direction
        if event.angleDelta().y() > 0:  # Scrolling up
            if self.currentIndex() == 0:
                self.setCurrentIndex(self.count() - 1)  # Loop to the last item
            else:
                self.setCurrentIndex(self.currentIndex() - 1)  # Move up
        else:  # Scrolling down
            if self.currentIndex() == self.count() - 1:
                self.setCurrentIndex(0)  # Loop to the first item
            else:
                self.setCurrentIndex(self.currentIndex() + 1)  # Move down

        # No need to call the base implementation; we handle everything here
        event.accept()  # Accept the event so it doesn't propagate further