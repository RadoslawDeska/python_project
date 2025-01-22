from PyQt5.QtWidgets import (
    QWidget, QApplication, QProgressBar, QMainWindow,
    QHBoxLayout, QPushButton
)

from PyQt5.QtCore import (
    Qt, QObject, pyqtSignal, pyqtSlot, QRunnable, QThreadPool
)
import time


class WorkerSignals(QObject):
    progress = pyqtSignal(int)


class JobRunner(QRunnable):
    is_paused = False
    is_killed = False
        
    def __init__(self, func, *args, **kwargs):
        super().__init__()
        self.func = func
        self.args = args
        self.kwargs = kwargs

        self.signals = WorkerSignals()
        self.is_paused = False
        self.is_killed = False

    @pyqtSlot()
    def run(self):
        # Call the passed function (assumed to be a generator) and handle its output
        try:
            self.func(*self.args, **self.kwargs)                
        except Exception as e:
            print(f"Error in worker: {e}")
        print('Done')

    def pause(self):
        self.is_paused = True

    def resume(self):
        self.is_paused = False

    def kill(self):
        self.is_killed = True


class TaskHandler:
    def sample_task(self, is_paused):
        """ Sample task function to demonstrate usage. """
        while is_paused:
            time.sleep(0)
            
            if self.is_killed:
                break


class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()

        # Instantiate the task handler
        self.task_handler = TaskHandler()

        # Some buttons
        w = QWidget()
        l = QHBoxLayout()
        w.setLayout(l)

        btn_stop = QPushButton("Stop")
        btn_pause = QPushButton("Pause")
        btn_resume = QPushButton("Resume")

        l.addWidget(btn_stop)
        l.addWidget(btn_pause)
        l.addWidget(btn_resume)

        self.setCentralWidget(w)

        # Create a statusbar.
        self.status = self.statusBar()
        self.progress = QProgressBar()
        self.status.addPermanentWidget(self.progress)

        # Thread runner
        self.threadpool = QThreadPool()

        # Create a runner with a sample task
        self.runner = JobRunner(self.task_handler.sample_task, JobRunner.is_paused)
        self.runner.signals.progress.connect(self.update_progress)
        self.threadpool.start(self.runner)

        btn_stop.pressed.connect(self.runner.kill)
        btn_pause.pressed.connect(self.runner.pause)
        btn_resume.pressed.connect(self.runner.resume)

        self.show()

    def update_progress(self, n):
        self.progress.setValue(n)


if __name__ == "__main__":
    app = QApplication([])
    w = MainWindow()
    app.exec_()
