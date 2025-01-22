from PyQt5.QtCore import QObject, pyqtSignal, QRunnable, pyqtSlot

import traceback
import sys

'''
https://www.pythonguis.com/tutorials/multithreading-pyqt-applications-qthreadpool/
'''

class WorkerSignals(QObject):
    '''
    Defines the signals available from a running worker thread.

    Supported signals are:

    finished
        No data

    error
        tuple (exctype, value, traceback.format_exc() )

    result
        object data returned from processing, anything

    progress
        int indicating % progress
        
    data_ready
        (int, int, object, str, bool) represents data (object)
        when ready and additional info about signal occurence

    '''
    finished = pyqtSignal()
    error = pyqtSignal(tuple)
    result = pyqtSignal(object)
    progress = pyqtSignal(int)
    job_starting = pyqtSignal()
    data_ready = pyqtSignal(int, int, object)
    acquisition_complete = pyqtSignal(bool, int)
    motion_finished = pyqtSignal(bool)


class Worker(QRunnable):
    '''
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thread. Supplied args and
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function

    '''

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()
        
        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        # Add the callback to our kwargs
        self.kwargs['progress_callback'] = self.signals.progress
        self.kwargs['data_ready_callback'] = self.signals.data_ready
        self.kwargs['acquisition_complete'] = self.signals.acquisition_complete
        self.kwargs['motion_finished_callback'] = self.signals.motion_finished

    @pyqtSlot()
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''

        # Retrieve args/kwargs here; and fire processing using them
        try:
            self.signals.job_starting.emit()
            result = self.fn(*self.args, **self.kwargs)
        except Exception:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)  # Return the result of the processing
        finally:
            self.signals.finished.emit()  # Done
    
    def resume(self):
        if 'mpositioner' in self.kwargs:
            self.kwargs['mpositioner'].resumed = True