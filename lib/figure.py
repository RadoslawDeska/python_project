from cycler import cycler
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.ticker import FormatStrFormatter
matplotlib.rcParams['axes.prop_cycle'] = cycler(color=['orange','#3c7ffc','magenta'])

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_axes([0.07,0.185,.922,.78]) # [left, bottom, width, height]
        self.axes.set_ylabel('Amplitude')
        self.axes.set_xlabel('Position (mm)')
        self.axes.minorticks_on()
        mpl_version = [int(i) for i in matplotlib.__version__.split('.')]
        if mpl_version[0] >= 3 and mpl_version[1] < 5:
            self.axes.grid(b=True, which='both', axis='both')
        else:
            self.axes.grid(visible=True, which='both', axis='both')
        self.axes.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))
        super(MplCanvas, self).__init__(self.fig)