# INITIALIZE CHARTS
from matplotlib.offsetbox import AnchoredText
from lib.figure import MplCanvas
import matplotlib.pyplot as plt

# Set antialiasing for Matplotlib
plt.rcParams['lines.antialiased'] = True
plt.rcParams['patch.antialiased'] = True

def initialize_measurement_charts(self):
    '''Initializes charts in the 'Measurement' Tab
    and sets their scales using 'rescale_measurement_plots()' method.'''
    # "Measurement charts"
    self.rel_chart = MplCanvas(self)
    layout_rel = self.relative_layout
    layout_rel.addWidget(self.rel_chart)

    self.abs_chart = MplCanvas(self)
    self.abs_chart.axes.set_ylabel('Amplitude (V)')
    layout_abs = self.absolute_layout
    layout_abs.addWidget(self.abs_chart)
    
    self.rms_text = AnchoredText(f"RMS noise = {self.rms_value*100:.3f}%",
                prop=dict(size=8),frameon=True, loc='upper right')
    self.rms_text.patch.set_alpha(0.5)
    self.rms_text.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    self.abs_chart.axes.add_artist(self.rms_text)
    
    self.charts = {"relative": self.rel_chart, "absolute": self.abs_chart}
    # initialize empty lines and data dictionaries
    self.measurement_lines = {"relative": {}, "absolute": {}}

def initialize_fitting_charts(self):
    # Silica chart
    self.silica_figure = MplCanvas(self)
    self.silica_figure.axes.plot([],[], marker='o', ms=5, linestyle='', color='tab:blue')
    self.silica_figure.axes.set_title("Silica - closed aperture")
    self.silica_figure.axes.set_ylabel("Norm. trasmittance")
    self.silica_figure.axes.set_position([0.135,0.105,.825,.825])
    layout_ca_silica = self.silicaCA_layout
    layout_ca_silica.addWidget(self.silica_figure)
    
    self.silicaOA_figure = MplCanvas(self)
    self.silicaOA_figure.axes.plot([],[], marker='o', ms=5, linestyle='', color='tab:blue')
    self.silicaOA_figure.axes.set_title("Silica - open aperture")
    self.silicaOA_figure.axes.set_ylabel("Norm. trasmittance")
    self.silicaOA_figure.axes.set_position([0.135,0.105,.825,.825])
    layout_oa_silica = self.silicaOA_layout
    layout_oa_silica.addWidget(self.silicaOA_figure)

    # Solvent charts
    self.solventCA_figure = MplCanvas(self)
    self.solventCA_figure.axes.plot([],[], marker='o', ms=5, linestyle='', color='tab:blue')
    self.solventCA_figure.axes.set_title("Solvent - closed aperture")
    self.solventCA_figure.axes.set_ylabel("Norm. trasmittance")
    self.solventCA_figure.axes.set_position([0.135,0.105,.825,.825])
    layout_ca_solvent = self.solventCA_layout
    layout_ca_solvent.addWidget(self.solventCA_figure)
    
    self.solventOA_figure = MplCanvas(self)
    self.solventOA_figure.axes.plot([],[], marker='o', ms=5, linestyle='', color='tab:blue')
    self.solventOA_figure.axes.set_title("Solvent - open aperture")
    self.solventOA_figure.axes.set_ylabel("Norm. trasmittance")
    self.solventOA_figure.axes.set_position([0.135,0.105,.825,.825])
    layout_oa_solvent = self.solventOA_layout
    layout_oa_solvent.addWidget(self.solventOA_figure)

    # Sample charts
    self.sampleCA_figure = MplCanvas(self)
    self.sampleCA_figure.axes.plot([],[], marker='o', ms=5, linestyle='', color='tab:blue')
    self.sampleCA_figure.axes.set_title("Sample - closed aperture")
    self.sampleCA_figure.axes.set_ylabel("Norm. trasmittance")
    self.sampleCA_figure.axes.set_position([0.135,0.105,.825,.825]) # [left, bottom, width, height]
    layout_ca_sample = self.sampleCA_layout
    layout_ca_sample.addWidget(self.sampleCA_figure)
    
    self.sampleOA_figure = MplCanvas(self)
    self.sampleOA_figure.axes.plot([],[], marker='o', ms=5, linestyle='', color='tab:blue')
    self.sampleOA_figure.axes.set_title("Sample - open aperture")
    self.sampleOA_figure.axes.set_ylabel("Norm. trasmittance")
    self.sampleOA_figure.axes.set_position([0.135,0.105,.825,.825])
    layout_oa_sample = self.sampleOA_layout
    layout_oa_sample.addWidget(self.sampleOA_figure)

    self.fitting_charts = {"Silica": {"CA": self.silica_figure, "OA": self.silicaOA_figure},
                            "Solvent": {"CA": self.solventCA_figure, "OA": self.solventOA_figure},
                            "Sample": {"CA": self.sampleCA_figure, "OA": self.sampleOA_figure}}

def configure_charts(self):
    for signal_type, chart in self.charts.items():             # e.g.: take the tuple ("relative", "self.rel_chart")
        for chan_no, color in zip(range(self.number_of_channels_used), ["#ffcc6e","#a6e4fd","#ff9cff"]):
            self.data[signal_type].update({chan_no: []})       # then fill "relative" dictionary in "self.data" dictionary with pairs of channel number and list of values.
            line, = chart.axes.plot(self.data["positions"],self.data[signal_type][chan_no], marker='.', zorder=2) # add empty line for each channel on the "relative" chart
            mean_line, = chart.axes.plot(self.data["positions"],self.data[signal_type][chan_no], color=color, linewidth=2.0)
            self.measurement_lines[signal_type].update({chan_no: {"current": line, "mean": mean_line}})    # update the "lines" dictionary with line for each channel