"""Chart Manager - Handles chart and visualization management.

Manages matplotlib figures, lines, and chart updates for measurement and fitting displays.
"""

from typing import Dict, List, Optional

from numpy.typing import NDArray
from PyQt5.QtWidgets import QWidget

from lib.figure import MplCanvas


class ChartManager:
    """Manages charts and visualizations for the Z-scan application."""

    def __init__(self, window: QWidget):
        """Initialize chart manager.

        Args:
            window: The main PyQt5 window
        """
        self.window = window
        self.charts: Dict[str, MplCanvas] = {}
        self.measurement_lines: Dict[str, List] = {}
        self.fitting_lines: Dict[str, List] = {}
        self.rms_text = None

    def setup_measurement_charts(self, chart_names: List[str]):
        """Setup measurement chart canvases.

        Args:
            chart_names: List of chart names (e.g., ["absolute", "relative"])
        """
        for chart_name in chart_names:
            # Create MplCanvas and add to appropriate layout
            canvas = MplCanvas(self.window, width=4, height=3, dpi=100)
            self.charts[chart_name] = canvas

            # Setup lines for each channel (typically 4)
            self.measurement_lines[chart_name] = []
            ax = canvas.axes
            for chan_no in range(4):
                (line,) = ax.plot([], [], label=f"CH{chan_no}")
                self.measurement_lines[chart_name].append(line)

            ax.legend()
            ax.set_xlabel("Position (mm)")
            ax.set_ylabel(
                "Intensity (V)"
                if chart_name == "absolute"
                else "Normalized Intensity"
            )

    def setup_fitting_charts(self, chart_names: List[str]):
        """Setup fitting chart canvases.

        Args:
            chart_names: List of fitting chart names
        """
        for chart_name in chart_names:
            canvas = MplCanvas(self.window, width=5, height=3, dpi=100)
            self.charts[chart_name] = canvas
            self.fitting_lines[chart_name] = []

    def update_measurement_lines(
        self,
        chart_name: str,
        positions: List[float],
        data: Dict[int, List[float]],
    ):
        """Update measurement line data.

        Args:
            chart_name: Name of the chart ("absolute" or "relative")
            positions: X-axis data (positions)
            data: Dictionary mapping channel numbers to Y-axis data
        """
        if chart_name not in self.measurement_lines:
            return

        for chan_no, line in enumerate(self.measurement_lines[chart_name]):
            if chan_no in data:
                line.set_xdata(positions)
                line.set_ydata(data[chan_no])

    def update_fitting_line(
        self,
        chart_name: str,
        x_data: NDArray,
        y_data: NDArray,
        label: Optional[str] = None,
    ):
        """Update or add a fitting line.

        Args:
            chart_name: Name of the fitting chart
            x_data: X-axis data
            y_data: Y-axis data
            label: Optional line label
        """
        if chart_name not in self.charts:
            return

        ax = self.charts[chart_name].axes
        (line,) = ax.plot(x_data, y_data, label=label)
        self.fitting_lines[chart_name].append(line)

    def clear_fitting_chart(self, chart_name: str):
        """Clear all lines from a fitting chart.

        Args:
            chart_name: Name of the fitting chart
        """
        if chart_name not in self.charts:
            return

        ax = self.charts[chart_name].axes
        ax.clear()
        self.fitting_lines[chart_name] = []

    def rescale_chart(self, chart_name: str, auto: bool = True):
        """Rescale chart axes.

        Args:
            chart_name: Name of the chart
            auto: Whether to use autoscale
        """
        if chart_name not in self.charts:
            return

        canvas = self.charts[chart_name]
        ax = canvas.axes

        if auto:
            ax.relim()
            ax.autoscale_view()

        canvas.draw()

    def draw_all_charts(self):
        """Redraw all active charts."""
        for canvas in self.charts.values():
            canvas.draw()

    def set_rms_text(self, text: str):
        """Set RMS noise display text.

        Args:
            text: Text to display
        """
        if self.rms_text:
            self.rms_text.set_text(text)
