"""
Z-Scan Data Plotter with Real Workflow Fitting
Uses actual lib.fitting and lib.integration for physics-based fitting
"""

import sys
import numpy as np
from pathlib import Path
from typing import Optional

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QLineEdit, QPushButton, QComboBox, QFileDialog,
    QGroupBox, QGridLayout, QDoubleSpinBox, QStatusBar
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# Import your actual physics modules
from lib.fitting import Fitting
from lib.integration import Integration
from lib.config import N_COMPONENTS, INTEGRATION_STEPS, SILICA_BETA


class ParameterEstimator:
    """Estimate fitting parameters from data geometry and physics"""
    
    @staticmethod
    def infer_parameters_from_data(data: np.ndarray, wavelength_nm: float) -> dict:
        """
        Infer DPhi0 and beamwaist from CA data geometry.
        
        Uses the characteristic Z-scan signature:
        - Peak-to-valley separation relates to Rayleigh range z0
        - Peak-to-valley amplitude relates to DPhi0
        """
        wavelength = wavelength_nm * 1e-9
        
        # Get CA antisymmetric data
        ca_data = data['ca_antisym']
        positions = data['position_centered']  # in mm
        
        # Find peak (maximum) and valley (minimum)
        max_idx = np.argmax(ca_data)
        min_idx = np.argmin(ca_data)
        
        max_val = ca_data[max_idx]
        min_val = ca_data[min_idx]
        max_pos = positions[max_idx] * 1e-3  # Convert to m
        min_pos = positions[min_idx] * 1e-3  # Convert to m
        
        # Peak-to-valley separation in z
        z_pv = np.abs(max_pos - min_pos)
        
        # Peak-to-valley amplitude
        pv_amplitude = max_val - min_val
        
        # For a CA Z-scan: z0 ≈ z_pv / 1.7 (empirical relation)
        # This comes from the characteristic shape of dφ0 / (1 + z²/z0²)
        z0_estimate = z_pv / 1.7
        
        # Beamwaist: z0 = π * w0² / λ
        # w0 = sqrt(z0 * λ / π)
        beamwaist_estimate = np.sqrt(z0_estimate * wavelength / np.pi)
        
        # DPhi0 estimate from amplitude
        # The peak amplitude ≈ 0.87 * DPhi0 for the antisymmetric CA signature
        dphi0_estimate = pv_amplitude / 0.87
        
        return {
            'dphi0': dphi0_estimate,
            'beamwaist': beamwaist_estimate,
            'z0': z0_estimate,
            'pv_amplitude': pv_amplitude,
            'z_pv': z_pv,
        }


class ZScanDataProcessor:
    """Parse and process Z-scan data files"""
    
    @staticmethod
    def parse_file(filepath: str) -> Optional[dict]:
        """Parse Z-scan data file and extract header info + data"""
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            # Parse extended header
            start_pos = 40
            end_pos = 80
            wavelength = 800  # default in nm
            sample_type = 'Unknown'
            code = 'Unknown'
            header_end = 0
            
            for i, line in enumerate(lines):
                if 'Sample type:' in line:
                    sample_type = line.split(':')[1].strip()
                if 'Code:' in line:
                    code = line.split(':')[1].strip()
                if 'Wavelength:' in line:
                    wavelength_str = line.split(':')[1].strip()
                    wavelength = float(wavelength_str.split()[0])
                if 'Starting pos:' in line:
                    start_pos = int(line.split(':')[1].strip())
                if 'Ending pos:' in line:
                    end_pos = int(line.split(':')[1].strip())
                if 'SNo.' in line:
                    header_end = i + 2
                    break
            
            # Parse data rows
            data = []
            for line in lines[header_end:]:
                line = line.strip()
                if not line or line.startswith('-'):
                    continue
                
                try:
                    values = [float(v) for v in line.split() if v]
                    if len(values) >= 4:
                        data.append({
                            'sno': values[0],
                            'ca_raw': values[1],
                            'ref': values[2],
                            'oa_raw': values[3],
                        })
                except ValueError:
                    continue
            
            if not data:
                return None
            
            return {
                'start_pos': start_pos,
                'end_pos': end_pos,
                'wavelength': wavelength,
                'sample_type': sample_type,
                'code': code,
                'data': data,
            }
        
        except Exception as e:
            raise ValueError(f"Error parsing file: {str(e)}")
    
    @staticmethod
    def process_data(raw_data: list, start_pos: float, end_pos: float) -> np.ndarray:
        """
        Process raw data into normalized CA and OA signals
        Returns structured array with processed data
        """
        n = len(raw_data)
        processed = np.zeros(n, dtype=[
            ('index', int),
            ('position_mm', float),
            ('position_centered', float),
            ('ca_raw', float),
            ('ref', float),
            ('oa_raw', float),
            ('ca', float),
            ('oa', float),
            ('ca_antisym', float),
        ])
        
        # Calculate positions and normalize
        for idx, row in enumerate(raw_data):
            position_mm = start_pos + (idx / (n - 1)) * (end_pos - start_pos)
            range_mm = end_pos - start_pos
            position_centered = position_mm - (start_pos + range_mm / 2)
            
            ca = row['ca_raw'] / row['ref']
            oa = row['oa_raw'] / row['ref']
            
            processed['index'][idx] = idx
            processed['position_mm'][idx] = position_mm
            processed['position_centered'][idx] = position_centered
            processed['ca_raw'][idx] = row['ca_raw']
            processed['ref'][idx] = row['ref']
            processed['oa_raw'][idx] = row['oa_raw']
            processed['ca'][idx] = ca
            processed['oa'][idx] = oa
        
        # Normalize to zero level = 1.0 (shift so data oscillates around 1, not 0)
        ca_mean = np.mean(processed['ca'])
        oa_mean = np.mean(processed['oa'])
        
        processed['ca'] = 1.0 + (processed['ca'] - ca_mean)
        processed['oa'] = 1.0 + (processed['oa'] - oa_mean)
        
        # Antisymmetrize CA around center while preserving the 1.0 baseline
        center_idx = n // 2
        
        for idx in range(n):
            dist_from_center = idx - center_idx
            mirror_idx = center_idx - dist_from_center
            
            if 0 <= mirror_idx < n:
                mirror_value = processed['ca'][mirror_idx]
                # Antisymmetric deviation from 1.0
                processed['ca_antisym'][idx] = 1.0 + (processed['ca'][idx] - 1.0) - (mirror_value - 1.0)
            else:
                processed['ca_antisym'][idx] = processed['ca'][idx]
        
        return processed


class PlotCanvas(FigureCanvas):
    """Matplotlib canvas for plotting"""
    
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(10, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        super().__init__(self.fig)
        self.setParent(parent)
    
    def plot_data(
        self,
        data: np.ndarray,
        fit_curve: Optional[np.ndarray] = None,
        scan_type: str = 'CA',
    ):
        """Plot Z-scan data and optional fit"""
        self.ax.clear()
        
        if scan_type == 'CA':
            y_data = data['ca_antisym']
            y_label = 'ΔT/T₀ (Normalized)'
            title = 'Closed Aperture Z-Scan (Antisymmetric)'
        else:
            y_data = data['oa']
            y_label = '1 - T (Normalized)'
            title = 'Open Aperture Z-Scan'
        
        # Plot data points
        self.ax.scatter(
            data['position_centered'],
            y_data,
            color='#3b82f6',
            s=30,
            alpha=0.7,
            label='Data',
            zorder=3
        )
        
        # Plot fit if available
        if fit_curve is not None:
            self.ax.plot(
                data['position_centered'],
                fit_curve,
                color='#ef4444',
                linewidth=2.5,
                label='Fit (Physics)',
                zorder=5
            )
        
        self.ax.set_xlabel('Centered Position (mm)', fontsize=11)
        self.ax.set_ylabel(y_label, fontsize=11)
        self.ax.set_title(title, fontsize=12, fontweight='bold')
        self.ax.grid(True, alpha=0.3)
        self.ax.legend(fontsize=10)
        self.fig.tight_layout()
        self.draw()


class ZScanPlotterApp(QMainWindow):
    """Main application window"""
    
    def __init__(self):
        super().__init__()
        self.data = None
        self.file_wavelength = 800
        self.file_code = 'Unknown'
        self.z_range_mm = 40  # Default
        self.init_ui()
    
    def init_ui(self):
        """Initialize the UI"""
        self.setWindowTitle('Z-Scan Physics-Based Plotter')
        self.setGeometry(100, 100, 1400, 800)
        
        main_widget = QWidget()
        main_layout = QHBoxLayout()
        
        left_panel = self.create_control_panel()
        main_layout.addWidget(left_panel, 0)
        
        self.canvas = PlotCanvas(self)
        main_layout.addWidget(self.canvas, 1)
        
        main_widget.setLayout(main_layout)
        self.setCentralWidget(main_widget)
        
        self.status = QStatusBar()
        self.setStatusBar(self.status)
        self.status.showMessage('Ready')
    
    def create_control_panel(self) -> QWidget:
        """Create the left control panel"""
        panel = QWidget()
        layout = QVBoxLayout()
        
        # File loading
        file_group = QGroupBox('Data File')
        file_layout = QVBoxLayout()
        self.file_button = QPushButton('Load Z-Scan Data')
        self.file_button.clicked.connect(self.load_file)
        self.file_label = QLabel('No file loaded')
        self.file_label.setWordWrap(True)
        file_layout.addWidget(self.file_button)
        file_layout.addWidget(self.file_label)
        file_group.setLayout(file_layout)
        layout.addWidget(file_group)
        
        # Position range
        pos_group = QGroupBox('Position Range')
        pos_layout = QGridLayout()
        pos_layout.addWidget(QLabel('Start (mm):'), 0, 0)
        self.start_display = QLineEdit()
        self.start_display.setReadOnly(True)
        pos_layout.addWidget(self.start_display, 0, 1)
        pos_layout.addWidget(QLabel('End (mm):'), 1, 0)
        self.end_display = QLineEdit()
        self.end_display.setReadOnly(True)
        pos_layout.addWidget(self.end_display, 1, 1)
        pos_group.setLayout(pos_layout)
        layout.addWidget(pos_group)
        
        # Scan configuration
        scan_group = QGroupBox('Scan Configuration')
        scan_layout = QGridLayout()
        scan_layout.addWidget(QLabel('Scan Type:'), 0, 0)
        self.scan_type_combo = QComboBox()
        self.scan_type_combo.addItems(['CA', 'OA'])
        scan_layout.addWidget(self.scan_type_combo, 0, 1)
        scan_group.setLayout(scan_layout)
        layout.addWidget(scan_group)
        
        # Physical parameters
        param_group = QGroupBox('Physical Parameters')
        param_layout = QGridLayout()
        
        # DPhi0 / T value
        param_layout.addWidget(QLabel('DPhi0 / T:'), 0, 0)
        self.param_amplitude = QDoubleSpinBox()
        self.param_amplitude.setRange(-10, 10)
        self.param_amplitude.setSingleStep(0.01)
        self.param_amplitude.setValue(0.5)
        param_layout.addWidget(self.param_amplitude, 0, 1)
        
        # Beamwaist
        param_layout.addWidget(QLabel('Beamwaist (µm):'), 1, 0)
        self.param_beamwaist = QDoubleSpinBox()
        self.param_beamwaist.setRange(0.1, 100)
        self.param_beamwaist.setSingleStep(0.1)
        self.param_beamwaist.setValue(20)
        param_layout.addWidget(self.param_beamwaist, 1, 1)
        
        # Zero level
        param_layout.addWidget(QLabel('Zero Level:'), 2, 0)
        self.param_zero_level = QDoubleSpinBox()
        self.param_zero_level.setRange(0.5, 1.5)
        self.param_zero_level.setSingleStep(0.01)
        self.param_zero_level.setValue(1.0)
        param_layout.addWidget(self.param_zero_level, 2, 1)
        
        # n2 (for CA)
        param_layout.addWidget(QLabel('n2 (auto from λ):'), 3, 0)
        self.param_n2_display = QLineEdit()
        self.param_n2_display.setReadOnly(True)
        self.param_n2_display.setText('—')
        param_layout.addWidget(self.param_n2_display, 3, 1)
        
        # d0 (aperture distance)
        param_layout.addWidget(QLabel('d0 (mm):'), 4, 0)
        self.param_d0 = QDoubleSpinBox()
        self.param_d0.setRange(0.1, 500)
        self.param_d0.setSingleStep(1)
        self.param_d0.setValue(100)
        param_layout.addWidget(self.param_d0, 4, 1)
        
        # Aperture radius
        param_layout.addWidget(QLabel('Aperture (mm):'), 5, 0)
        self.param_ra = QDoubleSpinBox()
        self.param_ra.setRange(0.1, 10)
        self.param_ra.setSingleStep(0.1)
        self.param_ra.setValue(1)
        param_layout.addWidget(self.param_ra, 5, 1)
        
        param_group.setLayout(param_layout)
        layout.addWidget(param_group)
        
        # Generate fit button
        self.fit_button = QPushButton('Generate Physics Fit')
        self.fit_button.setFont(QFont('Arial', 11, QFont.Bold))
        self.fit_button.setEnabled(False)
        self.fit_button.clicked.connect(self.generate_fit)
        layout.addWidget(self.fit_button)
        
        # Summary
        summary_group = QGroupBox('Summary')
        summary_layout = QVBoxLayout()
        self.summary_label = QLabel('No data loaded')
        self.summary_label.setWordWrap(True)
        summary_layout.addWidget(self.summary_label)
        summary_group.setLayout(summary_layout)
        layout.addWidget(summary_group)
        
        layout.addStretch()
        panel.setLayout(layout)
        return panel
    
    def load_file(self):
        """Load a Z-scan data file"""
        filepath, _ = QFileDialog.getOpenFileName(
            self, 'Open Z-Scan Data', '', 'Text Files (*.txt);;CSV Files (*.csv)'
        )
        
        if not filepath:
            return
        
        try:
            result = ZScanDataProcessor.parse_file(filepath)
            if not result:
                self.status.showMessage('Error: No valid data found')
                return
            
            self.data = ZScanDataProcessor.process_data(
                result['data'],
                result['start_pos'],
                result['end_pos']
            )
            
            self.file_wavelength = result['wavelength']
            self.file_code = result['code']
            self.z_range_mm = result['end_pos'] - result['start_pos']
            
            # Calculate and display n2 for silica at this wavelength
            wavelength_m = self.file_wavelength * 1e-9
            silica_n2 = 2.8203e-20 - 3e-27 / wavelength_m + 2e-33 / (wavelength_m ** 2)
            self.param_n2_display.setText(f'{silica_n2:.3e} m²/W')
            
            self.start_display.setText(str(result['start_pos']))
            self.end_display.setText(str(result['end_pos']))
            
            # Update summary
            summary_text = (
                f'Data Points: {len(self.data)}\n'
                f'Range: {self.z_range_mm:.1f} mm\n'
                f'Wavelength: {self.file_wavelength:.0f} nm\n'
                f'Code: {self.file_code}\n'
                f'Scan Type: {self.scan_type_combo.currentText()}'
            )
            self.summary_label.setText(summary_text)
            
            # Infer initial parameters from data geometry
            if self.scan_type_combo.currentText() == 'CA':
                params = ParameterEstimator.infer_parameters_from_data(
                    self.data, self.file_wavelength
                )
                self.param_amplitude.setValue(params['dphi0'])
                self.param_beamwaist.setValue(params['beamwaist'] * 1e6)  # to µm
            
            self.fit_button.setEnabled(True)
            self.canvas.plot_data(self.data, scan_type=self.scan_type_combo.currentText())
            
            self.file_label.setText(f'Loaded: {Path(filepath).name}')
            self.status.showMessage(f'Loaded {len(self.data)} data points')
        
        except Exception as e:
            self.status.showMessage(f'Error: {str(e)}')
    
    def generate_fit(self):
        """Generate physics-based fit using actual Integration and Fitting classes"""
        if self.data is None:
            self.status.showMessage('Error: No data loaded')
            return
        
        try:
            # Get parameters
            amplitude = self.param_amplitude.value()
            beamwaist_um = self.param_beamwaist.value()
            beamwaist = beamwaist_um * 1e-6  # Convert to meters
            zero_level = self.param_zero_level.value()
            d0 = self.param_d0.value() * 1e-3  # mm to m
            ra = self.param_ra.value() * 1e-3  # mm to m
            wavelength = self.file_wavelength * 1e-9  # nm to m
            
            # Calculate silica n2 from wavelength-dependent formula
            silica_n2 = (
                2.8203e-20 - 3e-27 / wavelength + 2e-33 / (wavelength ** 2)
            )  # [m²/W]
            
            z_range = self.z_range_mm * 1e-3  # mm to m
            scan_type = self.scan_type_combo.currentText()
            
            # Create z-positions array (same as data)
            n_points = len(self.data)
            z_positions = np.array(self.data['position_centered']) * 1e-3  # mm to m
            
            # Create Integration object - use silica n2 for CA, 0 for OA
            beta = amplitude if scan_type == 'OA' else 0
            dphi0 = amplitude if scan_type == 'CA' else 0
            n2 = silica_n2 if scan_type == 'CA' else 0
            
            integration = Integration(
                beta=beta,
                n2=n2,
                DPhi0=dphi0,
                positions=z_positions,
                d0=d0,
                aperture_radius=ra,
                wavelength=wavelength,
                beamwaist=beamwaist,
                n_components=N_COMPONENTS,
                integration_steps=INTEGRATION_STEPS,
                stype=scan_type,
            )
            
            # Create Fitting object
            y_data = self.data[f'{scan_type.lower()}_antisym' if scan_type == 'CA' else f'{scan_type.lower()}']
            
            fitter = Fitting(
                sample_type=integration,
                amplitude=amplitude,
                beamwaist=beamwaist,
                zero_level=zero_level,
                centerpoint=0,
                nop=n_points,
                data=y_data,
            )
            
            # Generate fit curve
            fit_curve = fitter.manual(
                zero_level=zero_level,
                centerpoint=0,
                amplitude=amplitude,
                beamwaist=beamwaist,
                z_range=z_range,
                window=None,  # Headless mode
                stype=scan_type,
            )
            
            # Plot with fit
            self.canvas.plot_data(self.data, fit_curve, scan_type=scan_type)
            
            # Show calculated n2 in status
            if scan_type == 'CA':
                self.status.showMessage(
                    f'Physics-based fit generated | n2={silica_n2:.3e} m²/W | DPhi0={amplitude:.4f} rad'
                )
            else:
                self.status.showMessage('Physics-based fit generated')
        
        except Exception as e:
            import traceback
            self.status.showMessage(f'Error: {str(e)}')
            traceback.print_exc()


def main():
    app = QApplication(sys.argv)
    window = ZScanPlotterApp()
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()