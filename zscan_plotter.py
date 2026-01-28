"""
Complete Z-Scan Application: window.ui + Physics-Based Fitting + Auto-Optimization
Full measurement analysis system with least-squares optimization
"""

import sys
import numpy as np
from pathlib import Path
from typing import Optional, Dict, Tuple, Callable
from dataclasses import dataclass
from threading import Thread
from enum import Enum

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QFileDialog, QVBoxLayout, 
    QWidget, QMessageBox, QProgressBar, QLabel
)
from PyQt5.QtCore import Qt, QTimer, pyqtSignal, QObject, QThread
from PyQt5.uic import loadUi

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

try:
    from lib.integration_headless import Integration, Fitting
    from lib.config import N_COMPONENTS, INTEGRATION_STEPS
except ImportError:
    N_COMPONENTS = 50
    INTEGRATION_STEPS = 100


class SampleType(Enum):
    SILICA = "silica"
    SOLVENT = "solvent"
    SAMPLE = "sample"


@dataclass
class ZScanFileData:
    """Container for parsed Z-scan data"""
    position_mm: np.ndarray
    ca_raw: np.ndarray
    ref: np.ndarray
    oa_raw: np.ndarray
    start_pos: float
    end_pos: float
    wavelength: float
    sample_type: str
    code: str


@dataclass
class FittingParams:
    """Fitting parameters container"""
    amplitude: float  # DPhi0 for CA, T for OA
    beamwaist: float  # meters
    zero_level: float
    centerpoint: float  # data points
    d0: float  # meters
    ra: float  # meters
    
    def to_dict(self) -> Dict[str, float]:
        return {
            'amplitude': self.amplitude,
            'beamwaist': self.beamwaist,
            'zero_level': self.zero_level,
            'centerpoint': self.centerpoint,
            'd0': self.d0,
            'ra': self.ra,
        }


@dataclass
class FittingResult:
    """Fitting result container"""
    params: FittingParams
    fit_curve: np.ndarray
    residuals: np.ndarray
    r_squared: float
    chi_squared: float
    n2: Optional[float] = None
    errors: Optional[Dict[str, float]] = None


class FitWorker(QObject):
    """Worker for running fits in background thread"""
    
    progress = pyqtSignal(str)
    finished = pyqtSignal(FittingResult)
    error = pyqtSignal(str)
    
    def __init__(self, fit_func: Callable, **kwargs):
        super().__init__()
        self.fit_func = fit_func
        self.kwargs = kwargs
    
    def run(self):
        try:
            result = self.fit_func(**self.kwargs)
            self.finished.emit(result)
        except Exception as e:
            self.error.emit(str(e))


class ZScanDataParser:
    """Parse Z-scan data files"""
    
    @staticmethod
    def parse(filepath: str) -> Optional[ZScanFileData]:
        """Parse Z-scan data file"""
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            metadata = {
                'start_pos': 40.0,
                'end_pos': 80.0,
                'wavelength': 800.0,
                'sample_type': 'Unknown',
                'code': 'Unknown',
            }
            
            header_end = 0
            for i, line in enumerate(lines):
                if 'Sample type:' in line:
                    metadata['sample_type'] = line.split(':')[1].strip()
                if 'Code:' in line:
                    metadata['code'] = line.split(':')[1].strip()
                if 'Wavelength:' in line:
                    wl_str = line.split(':')[1].strip().split()[0]
                    metadata['wavelength'] = float(wl_str)
                if 'Starting pos:' in line:
                    metadata['start_pos'] = float(line.split(':')[1].strip())
                if 'Ending pos:' in line:
                    metadata['end_pos'] = float(line.split(':')[1].strip())
                if 'SNo.' in line:
                    header_end = i + 2
                    break
            
            data_rows = []
            for line in lines[header_end:]:
                line = line.strip()
                if not line or line.startswith('-'):
                    continue
                try:
                    values = [float(v) for v in line.split() if v]
                    if len(values) >= 4:
                        data_rows.append(values)
                except ValueError:
                    continue
            
            if not data_rows:
                return None
            
            data_array = np.array(data_rows)
            n_points = len(data_array)
            
            pos_mm = np.linspace(
                metadata['start_pos'],
                metadata['end_pos'],
                n_points
            )
            
            return ZScanFileData(
                position_mm=pos_mm,
                ca_raw=data_array[:, 1],
                ref=data_array[:, 2],
                oa_raw=data_array[:, 3],
                start_pos=metadata['start_pos'],
                end_pos=metadata['end_pos'],
                wavelength=metadata['wavelength'],
                sample_type=metadata['sample_type'],
                code=metadata['code'],
            )
        
        except Exception as e:
            print(f"Parse error: {e}")
            return None


class ZScanProcessor:
    """Process raw Z-scan data"""
    
    @staticmethod
    def normalize(file_data: ZScanFileData) -> Dict[str, np.ndarray]:
        """Normalize raw data"""
        n = len(file_data.ca_raw)
        center_idx = n // 2
        
        ca = file_data.ca_raw / file_data.ref
        oa = file_data.oa_raw / file_data.ref
        
        ca = 1.0 + (ca - np.mean(ca))
        oa = 1.0 + (oa - np.mean(oa))
        
        ca_antisym = np.zeros(n)
        for i in range(n):
            mirror_i = 2 * center_idx - i
            if 0 <= mirror_i < n:
                ca_antisym[i] = 1.0 + (ca[i] - 1.0) - (ca[mirror_i] - 1.0)
            else:
                ca_antisym[i] = ca[i]
        
        z_range = file_data.end_pos - file_data.start_pos
        center_pos = file_data.start_pos + z_range / 2.0
        position_centered = file_data.position_mm - center_pos
        
        return {
            'position_mm': file_data.position_mm,
            'position_centered': position_centered,
            'ca_raw': file_data.ca_raw,
            'ca': ca,
            'ca_antisym': ca_antisym,
            'oa': oa,
            'ref': file_data.ref,
        }


class FitPlotCanvas(FigureCanvas):
    """Matplotlib canvas for Z-scan plots"""
    
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.ax = self.fig.add_subplot(111)
        super().__init__(self.fig)
        self.setParent(parent)
        self.fig.tight_layout()
    
    def plot(self, x: np.ndarray, y: np.ndarray, 
             y_fit: Optional[np.ndarray] = None, 
             title: str = '', xlabel: str = '', ylabel: str = '',
             residuals: Optional[np.ndarray] = None):
        """Plot data with optional fit and residuals"""
        if residuals is not None:
            # Create subplot with residuals
            self.fig.clear()
            ax1 = self.fig.add_subplot(2, 1, 1)
            ax2 = self.fig.add_subplot(2, 1, 2)
            
            ax1.scatter(x, y, color='#3b82f6', s=25, alpha=0.8, 
                       label='Data', zorder=3)
            if y_fit is not None:
                ax1.plot(x, y_fit, color='#ef4444', linewidth=2.5, 
                        label='Physics Fit', zorder=5)
            
            ax1.set_ylabel(ylabel, fontsize=10)
            ax1.set_title(title, fontsize=11, fontweight='bold')
            ax1.grid(True, alpha=0.25, linestyle='--')
            ax1.legend(fontsize=9)
            
            ax2.scatter(x, residuals, color='#8b5cf6', s=20, alpha=0.7)
            ax2.axhline(y=0, color='red', linestyle='--', linewidth=1)
            ax2.set_xlabel(xlabel, fontsize=10)
            ax2.set_ylabel('Residuals', fontsize=10)
            ax2.grid(True, alpha=0.25, linestyle='--')
            
            self.ax = ax1
        else:
            self.ax.clear()
            self.ax.scatter(x, y, color='#3b82f6', s=25, alpha=0.8, 
                           label='Data', zorder=3)
            if y_fit is not None:
                self.ax.plot(x, y_fit, color='#ef4444', linewidth=2.5, 
                            label='Physics Fit', zorder=5)
            
            self.ax.set_xlabel(xlabel, fontsize=10)
            self.ax.set_ylabel(ylabel, fontsize=10)
            self.ax.set_title(title, fontsize=11, fontweight='bold')
            self.ax.grid(True, alpha=0.25, linestyle='--')
            self.ax.legend(fontsize=9)
        
        self.fig.tight_layout()
        self.draw()


class ZScanPhysicsEngine:
    """Physics-based fitting with automatic optimization"""
    
    @staticmethod
    def estimate_silica_n2(wavelength_nm: float) -> float:
        """Estimate n2 for silica at given wavelength"""
        wavelength = wavelength_nm * 1e-9
        return 2.8203e-20 - 3e-27 / wavelength + 2e-33 / (wavelength ** 2)
    
    @staticmethod
    def estimate_solvent_n2(solvent_name: str) -> float:
        """Get n2 for common solvents"""
        solvent_n2 = {
            'water': 1.3e-20,
            'toluene': 2.8e-20,
            'chloroform': 1.8e-20,
            'acetone': 1.2e-20,
            'dmso': 2.5e-20,
        }
        return solvent_n2.get(solvent_name.lower(), 1.0e-20)
    
    @staticmethod
    def infer_initial_params(data: Dict[str, np.ndarray], 
                            wavelength_nm: float) -> FittingParams:
        """Infer DPhi0 and beamwaist from data geometry"""
        ca_antisym = data['ca_antisym']
        pos = data['position_centered'] * 1e-3
        
        max_idx = np.argmax(ca_antisym)
        min_idx = np.argmin(ca_antisym)
        
        z_pv = abs(pos[max_idx] - pos[min_idx])
        pv_amp = ca_antisym[max_idx] - ca_antisym[min_idx]
        
        wavelength = wavelength_nm * 1e-9
        z0 = z_pv / 1.7
        beamwaist = np.sqrt(z0 * wavelength / np.pi)
        dphi0 = pv_amp / 0.87
        
        return FittingParams(
            amplitude=dphi0,
            beamwaist=beamwaist,
            zero_level=1.0,
            centerpoint=0.0,
            d0=0.26,
            ra=0.001,
        )
    
    @staticmethod
    def calculate_residuals(y_true: np.ndarray, y_pred: np.ndarray) -> np.ndarray:
        """Calculate residuals"""
        return y_pred - y_true
    
    @staticmethod
    def calculate_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Tuple[float, float]:
        """Calculate R² and χ² metrics"""
        residuals = y_pred - y_true
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        chi_squared = ss_res / len(y_true)
        return r_squared, chi_squared
    
    @staticmethod
    def fit_ca_manual(data: Dict[str, np.ndarray],
                      wavelength_nm: float,
                      params: FittingParams,
                      n2: Optional[float] = None) -> Optional[np.ndarray]:
        """Manual CA fit with given parameters"""
        try:
            wavelength = wavelength_nm * 1e-9
            z_range = (data['position_mm'][-1] - data['position_mm'][0]) * 1e-3
            positions = data['position_centered'] * 1e-3
            
            if n2 is None:
                n2 = ZScanPhysicsEngine.estimate_silica_n2(wavelength_nm)
            
            integration = Integration(
                beta=0,
                n2=n2,
                DPhi0=params.amplitude,
                positions=positions,
                d0=params.d0,
                aperture_radius=params.ra,
                wavelength=wavelength,
                beamwaist=params.beamwaist,
                n_components=N_COMPONENTS,
                integration_steps=INTEGRATION_STEPS,
                stype='CA',
            )
            
            fitter = Fitting(
                integration=integration,
                amplitude=params.amplitude,
                beamwaist=params.beamwaist,
                zero_level=params.zero_level,
                centerpoint=params.centerpoint,
                nop=len(data['ca_antisym']),
                y_data=data['ca_antisym'],
            )
            
            fit_curve = fitter.manual(
                zero_level=params.zero_level,
                centerpoint=params.centerpoint,
                amplitude=params.amplitude,
                beamwaist=params.beamwaist,
                z_range=z_range,
                d0=params.d0,
                ra=params.ra,
                stype='CA',
            )
            
            return np.asarray(fit_curve, dtype=float)
        
        except Exception as e:
            print(f"Fit error: {e}")
            return None
    
    @staticmethod
    def fit_oa_manual(data: Dict[str, np.ndarray],
                      wavelength_nm: float,
                      params: FittingParams,
                      n2: Optional[float] = None) -> Optional[np.ndarray]:
        """Manual OA fit with given parameters"""
        try:
            wavelength = wavelength_nm * 1e-9
            z_range = (data['position_mm'][-1] - data['position_mm'][0]) * 1e-3
            positions = data['position_centered'] * 1e-3
            
            if n2 is None:
                n2 = ZScanPhysicsEngine.estimate_silica_n2(wavelength_nm)
            
            integration = Integration(
                beta=params.amplitude,
                n2=0,
                DPhi0=0,
                positions=positions,
                d0=params.d0,
                aperture_radius=params.ra,
                wavelength=wavelength,
                beamwaist=params.beamwaist,
                n_components=N_COMPONENTS,
                integration_steps=INTEGRATION_STEPS,
                stype='OA',
            )
            
            fitter = Fitting(
                integration=integration,
                amplitude=params.amplitude,
                beamwaist=params.beamwaist,
                zero_level=params.zero_level,
                centerpoint=params.centerpoint,
                nop=len(data['oa']),
                y_data=data['oa'],
            )
            
            fit_curve = fitter.manual(
                zero_level=params.zero_level,
                centerpoint=params.centerpoint,
                amplitude=params.amplitude,
                beamwaist=params.beamwaist,
                z_range=z_range,
                d0=params.d0,
                ra=params.ra,
                stype='OA',
            )
            
            return np.asarray(fit_curve, dtype=float)
        
        except Exception as e:
            print(f"OA Fit error: {e}")
            return None
    
    @staticmethod
    def auto_fit_ca(data: Dict[str, np.ndarray],
                    wavelength_nm: float,
                    initial_params: FittingParams,
                    n2: Optional[float] = None) -> FittingResult:
        """Automatic CA fitting with scipy optimization"""
        from scipy.optimize import minimize
        
        wavelength = wavelength_nm * 1e-9
        z_range = (data['position_mm'][-1] - data['position_mm'][0]) * 1e-3
        positions = data['position_centered'] * 1e-3
        y_true = data['ca_antisym']
        
        if n2 is None:
            n2 = ZScanPhysicsEngine.estimate_silica_n2(wavelength_nm)
        
        def objective(param_vec):
            """Objective function to minimize"""
            amplitude, bw, zero_lvl, center = param_vec
            
            if bw <= 0 or amplitude < -5 or amplitude > 5:
                return 1e6
            
            try:
                integration = Integration(
                    beta=0, n2=n2, DPhi0=amplitude,
                    positions=positions, d0=initial_params.d0,
                    aperture_radius=initial_params.ra,
                    wavelength=wavelength, beamwaist=bw,
                    n_components=N_COMPONENTS,
                    integration_steps=INTEGRATION_STEPS, stype='CA',
                )
                
                fitter = Fitting(
                    integration=integration, amplitude=amplitude,
                    beamwaist=bw, zero_level=zero_lvl,
                    centerpoint=center, nop=len(y_true), y_data=y_true,
                )
                
                y_fit = fitter.manual(
                    zero_level=zero_lvl, centerpoint=center,
                    amplitude=amplitude, beamwaist=bw,
                    z_range=z_range, d0=initial_params.d0,
                    ra=initial_params.ra, stype='CA',
                )
                
                residuals = y_fit - y_true
                return np.sum(residuals ** 2)
            except:
                return 1e6
        
        # Optimize
        x0 = [initial_params.amplitude, initial_params.beamwaist,
              initial_params.zero_level, initial_params.centerpoint]
        
        bounds = [(-5, 5), (1e-6, 100e-6), (0.8, 1.2), (-50, 50)]
        
        result = minimize(objective, x0, method='L-BFGS-B',
                         bounds=bounds, options={'maxiter': 100})
        
        amplitude, bw, zero_lvl, center = result.x
        
        fitted_params = FittingParams(
            amplitude=amplitude, beamwaist=bw,
            zero_level=zero_lvl, centerpoint=center,
            d0=initial_params.d0, ra=initial_params.ra,
        )
        
        y_fit = ZScanPhysicsEngine.fit_ca_manual(
            data, wavelength_nm, fitted_params, n2
        )
        
        residuals = ZScanPhysicsEngine.calculate_residuals(y_true, y_fit)
        r2, chi2 = ZScanPhysicsEngine.calculate_metrics(y_true, y_fit)
        
        return FittingResult(
            params=fitted_params,
            fit_curve=y_fit,
            residuals=residuals,
            r_squared=r2,
            chi_squared=chi2,
            n2=n2,
        )
    
    @staticmethod
    def auto_fit_oa(data: Dict[str, np.ndarray],
                    wavelength_nm: float,
                    initial_params: FittingParams,
                    n2: Optional[float] = None) -> FittingResult:
        """Automatic OA fitting with scipy optimization"""
        from scipy.optimize import minimize
        
        wavelength = wavelength_nm * 1e-9
        z_range = (data['position_mm'][-1] - data['position_mm'][0]) * 1e-3
        positions = data['position_centered'] * 1e-3
        y_true = data['oa']
        
        if n2 is None:
            n2 = ZScanPhysicsEngine.estimate_silica_n2(wavelength_nm)
        
        def objective(param_vec):
            """Objective function"""
            amplitude, bw, zero_lvl, center = param_vec
            
            if bw <= 0 or amplitude < -5 or amplitude > 5:
                return 1e6
            
            try:
                integration = Integration(
                    beta=amplitude, n2=0, DPhi0=0,
                    positions=positions, d0=initial_params.d0,
                    aperture_radius=initial_params.ra,
                    wavelength=wavelength, beamwaist=bw,
                    n_components=N_COMPONENTS,
                    integration_steps=INTEGRATION_STEPS, stype='OA',
                )
                
                fitter = Fitting(
                    integration=integration, amplitude=amplitude,
                    beamwaist=bw, zero_level=zero_lvl,
                    centerpoint=center, nop=len(y_true), y_data=y_true,
                )
                
                y_fit = fitter.manual(
                    zero_level=zero_lvl, centerpoint=center,
                    amplitude=amplitude, beamwaist=bw,
                    z_range=z_range, d0=initial_params.d0,
                    ra=initial_params.ra, stype='OA',
                )
                
                residuals = y_fit - y_true
                return np.sum(residuals ** 2)
            except:
                return 1e6
        
        x0 = [initial_params.amplitude, initial_params.beamwaist,
              initial_params.zero_level, initial_params.centerpoint]
        
        bounds = [(-5, 5), (1e-6, 100e-6), (0.8, 1.2), (-50, 50)]
        
        result = minimize(objective, x0, method='L-BFGS-B',
                         bounds=bounds, options={'maxiter': 100})
        
        amplitude, bw, zero_lvl, center = result.x
        
        fitted_params = FittingParams(
            amplitude=amplitude, beamwaist=bw,
            zero_level=zero_lvl, centerpoint=center,
            d0=initial_params.d0, ra=initial_params.ra,
        )
        
        y_fit = ZScanPhysicsEngine.fit_oa_manual(
            data, wavelength_nm, fitted_params, n2
        )
        
        residuals = ZScanPhysicsEngine.calculate_residuals(y_true, y_fit)
        r2, chi2 = ZScanPhysicsEngine.calculate_metrics(y_true, y_fit)
        
        return FittingResult(
            params=fitted_params,
            fit_curve=y_fit,
            residuals=residuals,
            r_squared=r2,
            chi_squared=chi2,
            n2=n2,
        )


class MainWindow(QMainWindow):
    """Main application window with full fitting capabilities"""
    
    def __init__(self, ui_path: str = "window.ui"):
        super().__init__()
        
        try:
            self.ui = loadUi(ui_path, self)
        except Exception as e:
            print(f"Warning: Could not load UI: {e}")
            self.create_minimal_ui()
        
        self.setWindowTitle('Z-Scan Physics Analysis System - Full Auto-Fit Edition')
        
        # Data storage
        self.file_data: Dict[SampleType, Optional[ZScanFileData]] = {
            t: None for t in SampleType
        }
        self.processed_data: Dict[SampleType, Optional[Dict]] = {
            t: None for t in SampleType
        }
        self.fitting_results: Dict[Tuple[SampleType, str], Optional[FittingResult]] = {}
        
        # Threading
        self.fit_thread: Optional[QThread] = None
        self.fit_worker: Optional[FitWorker] = None
        
        self._setup_plot_canvases()
        self._connect_ui_elements()
        
        if hasattr(self, 'statusbar'):
            self.statusbar.showMessage('Ready - Load data files to begin')
    
    def create_minimal_ui(self):
        """Create minimal UI if file loading fails"""
        from PyQt5.QtWidgets import (
            QWidget, QVBoxLayout, QPushButton, QLabel, QStatusBar
        )
        
        central = QWidget()
        layout = QVBoxLayout()
        layout.addWidget(QLabel('Z-Scan Physics Fitting System'))
        
        self.load_silica_btn = QPushButton('Load Silica Data')
        self.load_silica_btn.clicked.connect(lambda: self.load_data_file(SampleType.SILICA))
        layout.addWidget(self.load_silica_btn)
        
        self.load_solvent_btn = QPushButton('Load Solvent Data')
        self.load_solvent_btn.clicked.connect(lambda: self.load_data_file(SampleType.SOLVENT))
        layout.addWidget(self.load_solvent_btn)
        
        self.load_sample_btn = QPushButton('Load Sample Data')
        self.load_sample_btn.clicked.connect(lambda: self.load_data_file(SampleType.SAMPLE))
        layout.addWidget(self.load_sample_btn)
        
        self.fit_silica_ca_btn = QPushButton('Auto-Fit Silica CA')
        self.fit_silica_ca_btn.setEnabled(False)
        self.fit_silica_ca_btn.clicked.connect(lambda: self.auto_fit(SampleType.SILICA, 'CA'))
        layout.addWidget(self.fit_silica_ca_btn)
        
        self.fit_silica_oa_btn = QPushButton('Auto-Fit Silica OA')
        self.fit_silica_oa_btn.setEnabled(False)
        self.fit_silica_oa_btn.clicked.connect(lambda: self.auto_fit(SampleType.SILICA, 'OA'))
        layout.addWidget(self.fit_silica_oa_btn)
        
        layout.addStretch()
        central.setLayout(layout)
        self.setCentralWidget(central)
        
        self.statusbar = QStatusBar()
        self.setStatusBar(self.statusbar)
    
    def _setup_plot_canvases(self):
        """Add matplotlib canvases to tabs"""
        layouts = {
            'silica': ('silicaCA_layout', 'silicaOA_layout'),
            'solvent': ('solventCA_layout', 'solventOA_layout'),
            'sample': ('sampleCA_layout', 'sampleOA_layout'),
        }
        
        self.canvases = {}
        for sample_type, (ca_layout_name, oa_layout_name) in layouts.items():
            try:
                if hasattr(self, ca_layout_name):
                    ca_canvas = FitPlotCanvas(self, width=4, height=3)
                    getattr(self, ca_layout_name).addWidget(ca_canvas)
                    self.canvases[f'{sample_type}_ca'] = ca_canvas
                
                if hasattr(self, oa_layout_name):
                    oa_canvas = FitPlotCanvas(self, width=4, height=3)
                    getattr(self, oa_layout_name).addWidget(oa_canvas)
                    self.canvases[f'{sample_type}_oa'] = oa_canvas
            except Exception as e:
                print(f"Could not set up {sample_type} canvases: {e}")
    
    def _connect_ui_elements(self):
        """Connect UI signals"""
        file_buttons = {
            'customSilicaFile_pushButton': SampleType.SILICA,
            'customSolventFile_pushButton': SampleType.SOLVENT,
            'customSampleFile_pushButton': SampleType.SAMPLE,
        }
        
        for btn_name, sample_type in file_buttons.items():
            if hasattr(self, btn_name):
                getattr(self, btn_name).clicked.connect(
                    lambda checked, st=sample_type: self.load_data_file(st)
                )
        
        # Fit buttons
        fit_buttons = {
            'silicaCA_fit_pushButton': (SampleType.SILICA, 'CA'),
            'silicaOA_fit_pushButton': (SampleType.SILICA, 'OA'),
            'solventCA_fit_pushButton': (SampleType.SOLVENT, 'CA'),
            'solventOA_fit_pushButton': (SampleType.SOLVENT, 'OA'),
            'sampleCA_fit_pushButton': (SampleType.SAMPLE, 'CA'),
            'sampleOA_fit_pushButton': (SampleType.SAMPLE, 'OA'),
        }
        
        for btn_name, (sample_type, scan_type) in fit_buttons.items():
            if hasattr(self, btn_name):
                getattr(self, btn_name).clicked.connect(
                    lambda checked, st=sample_type, sc=scan_type: self.auto_fit(st, sc)
                )
    
    def load_data_file(self, sample_type: SampleType):
        """Load a Z-scan data file"""
        filepath, _ = QFileDialog.getOpenFileName(
            self, f'Open {sample_type.value.title()} Data', '',
            'Text Files (*.txt);;CSV Files (*.csv)'
        )
        
        if not filepath:
            return
        
        file_data = ZScanDataParser.parse(filepath)
        if not file_data:
            self.statusbar.showMessage(f'Error: Could not parse {filepath}')
            return
        
        self.file_data[sample_type] = file_data
        self.processed_data[sample_type] = ZScanProcessor.normalize(file_data)
        
        # Update UI
        filename_widget_name = f'{sample_type.value}Filename_lineEdit'
        if hasattr(self, filename_widget_name):
            getattr(self, filename_widget_name).setText(Path(filepath).name)
        
        # Update wavelength
        if sample_type == SampleType.SILICA and hasattr(self, 'wavelength_dataFittingTab_doubleSpinBox'):
            self.wavelength_dataFittingTab_doubleSpinBox.setValue(file_data.wavelength)
        
        self.statusbar.showMessage(
            f'Loaded {len(self.processed_data[sample_type]["ca"])} points: {Path(filepath).name}'
        )
        
        # Enable fit buttons
        for scan_type in ['CA', 'OA']:
            btn_name = f'{sample_type.value}{scan_type}_fit_pushButton'
            if hasattr(self, btn_name):
                getattr(self, btn_name).setEnabled(True)
        
        # Plot initial data
        self.plot_data(sample_type, 'CA')
    
    def plot_data(self, sample_type: SampleType, scan_type: str):
        """Plot data with fit if available"""
        data = self.processed_data[sample_type]
        if data is None:
            return
        
        canvas_key = f'{sample_type.value}_{scan_type.lower()}'
        canvas = self.canvases.get(canvas_key)
        
        if canvas is None:
            return
        
        x = data['position_centered']
        
        if scan_type == 'CA':
            y = data['ca_antisym']
            ylabel = 'ΔT/T₀ (Normalized)'
            title = f'{sample_type.value.title()} CA (Antisymmetric)'
        else:
            y = data['oa']
            ylabel = '1 - T (Normalized)'
            title = f'{sample_type.value.title()} OA'
        
        result_key = (sample_type, scan_type)
        result = self.fitting_results.get(result_key)
        
        if result:
            canvas.plot(x, y, y_fit=result.fit_curve, residuals=result.residuals,
                       title=title, xlabel='Position (mm)', ylabel=ylabel)
        else:
            canvas.plot(x, y, title=title, xlabel='Position (mm)', ylabel=ylabel)
    
    def auto_fit(self, sample_type: SampleType, scan_type: str):
        """Run automatic fitting"""
        data = self.processed_data[sample_type]
        file_data = self.file_data[sample_type]
        
        if data is None or file_data is None:
            self.statusbar.showMessage(f'Error: No {sample_type.value} data loaded')
            return
        
        self.statusbar.showMessage(f'Auto-fitting {sample_type.value} {scan_type}...')
        
        # Infer initial parameters
        initial_params = ZScanPhysicsEngine.infer_initial_params(
            data, file_data.wavelength
        )
        
        # Run fit in background
        if scan_type == 'CA':
            fit_func = ZScanPhysicsEngine.auto_fit_ca
        else:
            fit_func = ZScanPhysicsEngine.auto_fit_oa
        
        self.fit_thread = QThread()
        self.fit_worker = FitWorker(
            fit_func,
            data=data,
            wavelength_nm=file_data.wavelength,
            initial_params=initial_params,
        )
        
        self.fit_worker.moveToThread(self.fit_thread)
        self.fit_thread.started.connect(self.fit_worker.run)
        self.fit_worker.finished.connect(
            lambda result: self.on_fit_complete(sample_type, scan_type, result)
        )
        self.fit_worker.error.connect(self.on_fit_error)
        
        self.fit_thread.start()
    
    def on_fit_complete(self, sample_type: SampleType, scan_type: str, result: FittingResult):
        """Handle fit completion"""
        self.fitting_results[(sample_type, scan_type)] = result
        
        self.plot_data(sample_type, scan_type)
        
        # Update summary displays
        self.update_fit_summary(sample_type, scan_type, result)
        
        self.statusbar.showMessage(
            f'{sample_type.value.title()} {scan_type} fitted | R²={result.r_squared:.4f} | χ²={result.chi_squared:.6f}'
        )
        
        self.fit_thread.quit()
        self.fit_thread.wait()
    
    def on_fit_error(self, error_msg: str):
        """Handle fit error"""
        self.statusbar.showMessage(f'Fit error: {error_msg}')
        self.fit_thread.quit()
        self.fit_thread.wait()
    
    def update_fit_summary(self, sample_type: SampleType, scan_type: str, result: FittingResult):
        """Update UI with fit results"""
        param_prefix = f'{sample_type.value}{scan_type}'
        
        # Update amplitude
        amp_widget_name = f'{param_prefix}_deltaPhi0Summary_doubleSpinBox' if scan_type == 'CA' else f'{param_prefix}_TSummary_doubleSpinBox'
        if hasattr(self, amp_widget_name):
            getattr(self, amp_widget_name).setValue(result.params.amplitude)
        
        # Update beamwaist
        bw_widget_name = f'{param_prefix}_beamwaistSummary_doubleSpinBox'
        if hasattr(self, bw_widget_name):
            getattr(self, bw_widget_name).setValue(result.params.beamwaist * 1e6)
        
        # Update zero level
        zl_widget_name = f'{param_prefix}_zeroLevel_doubleSpinBox'
        if hasattr(self, zl_widget_name):
            getattr(self, zl_widget_name).setValue(result.params.zero_level)
        
        # Calculate Rayleigh range
        z0 = result.params.beamwaist / np.sqrt(2)
        rl_widget_name = f'{param_prefix}_rayleighRangeSummary_doubleSpinBox'
        if hasattr(self, rl_widget_name):
            getattr(self, rl_widget_name).setValue(z0 * 1000)


def main():
    app = QApplication(sys.argv)
    
    ui_path = Path(__file__).parent / "window.ui"
    
    if ui_path.exists():
        window = MainWindow(str(ui_path))
    else:
        print(f"UI file not found, using minimal UI")
        window = MainWindow(None)
    
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
