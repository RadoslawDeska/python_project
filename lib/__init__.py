"""
Z-Scan library package.

Core modules for Z-scan data acquisition and analysis:

Physics & Fitting:
- integration: Sheik-Bahae theoretical framework
- fitting: Automatic and manual curve fitting
- sample_fitting: Sample measurement workflow

Managers (for modular architecture):
- state_manager: Application state management
- data_manager: Data storage and calculations
- ui_manager: UI initialization and styling
- chart_manager: Chart and visualization management
- fitting_manager: Fitting workflow management
- measurement_manager: Hardware measurement and motor control

Utilities:
- config: Centralized configuration
- scientific_rounding: Error-based rounding
- figure: Matplotlib integration
- worker: Threading utilities
- cursors: Interactive plot tools
"""

from lib.integration import Integration
from lib.fitting import Fitting
from lib.sample_fitting import SampleFitter
from lib.scientific_rounding import error_rounding
from lib.state_manager import StateManager
from lib.data_manager import DataManager
from lib.ui_manager import UIManager
from lib.chart_manager import ChartManager
from lib.fitting_manager import FittingManager
from lib.measurement_manager import MeasurementManager

__all__ = [
    # Physics & Fitting
    'Integration',
    'Fitting',
    'SampleFitter',
    'error_rounding',
    # Managers
    'StateManager',
    'DataManager',
    'UIManager',
    'ChartManager',
    'FittingManager',
    'MeasurementManager',
]

__version__ = '0.3.0'
__author__ = 'Z-Scan Development Team'
