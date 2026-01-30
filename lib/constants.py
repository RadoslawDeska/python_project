"""
Central repository for all magic numbers and constants used in the application.
This makes code more maintainable and easier to test.
"""

from typing import Dict
# ============================================================================
# PHYSICAL CONSTANTS AND Z-SCAN PARAMETERS
# ============================================================================

SILICA_BETA = 0
"""Non-linear absorption coefficient for silica reference (normalized)"""

N_COMPONENTS = 8
"""Number of Gaussian field components for decomposition"""

INTEGRATION_STEPS = 30
"""Number of radial integration steps for transmittance calculation"""

CUVETTE_PATH_LENGTH = 0.001  # meters
"""Path length through the cuvette for transmission"""

MAX_DPHI0 = 3.142
"""Maximum on-axis phase shift for silica (used for slider limits)"""

SOLVENT_T_SLIDER_MAX = 1
"""Maximum transmittance value for solvent OA fitting slider"""


# Slider Constants
SOLVENT_CA_SLIDER_CENTER_OFFSET = 50
SOLVENT_OA_SLIDER_CENTER_OFFSET = 50
SOLVENT_OA_SLIDER_MIN = 0
SOLVENT_OA_SLIDER_MAX = 100
SAMPLE_OA_SLIDER_CENTER_OFFSET = 50

# Fitting Parameter Scales
DPHI0_SCALE_FACTOR = 1.0 / 3.142  # DPhi0 / π
BEAMWAIST_UM_SCALE = 1e6  # Convert meters to micrometers
RAYLEIGH_LENGTH_MM_SCALE = 1e3  # Convert meters to millimeters
LASER_INTENSITY_GW_SCALE = 1e-13  # Convert to GW/cm²

# Timing Constants
TIMER_UPDATE_INTERVAL = 100  # milliseconds

# ============================================================================
# DATA ACQUISITION PARAMETERS
# ============================================================================

DAQ_SAMPLING_RATE = 1000  # Hz (1 kHz, matching laser repetition rate)
"""National Instruments DAQ sampling rate"""

DAQ_TIMING_WAIT = 0.2  # seconds
"""Time to wait for DAQ to collect all samples"""

# ============================================================================
# MOTOR CONTROL PARAMETERS
# ============================================================================

MOTOR_INIT_WAIT = 0.2  # seconds
"""Time to wait after motor initialization"""

MOTOR_HOMING_WAIT = 0.2  # seconds
"""Time to wait after motor homing completes"""

MOTOR_POSITION_MAX = 100  # mm (arbitrary limit for custom positions)
"""Maximum allowed custom motor position"""

# ============================================================================
# AUDIO FEEDBACK
# ============================================================================

COMPLETION_SOUND_FREQUENCY = 800  # Hz
"""Frequency of completion beep"""

COMPLETION_SOUND_DURATION = 500  # milliseconds
"""Duration of completion beep"""

COMPLETION_SOUND_REPEATS = 3
"""Number of times to repeat completion beep"""

COMPLETION_SOUND_INTERVAL = 0.05  # seconds
"""Time between beeps"""


# Normalization
RMS_NOISE_EDGE_WINDOW = 10  # Use first/last 10 points for normalization
DATA_NORMALIZATION_WINDOW = 10


# Default values
DEFAULT_BEAMWAIST = 50e-6  # meters
DEFAULT_ZERO_LEVEL = 1.0
DEFAULT_CENTER_POINT = 0

# ============================================================================
# UI AND DISPLAY SETTINGS
# ============================================================================

PLOT_PADDING_VERTICAL = 0.01
"""Default matplotlib vertical padding for plots"""

PLOT_FONT_SIZE = 7
"""Default matplotlib font size for plots"""

# Measurement chart types
CHART_TYPES = {"absolute": "Absolute Signal", "relative": "Relative Signal"}

# Sample types
SAMPLE_TYPES = {
    "Silica": "Reference Silica",
    "Solvent": "Solvent Baseline",
    "Sample": "Sample Under Study",
}

# Measurement types
MEASUREMENT_TYPES = {"CA": "Closed Aperture", "OA": "Open Aperture"}

# Unit converters
POSITIONS_MM_SCALE = 1e3  # Convert meters to millimeters

# ============================================================================
# FILE PATHS AND DEFAULTS
# ============================================================================

DEFAULT_DATA_FOLDER = "data"
"""Default folder for data files"""

# ============================================================================
# FITTING PARAMETERS AND LIMITS
# ============================================================================

# Closed aperture (CA) fitting bounds
CA_FITTING_PARAMS: Dict[str, Dict[str, float]] = {
    "Zero": {"min": 0.5, "max": 1.5},
    "Center": {"min": -50, "max": 50},
    "DPhi0": {"min": -10, "max": 10},
    "Beamwaist": {"min": 1e-6, "max": 500e-6},
}

# Open aperture (OA) fitting bounds
OA_FITTING_PARAMS: Dict[str, Dict[str, float]] = {
    "Zero": {"min": 0.75, "max": 1.25},
    "Center": {"min": -50, "max": 50},
    "T": {"min": -10, "max": 10},
    "Beamwaist": {"min": 1e-6, "max": 500e-6},
}

MAX_FITTING_ITERATIONS = 1000
"""Maximum iterations for least-squares optimization"""

# ============================================================================
# ABSORPTION MODELS FOR OA FITTING
# ============================================================================

OA_ABSORPTION_MODELS = [
    "2PA",  # 2-photon absorption
    "3PA",  # 3-photon absorption
    "2PA+3PA",  # Combined 2-photon and 3-photon
    "RSA",  # Reverse saturable absorption
    "SA",  # Saturable absorption
    "2PA+SA",  # 2-photon + saturable absorption
]

# ============================================================================
# COLOR PALETTE (Dark Theme)
# ============================================================================

DARK_THEME = {
    "window_bg": (35, 35, 40),
    "window_text": (200, 200, 200),
    "base": (60, 60, 65),
    "text": (200, 200, 200),
    "button": (35, 35, 40),
    "button_text": (200, 200, 200),
    "highlight": (42, 130, 218),
    "disabled_text": (100, 100, 100),
    "grid": (60, 60, 65),
    "spine": (200, 200, 200),
}

# ============================================================================
# COLOR PALETTE (Light Theme)
# ============================================================================

LIGHT_THEME = {
    "window_bg": (255, 255, 255),
    "window_text": (0, 0, 0),
    "base": (255, 255, 255),
    "text": (0, 0, 0),
    "button": (240, 240, 240),
    "button_text": (0, 0, 0),
    "highlight": (42, 130, 218),
    "disabled_text": (100, 100, 100),
    "grid": (192, 192, 192),
    "spine": (0, 0, 0),
}

# ============================================================================
# PLOT CONFIGURATION
# ============================================================================

PLOT_FIGURE_SIZE = (5, 4)
"""Default matplotlib figure size (width, height) in inches"""

PLOT_DPI = 100
"""Default matplotlib DPI"""

PLOT_AXES_RECT = [0.07, 0.185, 0.922, 0.78]
"""Axes position [left, bottom, width, height] in figure fraction"""

# ============================================================================
# REGULAR EXPRESSIONS FOR FILE PARSING
# ============================================================================

# Pattern for extracting numeric values from file headers
NUMERIC_PATTERN = r"(([1-9][0-9]*\.?[0-9]*)|(\.[0-9]+))([Ee][+-]?[0-9]+)?"
"""Regex pattern for matching numeric values"""

CONCENTRATION_PATTERN = (
    r"(([0-9][0-9]*\.?[0-9]*\s*%)|(\.[0-9]()\s*%+))([Ee][+-]?[0-9]()\s*%+)?"
)
"""Regex pattern for concentration values with % symbol"""
