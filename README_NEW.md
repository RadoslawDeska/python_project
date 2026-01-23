# Z-Scan Measurement Analysis Application

A comprehensive Python application for automated Z-scan measurements and data analysis using ultrafast lasers, National Instruments DAQ hardware, and Thorlabs motion control components.

## Project Overview

Z-scan is a powerful technique for measuring nonlinear optical properties of materials. This application provides:

- **Real-time data acquisition** from NI DAQ devices
- **Motorized stage control** for precise positioning
- **Automated curve fitting** using Sheik-Bahae theory
- **Multi-sample analysis** (reference silica, solvent baseline, and unknown samples)
- **Closed-aperture (CA) and open-aperture (OA) measurements**
- **Nonlinear parameter extraction** (n₂, β, absorption coefficients)
- **Interactive fitting** with real-time visualization

## Project Structure

```
python_project/
├── zscan1.py                 # Main PyQt5 GUI application
├── window.ui                 # Qt Designer UI definition
├── solvents.json             # Solvent properties database
├── data/                     # Z-scan measurement data files
├── lib/                      # Core library modules
│   ├── __init__.py
│   ├── config.py             # Configuration and constants
│   ├── integration.py        # Sheik-Bahae integration theory
│   ├── fitting.py            # Curve fitting engine
│   ├── sample_fitting.py     # Sample-specific fitting utilities
│   ├── scientific_rounding.py # Error propagation rounding
│   ├── figure.py             # Matplotlib canvas integration
│   ├── worker.py             # Threading support
│   ├── cursors.py            # Interactive plot cursors
│   ├── mgmotor.py            # Thorlabs motor interface
│   └── cursors.py            # Interactive plot tools
├── testing/                  # Test suite
│   └── test_integration_class.py
└── README.md
```

## Core Modules

### Main Application (`zscan1.py`)
- PyQt5-based graphical user interface
- Real-time measurement control and visualization
- Data acquisition integration with NI DAQ
- Motor control and positioning
- Multi-tab interface for different measurement types

### Physics & Mathematics (`lib/integration.py`)
Implements the Sheik-Bahae theory for Z-scan measurements:
- Electric field propagation calculation
- Gaussian decomposition (up to 8 components)
- Transmittance calculation for CA and OA geometries
- Nonlinear absorption models (2PA, 3PA, RSA, SA, combinations)

### Curve Fitting (`lib/fitting.py`)
- Manual fitting with parameter sliders
- Automatic optimization using least-squares
- Weighted fitting based on cursor regions of interest
- Parameter uncertainty estimation via error propagation

### Sample Fitting (`lib/sample_fitting.py`)
Complete workflow for unknown sample analysis:
- Prerequisite validation (silica + solvent must be fitted first)
- Automatic and manual CA fitting
- Automatic and manual OA fitting
- Parameter extraction and organization

### Configuration (`lib/config.py`)
Centralized configuration including:
- Physical constants and beam parameters
- UI settings and color themes
- Hardware device parameters
- Fitting bounds and models
- File I/O settings

### Utilities
- `scientific_rounding.py`: Error-based rounding for publication-quality results
- `figure.py`: Matplotlib integration for real-time plotting
- `worker.py`: Threading for non-blocking GUI operations
- `cursors.py`: Interactive plot tools for data selection
- `mgmotor.py`: Thorlabs motion control interface

## Key Features

### Multi-Measurement Support
- **Silica Reference**: Calibration standard for extracting absolute n₂ and β
- **Solvent Baseline**: Reference for sample measurements
- **Unknown Samples**: Full characterization using solvent background subtraction

### Flexible Fitting Options
- **Closed Aperture (CA)**: Measures refractive nonlinearity (n₂)
- **Open Aperture (OA)**: Measures absorptive nonlinearity (β)
- Both CA and OA can be fitted manually or automatically

### Interactive Controls
- Real-time parameter adjustment via sliders
- Cursor-based region-of-interest selection
- Live curve preview during manual fitting
- Automatic fitting with cursor weighting

### Hardware Integration
- National Instruments DAQ for multi-channel data acquisition
- Thorlabs MG17 stepper motor for Z-positioning
- Multi-channel signal detection
- External trigger synchronization with laser

## Measurement Workflow

### Basic Procedure
1. **Setup Phase**
   - Load solvent properties from database
   - Configure measurement parameters (wavelength, z-range, etc.)
   - Set up hardware and establish connections

2. **Reference Measurement (Silica)**
   - Perform Z-scan on reference silica sample
   - Fit CA data to extract beam profile
   - Results used to normalize subsequent measurements

3. **Solvent Baseline**
   - Measure pure solvent (no sample)
   - Fit CA data to characterize solvent nonlinearity
   - CA and OA data can be analyzed

4. **Sample Measurement**
   - Insert sample into cuvette
   - Perform complete Z-scan (CA and OA)
   - Subtract solvent baseline
   - Extract sample nonlinear optical properties

### Data Processing
- Real-time averaging and noise reduction (median filtering)
- Automatic data normalization
- Error calculation via error propagation
- Results saved with full parameter set

## Configuration Guide

### Key Constants (in `lib/config.py`)
- `N_COMPONENTS = 8`: Gaussian decomposition order
- `INTEGRATION_STEPS = 30`: Radial integration resolution
- `CUVETTE_PATH_LENGTH = 0.001` m: Cuvette thickness
- `DAQ_SAMPLING_RATE = 1000` Hz: Acquisition rate

### Fitting Parameters
CA fitting bounds (see config.py):
- DPhi0: [-2, 2] rad (phase shift)
- Beamwaist: [15 μm, 150 μm]
- Zero level: [0.75, 1.25]
- Center: [-50, 50] points

### Nonlinear Absorption Models
Supported OA models:
- 2PA: 2-photon absorption
- 3PA: 3-photon absorption
- 2PA+3PA: Combined 2PA and 3PA
- RSA: Reverse saturable absorption
- SA: Saturable absorption
- 2PA+SA: Combined 2PA and SA

## Usage

### Running the Application
```bash
python zscan1.py
```

### Data Files
- Data is stored in `data/` directory
- File naming: `YYYY_MM_DD__HH_MM__samplename_params.txt`
- Headers contain measurement parameters for automatic loading

### Solvent Database
`solvents.json` contains:
```json
{
  "SolventName": {
    "density": 0.XXX,  // g/cm³
    "index": 1.XXXX    // refractive index at wavelength
  }
}
```

## Dependencies

### Core Requirements
- PyQt5 >= 5.15
- numpy >= 1.20
- scipy >= 1.7
- matplotlib >= 3.3
- lmfit >= 0.21
- nidaqmx >= 0.8 (for NI DAQ)

### Optional
- pytest (for testing)
- black (code formatting)
- pylint (code analysis)

## Architecture Notes

### Design Patterns
- **Model-View-Controller**: Separation of UI and logic
- **Factory Pattern**: Integration and Fitting object creation
- **Thread Pool**: Non-blocking hardware operations
- **Observer Pattern**: Signal-slot connections for UI updates

### Threading
- Main thread: UI responsiveness
- Worker threads: Long-running operations (motor movement, fitting)
- Signal-based communication between threads

### Error Handling
- Comprehensive logging throughout
- User-friendly error dialogs
- Hardware connection validation
- Parameter bounds checking

## Key Improvements Made

### Code Organization
1. **Separated concerns**: Physics, fitting, and UI logic now in distinct modules
2. **Configuration centralization**: All constants moved to `lib/config.py`
3. **Module documentation**: Comprehensive docstrings and module-level documentation
4. **Type hints**: Added type annotations throughout for clarity

### Completed Features
1. **Integration module** (`lib/integration.py`): Well-documented Sheik-Bahae implementation
2. **Fitting module** (`lib/fitting.py`): Unified manual and automatic fitting
3. **Sample fitting** (`lib/sample_fitting.py`): Complete sample measurement workflow
4. **Configuration** (`lib/config.py`): Centralized settings and constants

### Improved Maintainability
- Reduced main file complexity through modularization
- Consistent error handling patterns
- Clear separation between hardware and analysis code
- Reusable components for future extensions

## Future Improvements

Areas for enhancement:

### User Interface
- [ ] Progress indicators for fitting operations
- [ ] Real-time fitting statistics display
- [ ] Batch processing interface
- [ ] Results export tools

### Analysis
- [ ] Multiple curve comparison tools
- [ ] Statistical analysis suite
- [ ] Automated report generation
- [ ] Advanced data visualization

### Performance
- [ ] GPU acceleration for integration calculations
- [ ] Caching of computed curves
- [ ] Parallel fitting for multiple datasets
- [ ] Optimized numerical algorithms

### Testing
- [ ] Extended unit test coverage
- [ ] Integration tests for full workflows
- [ ] Hardware mock interfaces for CI/CD
- [ ] Validation against literature examples

## Testing

Run tests with:
```bash
pytest testing/
```

Current test coverage:
- Integration class basic functionality
- Parameter extraction and rounding

## References

Key theoretical references:
- Sheik-Bahae et al., IEEE J. Quantum Electron. 26, 760-769 (1990)
- Said et al., J. Opt. Soc. Am. B 9, 405-414 (1992)
- Xia et al., Opt. Lett. 19, 1849-1851 (1994)
