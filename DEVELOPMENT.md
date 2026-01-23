"""
DEVELOPMENT.md - Development and Integration Guide for Z-Scan Application

This document provides guidance for developers working on the Z-Scan application,
including integration of new modules, testing procedures, and best practices.
"""

# Development and Integration Guide

## Overview of Recent Refactoring

The Z-Scan application has been reorganized from a monolithic structure (2700+ lines) into
a modular architecture with separated concerns:

### Previous Structure
- Single `zscan1.py` file (2704 lines) containing:
  - PyQt5 GUI code
  - Physics/mathematics
  - Data fitting algorithms
  - Threading logic
  - Hardware control
  - All business logic

### New Structure
- `zscan1.py` (main GUI - to be refactored)
- `lib/config.py` - Configuration and constants
- `lib/integration.py` - Sheik-Bahae theory
- `lib/fitting.py` - Curve fitting engine
- `lib/sample_fitting.py` - Sample measurement workflow
- `lib/scientific_rounding.py` - Error propagation
- `lib/figure.py` - Matplotlib integration
- `lib/worker.py` - Threading
- `lib/cursors.py` - Interactive tools
- `lib/mgmotor.py` - Motor control

## Integration Checklist

To fully integrate the refactored code into `zscan1.py`:

### Phase 1: Import Updates
- [ ] Add imports for new modules at top of zscan1.py:
  ```python
  from lib.config import *  # Import all configuration constants
  from lib.integration import Integration
  from lib.fitting import Fitting
  from lib.sample_fitting import SampleFitter
  ```

### Phase 2: Replace Embedded Classes
- [ ] Remove `Integration` class from zscan1.py (now in lib/integration.py)
- [ ] Remove `Fitting` class from zscan1.py (now in lib/fitting.py)
- [ ] Remove `MotorPositioner` class (keep in lib/mgmotor.py)
- [ ] Remove all constants from top of zscan1.py (use lib/config.py)

### Phase 3: Update References to Constants
Replace all hardcoded constants with config imports:
- [ ] `SILICA_BETA` → `from lib.config import SILICA_BETA`
- [ ] `N_COMPONENTS` → `from lib.config import N_COMPONENTS`
- [ ] `INTEGRATION_STEPS` → `from lib.config import INTEGRATION_STEPS`
- [ ] And so on for all other constants

### Phase 4: Complete Sample Fitting
The sample fitting section currently has `pass` statements that need implementation:

#### Line 1074 (fit_automatically) - Sample case:
```python
case "Sample":
    # Implement using SampleFitter
    from lib.sample_fitting import SampleFitter
    fitter = SampleFitter()
    
    # Check prerequisites
    is_valid, msg = fitter.check_prerequisites(
        self.silica_autofit_done,
        self.solventCA_autofit_done
    )
    if not is_valid:
        self.showdialog('Info', msg)
        return
    
    # Fit data
    # ... implementation here
```

#### Lines involving Sample OA fitting:
- [ ] Uncomment and complete the sample OA fitting code
- [ ] Use SampleFitter.fit_automatically_oa() method
- [ ] Handle absorption model selection from UI

### Phase 5: Error Handling Improvements
- [ ] Add try-except blocks around critical sections
- [ ] Use logging module consistently
- [ ] Provide user-friendly error messages

### Phase 6: Documentation Updates
- [ ] Add method docstrings to Window class
- [ ] Document all signal-slot connections
- [ ] Add parameter descriptions

## Module Usage Examples

### Integration Module
```python
from lib.integration import Integration
import numpy as np

# Create integration object for closed aperture
positions = np.linspace(-5e-3, 5e-3, 50)  # -5 to +5 mm, 50 points
integration = Integration(
    beta=0,                    # No absorption (silica reference)
    n2=1e-15,                  # m^2/W
    DPhi0=0.5,                 # rad
    positions=positions,
    d0=0.3,                    # m
    aperture_radius=0.001,     # m
    wavelength=800e-9,         # m
    beamwaist=50e-6,           # m
    n_components=8,
    integration_steps=30,
    stype="CA"
)

# Get transmittance
T_ca = integration.closed_sum
T_oa = integration.open_sum
T_normalized = T_ca / T_oa
```

### Fitting Module
```python
from lib.fitting import Fitting
from lib.integration import Integration

# Create fitter
fitter = Fitting(
    sample_type=integration,
    amplitude=0.5,             # Initial DPhi0 or T
    beamwaist=50e-6,           # m
    zero_level=1.0,
    centerpoint=25,            # data point
    nop=50,                    # number of points
    data=measured_data
)

# Manual fitting (for interactive sliders)
result = fitter.manual(
    zero_level=1.0,
    centerpoint=25,
    amplitude=0.5,
    beamwaist=50e-6,
    z_range=0.01,              # m
    window=self,               # Main window reference
    stype="CA"
)

# Automatic fitting
minimizer, result = fitter.automatic(
    z_range=0.01,
    ftype="Sample",
    stype="CA",
    line_xydata=(x_data, y_data),
    window=self
)
```

### Sample Fitting Module
```python
from lib.sample_fitting import SampleFitter

fitter = SampleFitter()

# Check if prerequisites are met
is_valid, msg = fitter.check_prerequisites(
    silica_fitted=True,
    solvent_fitted=True
)

if is_valid:
    # Fit CA data
    minimizer, result = fitter.fit_automatically_ca(
        sample_curves=integration,
        initial_dphi0=0.5,
        initial_beamwaist=50e-6,
        initial_zero_level=1.0,
        initial_centerpoint=25,
        nop=50,
        z_range=0.01,
        line_data=(x_data, y_data),
        window=self
    )
    
    # Extract parameters
    params = fitter.extract_ca_parameters(minimizer)
    print(f"DPhi0 = {params['dphi0']['value']} ± {params['dphi0']['stderr']}")
```

## Configuration Management

### Modifying Constants
All constants are centralized in `lib/config.py`. To modify:

1. Locate the constant in the appropriate section
2. Update the value
3. Restart the application

### Adding New Constants
1. Add to appropriate section in config.py
2. Add documentation
3. Update imports in relevant modules
4. No changes needed to individual modules - they import from config

### Theme Customization
Edit `DARK_THEME` or `LIGHT_THEME` dicts in config.py:
```python
DARK_THEME = {
    'window_bg': (35, 35, 40),      # RGB tuple
    'highlight': (42, 130, 218),
    # ... more colors
}
```

## Testing Guidelines

### Test Structure
```python
# tests/test_module.py
import pytest
from lib.integration import Integration
import numpy as np

class TestIntegration:
    def setup_method(self):
        """Create test fixtures"""
        self.positions = np.linspace(-5e-3, 5e-3, 50)
    
    def test_ca_fitting_basic(self):
        """Test CA fitting with known parameters"""
        # Arrange
        integration = Integration(...)
        
        # Act
        result = integration.closed_sum
        
        # Assert
        assert result is not None
        assert len(result) == 50
```

### Running Tests
```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_integration.py

# Run with coverage
pytest --cov=lib tests/

# Verbose output
pytest -v tests/
```

## Best Practices

### 1. Error Handling
```python
try:
    result = some_calculation()
except ZeroDivisionError:
    logging.error("Division by zero in calculation")
    return None
except Exception as e:
    logging.error(f"Unexpected error: {str(e)}")
    self.showdialog('Error', 'An unexpected error occurred')
```

### 2. Type Hints
```python
def calculate_curve(
    parameters: dict,
    data: NDArray,
    z_range: float = 0.01
) -> NDArray:
    """Calculate fitted curve.
    
    Args:
        parameters: Fitting parameters dictionary
        data: Measured data points
        z_range: Total z-scan range in meters
    
    Returns:
        Calculated curve array
    """
    ...
```

### 3. Documentation
```python
def complex_function(arg1: float, arg2: int) -> tuple:
    """Short description.
    
    Longer description explaining the function's purpose,
    algorithm, and any important details.
    
    Args:
        arg1 (float): Description of arg1
        arg2 (int): Description of arg2
    
    Returns:
        tuple: (result1, result2) - Description of each
    
    Raises:
        ValueError: If arg1 is negative
        TypeError: If arg2 is not integer
    
    Examples:
        >>> result = complex_function(1.5, 3)
        >>> len(result)
        2
    """
    ...
```

### 4. Logging
```python
import logging

logger = logging.getLogger(__name__)

logger.debug("Detailed diagnostic information")
logger.info("General informational message")
logger.warning("Warning message")
logger.error("Error occurred")
logger.critical("Critical error")
```

### 5. Threading
Always use the `Worker` class for long operations:
```python
def perform_fitting(self):
    """Run fitting in background thread"""
    def fitting_work():
        result = complex_fitting_algorithm()
        return result
    
    worker = Worker(fitting_work)
    worker.signals.result.connect(self.on_fitting_complete)
    worker.signals.error.connect(self.on_fitting_error)
    self.threadpool.start(worker)
```

## File Organization Guidelines

### When to Create New Module
- If code exceeds 200 lines and has clear responsibility
- If code is reused in multiple places
- If code represents a distinct concept/component

### Module Structure Template
```python
"""
module_name.py - Brief description

Longer description of module purpose, key classes/functions,
and any important implementation details.

Classes:
    ClassName: Description
    AnotherClass: Description

Functions:
    function_name: Description
"""

from typing import Optional, List
import logging

logger = logging.getLogger(__name__)


class MainClass:
    """Description of class purpose and usage."""
    
    def __init__(self, param: int):
        """Initialize the class."""
        self.param = param
    
    def method(self) -> str:
        """Description of what method does.
        
        Returns:
            Description of return value
        """
        ...


def standalone_function(arg: str) -> bool:
    """Description.
    
    Args:
        arg: Description
    
    Returns:
        Description
    """
    ...
```

## Migration Path

To gradually integrate the refactored code without breaking existing functionality:

1. **Week 1**: Import new modules, keep old classes in zscan1.py
2. **Week 2**: Replace Integration class usage, test thoroughly
3. **Week 3**: Replace Fitting class usage, add sample fitting
4. **Week 4**: Complete sample OA fitting implementation
5. **Week 5**: Cleanup, testing, documentation

## Troubleshooting Common Issues

### Issue: Import circular dependency
**Solution**: Move shared utilities to a separate module, restructure imports

### Issue: GUI freezes during fitting
**Solution**: Ensure fitting runs in Worker thread, not main thread

### Issue: Numerical instability in integration
**Solution**: Increase INTEGRATION_STEPS, check parameter ranges

### Issue: Inconsistent results
**Solution**: Verify beam parameters, check data normalization

## Resources

- **NumPy Documentation**: https://numpy.org/doc/
- **SciPy Documentation**: https://docs.scipy.org/
- **lmfit Documentation**: https://lmfit.github.io/lmfit-py/
- **PyQt5 Documentation**: https://www.riverbankcomputing.com/static/Docs/PyQt5/
- **Z-Scan Theory**: Sheik-Bahae et al., IEEE J. Quantum Electron. 26, 760-769 (1990)
