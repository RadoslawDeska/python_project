# tests/conftest.py
"""
Shared pytest fixtures and configuration for all tests.
"""

import pytest
import numpy as np
from unittest.mock import Mock
from numpy.typing import NDArray


@pytest.fixture
def sample_z_positions() -> NDArray:
    """Create sample z-scan position array."""
    return np.linspace(-2, 2, 41)  # -2 to +2 mm in 41 steps


@pytest.fixture
def sample_ca_data(sample_z_positions: NDArray) -> NDArray:
    """Create realistic closed aperture data."""
    # Gaussian-like peak in center
    z = sample_z_positions
    return 1.0 + 0.1 * np.exp(-(z ** 2) / 0.5)


@pytest.fixture
def sample_oa_data(sample_z_positions: NDArray) -> NDArray:
    """Create realistic open aperture data."""
    # Dip at center (absorption)
    z = sample_z_positions
    return 0.9 - 0.05 * np.exp(-(z ** 2) / 1.0)


@pytest.fixture
def mock_fitting_manager():
    """Mock FittingManager for testing."""
    mock = Mock()
    mock.fit_ca_automatic = Mock(return_value=({"params": {}}, np.array([])))
    mock.fit_oa_automatic = Mock(return_value=({"params": {}}, np.array([])))
    return mock


@pytest.fixture
def mock_data_manager():
    """Mock DataManager for testing."""
    mock = Mock()
    mock.data = {
        "positions": [],
        "absolute": [[], [], []],
        "relative": [[], [], []]
    }
    mock.rms_value = 0.0
    return mock


@pytest.fixture
def mock_state_manager():
    """Mock StateManager for testing."""
    mock = Mock()
    mock.running = False
    mock.initialized = True
    mock.silica_autofit_done = False
    return mock