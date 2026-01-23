# tests/test_validation.py
import pytest
from lib.validation import ParameterValidator


class TestParameterValidator:
    
    def test_valid_general_parameters(self):
        valid, msg = ParameterValidator.validate_general_parameters(
            wavelength=800, silica_thickness=2.0, aperture_diameter=100, z_scan_range=5
        )
        assert valid is True
        assert msg == ""
    
    def test_negative_wavelength_fails(self):
        valid, msg = ParameterValidator.validate_general_parameters(
            wavelength=-800, silica_thickness=2.0, aperture_diameter=100, z_scan_range=5
        )
        assert valid is False
        assert "Wavelength" in msg
    
    def test_zero_silica_thickness_fails(self):
        valid, msg = ParameterValidator.validate_general_parameters(
            wavelength=800, silica_thickness=0, aperture_diameter=100, z_scan_range=5
        )
        assert valid is False
        assert "Silica thickness" in msg
    
    def test_valid_fitting_parameters(self):
        valid, msg = ParameterValidator.validate_fitting_parameters(
            dphi0=0.5, beamwaist=50e-6, zero_level=1.0, center_point=0
        )
        assert valid is True
    
    def test_zero_level_too_low_fails(self):
        valid, msg = ParameterValidator.validate_fitting_parameters(
            dphi0=0.5, beamwaist=50e-6, zero_level=0.3, center_point=0
        )
        assert valid is False
        assert "Zero level" in msg