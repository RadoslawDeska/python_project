# lib/validation.py
"""
Parameter validation to catch configuration errors early.

NOTE: This validator checks NUMERICAL VALIDITY and MEASUREMENT CONSTRAINTS,
not physics validity. Physical constraints depend on material properties
(e.g., n2 sign determines DPhi0 sign), so we allow a wide range except for
well-known materials like silica.
"""

from typing import Tuple


class ParameterValidator:
    """Validates Z-scan measurement and fitting parameters.
    
    Enforces numerical validity and measurement constraints, but NOT generic physics 
    constraints that depend on material properties.
    
    For well-known materials like silica with known physics (n2 > 0 in VIS/NIR),
    use material-specific validators like validate_silica_fitting_parameters().
    """
    
    @staticmethod
    def validate_general_parameters(
        wavelength: float,
        silica_thickness: float,
        aperture_diameter: float,
        z_scan_range: float
    ) -> Tuple[bool, str]:
        """Validate general measurement parameters.
        
        Checks numerical validity and experimental reasonableness.
        
        Args:
            wavelength: Laser wavelength in nm
            silica_thickness: Silica reference thickness in mm
            aperture_diameter: Aperture diameter in mm
            z_scan_range: Z-scan range in mm
            
        Returns:
            (is_valid, error_message)
        """
        if not isinstance(wavelength, (int, float)):
            return False, "Wavelength must be a number"
        if wavelength <= 0:
            return False, "Wavelength must be positive (>0 nm)"
        if wavelength > 2000:  # Reasonable upper limit for optical wavelengths
            return False, "Wavelength seems unrealistic (>2000 nm)"
        
        if not isinstance(silica_thickness, (int, float)):
            return False, "Silica thickness must be a number"
        if silica_thickness <= 0:
            return False, "Silica thickness must be positive (>0 mm)"
        if silica_thickness > 10:  # Reasonable upper limit
            return False, "Silica thickness seems unrealistic (>10 mm)"
        
        if not isinstance(aperture_diameter, (int, float)):
            return False, "Aperture diameter must be a number"
        if aperture_diameter <= 0:
            return False, "Aperture diameter must be positive (>0 mm)"
        
        if not isinstance(z_scan_range, (int, float)):
            return False, "Z-scan range must be a number"
        if z_scan_range <= 0:
            return False, "Z-scan range must be positive (>0 mm)"
        
        return True, ""
    
    @staticmethod
    def validate_silica_fitting_parameters(
        dphi0: float,
        beamwaist: float,
        zero_level: float,
        center_point: float
    ) -> Tuple[bool, str]:
        """Validate fitted parameters for SILICA specifically.
        
        Silica is a reference material with well-known physics:
        - In visible and NIR range, n2 > 0, so DPhi0 MUST be positive
        - For other materials, use validate_fitting_parameters() instead
        
        Args:
            dphi0: Phase shift parameter (MUST be positive for silica)
            beamwaist: Beam waist in meters
            zero_level: Baseline transmittance
            center_point: Center position in data points
            
        Returns:
            (is_valid, error_message)
        """
        # For silica in VIS/NIR, n2 > 0, so DPhi0 must be positive
        if not isinstance(dphi0, (int, float)):
            return False, "DPhi0 must be a number"
        if dphi0 <= 0:
            return False, f"DPhi0 must be positive for silica (n2 > 0 in VIS/NIR), got {dphi0}"
        
        if abs(dphi0) > 100:
            return False, f"DPhi0 magnitude seems unrealistic: {abs(dphi0)}"
        
        # Beamwaist must be physically positive and in reasonable range
        if not isinstance(beamwaist, (int, float)):
            return False, "Beamwaist must be a number"
        if beamwaist <= 0:
            return False, "Beamwaist must be positive"
        if beamwaist > 1e-4:  # >100 µm
            return False, "Beamwaist seems too large (unrealistic value >100 µm)"
        if beamwaist < 1e-7:  # <0.1 µm
            return False, "Beamwaist seems too small (unrealistic value <0.1 µm)"
        
        # Zero level is typically near 1, but allow reasonable variations
        if not isinstance(zero_level, (int, float)):
            return False, "Zero level must be a number"
        if zero_level <= 0:
            return False, "Zero level must be positive"
        if zero_level < 0.5 or zero_level > 2.0:
            return False, f"Zero level out of expected range [0.5, 2.0]: {zero_level}"
        
        # Center point should be within reasonable data range
        if not isinstance(center_point, (int, float)):
            return False, "Center point must be a number"
        
        return True, ""
    
    @staticmethod
    def validate_fitting_parameters(
        dphi0: float,
        beamwaist: float,
        zero_level: float,
        center_point: float
    ) -> Tuple[bool, str]:
        """Validate fitted parameter ranges for generic samples (solvent/sample).
        
        Checks numerical validity and reasonable bounds. Does NOT enforce
        physics constraints like sign of DPhi0 (which depends on material n2).
        For silica-specific validation, use validate_silica_fitting_parameters().
        
        Args:
            dphi0: Phase shift parameter (can be positive OR negative, material-dependent)
            beamwaist: Beam waist in meters
            zero_level: Baseline transmittance
            center_point: Center position in data points
            
        Returns:
            (is_valid, error_message)
        """
        # DPhi0 can be positive or negative (material-dependent)
        # but must be numerically reasonable (not infinite or NaN)
        if not isinstance(dphi0, (int, float)):
            return False, "DPhi0 must be a number"
        if abs(dphi0) > 100:  # Reasonable upper bound for Z-scan
            return False, f"DPhi0 magnitude seems unrealistic: {abs(dphi0)} (expected <100)"
        
        # Beamwaist must be physically positive and in reasonable range
        if not isinstance(beamwaist, (int, float)):
            return False, "Beamwaist must be a number"
        if beamwaist <= 0:
            return False, "Beamwaist must be positive"
        if beamwaist > 1e-4:  # >100 µm
            return False, "Beamwaist seems too large (unrealistic value >100 µm)"
        if beamwaist < 1e-7:  # <0.1 µm
            return False, "Beamwaist seems too small (unrealistic value <0.1 µm)"
        
        # Zero level is typically near 1, but allow reasonable variations
        if not isinstance(zero_level, (int, float)):
            return False, "Zero level must be a number"
        if zero_level <= 0:
            return False, "Zero level must be positive"
        if zero_level < 0.5 or zero_level > 2.0:
            return False, f"Zero level out of expected range [0.5, 2.0]: {zero_level}"
        
        # Center point should be within reasonable data range
        if not isinstance(center_point, (int, float)):
            return False, "Center point must be a number"
        
        return True, ""
    
    @staticmethod
    def validate_oa_parameters(t_value: float) -> Tuple[bool, str]:
        """Validate Open Aperture amplitude/transmittance parameter.
        
        T (transmittance) can be positive or negative in general calculations
        depending on the absorption model and material. Only check numerical
        validity and reasonable bounds.
        
        Args:
            t_value: Transmittance/amplitude parameter
            
        Returns:
            (is_valid, error_message)
        """
        if not isinstance(t_value, (int, float)):
            return False, "T value must be a number"
        
        # Check for reasonable bounds (not infinite or extremely large)
        # Allow both positive and negative values
        if abs(t_value) > 10:
            return False, f"OA T parameter magnitude seems unrealistic: {abs(t_value)}"
        
        return True, ""
    
    @staticmethod
    def validate_beamwaist(beamwaist: float) -> Tuple[bool, str]:
        """Validate beam waist radius.
        
        Beamwaist must be positive and in optical range.
        
        Args:
            beamwaist: Beam waist in meters
            
        Returns:
            (is_valid, error_message)
        """
        if not isinstance(beamwaist, (int, float)):
            return False, "Beamwaist must be a number"
        if beamwaist <= 0:
            return False, "Beamwaist must be positive (>0 m)"
        if beamwaist > 200e-6:  # >200 µm is unrealistic for focused beam
            return False, f"Beamwaist seems unrealistic: {beamwaist*1e6:.1f} µm (expected <200 µm)"
        if beamwaist < 1e-7:  # <0.1 µm
            return False, f"Beamwaist seems too small: {beamwaist*1e6:.4f} µm (expected >0.1 µm)"
        return True, ""
    
    @staticmethod
    def validate_zero_level(zero_level: float) -> Tuple[bool, str]:
        """Validate baseline/zero level.
        
        Zero level (transmittance at infinity) should be positive and near 1,
        but allow variations depending on setup.
        
        Args:
            zero_level: Normalized transmittance at infinity
            
        Returns:
            (is_valid, error_message)
        """
        if not isinstance(zero_level, (int, float)):
            return False, "Zero level must be a number"
        if zero_level <= 0:
            return False, "Zero level must be positive"
        # Allow wider range: some experimental setups may have different baselines
        if zero_level < 0.3 or zero_level > 3.0:
            return False, f"Zero level out of expected range [0.3, 3.0]: {zero_level}"
        return True, ""