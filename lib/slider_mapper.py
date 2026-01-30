"""
Encapsulates logic for mapping between slider positions and physical parameters.
This eliminates scattered, error-prone slider conversion code.
"""

from typing import Tuple


class ParameterSliderMapper:
    """Maps between slider position (0..max) and physical parameter range."""
    
    def __init__(self, param_min: float, param_max: float, slider_max: int):
        """
        Args:
            param_min: Minimum physical value
            param_max: Maximum physical value
            slider_max: Maximum slider position (usually 100)
        """
        self.param_min = param_min
        self.param_max = param_max
        self.slider_max = slider_max
        
        if slider_max <= 0:
            raise ValueError("slider_max must be positive")
        if param_min >= param_max:
            raise ValueError("param_min must be less than param_max")
    
    def slider_to_physical(self, slider_value: int) -> float:
        """Convert slider position to physical parameter value.
        
        Args:
            slider_value: Position on slider (0..slider_max)
            
        Returns:
            Physical value in range [param_min, param_max]
        """
        if not (0 <= slider_value <= self.slider_max):
            raise ValueError(
                f"Slider value {slider_value} out of range [0, {self.slider_max}]"
            )
        
        normalized = slider_value / self.slider_max
        return normalized * (self.param_max - self.param_min) + self.param_min
    
    def physical_to_slider(self, physical_value: float) -> int:
        """Convert physical parameter value to slider position.
        
        Args:
            physical_value: Value in range [param_min, param_max]
            
        Returns:
            Integer slider position (0..slider_max)
        """
        if not (self.param_min <= physical_value <= self.param_max):
            raise ValueError(
                f"Physical value {physical_value} out of range "
                f"[{self.param_min}, {self.param_max}]"
            )
        
        normalized = (physical_value - self.param_min) / (self.param_max - self.param_min)
        return int(round(normalized * self.slider_max))
    
    def get_range(self) -> Tuple[float, float]:
        """Get physical parameter range."""
        return self.param_min, self.param_max


class CenterPointSliderMapper(ParameterSliderMapper):
    """Special mapper for center point parameters that are symmetric around zero.
    
    Maps slider 0..100 to center -z_range/2..z_range/2
    The range depends on the actual z-scan range.
    """
    
    def __init__(self, slider_max: int, z_range_mm: float):
        """
        Args:
            slider_max: Maximum slider position (e.g., 200)
            z_range_mm: Total z-scan range in mm (e.g., 40.0)
        
        Maps slider 0..slider_max to centerpoint -z_range_mm/2..z_range_mm/2
        """
        # Centerpoint can be shifted ±half the z_range from origin
        center_range = z_range_mm / 2.0
        
        super().__init__(-center_range, center_range, slider_max)
        self.z_range_mm = z_range_mm
        self.center_range = center_range
    
    def slider_to_physical(self, slider_value: int) -> float:
        """Slider 0..slider_max maps to -z_range/2..z_range/2"""
        return super().slider_to_physical(slider_value)
    
    def physical_to_slider(self, physical_value: float) -> int:
        """Center value -z_range/2..z_range/2 maps to slider 0..slider_max"""
        return super().physical_to_slider(physical_value)