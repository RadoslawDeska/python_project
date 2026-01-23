# tests/test_slider_mapper.py
import pytest
from lib.slider_mapper import ParameterSliderMapper, CenterPointSliderMapper


class TestParameterSliderMapper:
    
    def test_slider_to_physical_at_min(self):
        mapper = ParameterSliderMapper(-2, 2, 100)
        assert mapper.slider_to_physical(0) == -2
    
    def test_slider_to_physical_at_max(self):
        mapper = ParameterSliderMapper(-2, 2, 100)
        assert mapper.slider_to_physical(100) == 2
    
    def test_slider_to_physical_at_midpoint(self):
        mapper = ParameterSliderMapper(-2, 2, 100)
        assert mapper.slider_to_physical(50) == 0
    
    def test_physical_to_slider_at_min(self):
        mapper = ParameterSliderMapper(-2, 2, 100)
        assert mapper.physical_to_slider(-2) == 0
    
    def test_physical_to_slider_at_max(self):
        mapper = ParameterSliderMapper(-2, 2, 100)
        assert mapper.physical_to_slider(2) == 100
    
    def test_physical_to_slider_roundtrip(self):
        mapper = ParameterSliderMapper(-2, 2, 100)
        original = 0.75
        slider = mapper.physical_to_slider(original)
        recovered = mapper.slider_to_physical(slider)
        assert abs(recovered - original) < 0.01
    
    def test_invalid_slider_value_raises(self):
        mapper = ParameterSliderMapper(-2, 2, 100)
        with pytest.raises(ValueError):
            mapper.slider_to_physical(101)
    
    def test_invalid_physical_value_raises(self):
        mapper = ParameterSliderMapper(-2, 2, 100)
        with pytest.raises(ValueError):
            mapper.physical_to_slider(5)


class TestCenterPointSliderMapper:
    
    def test_slider_0_maps_to_minus_50(self):
        mapper = CenterPointSliderMapper(100, center_range=50)
        assert mapper.slider_to_physical(0) == -50
    
    def test_slider_100_maps_to_plus_50(self):
        mapper = CenterPointSliderMapper(100, center_range=50)
        assert mapper.slider_to_physical(100) == 50
    
    def test_slider_50_maps_to_zero(self):
        mapper = CenterPointSliderMapper(100, center_range=50)
        assert mapper.slider_to_physical(50) == 0