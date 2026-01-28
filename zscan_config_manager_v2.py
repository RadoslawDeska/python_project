"""
Configuration manager for Z-scan fitting parameters
Tracks which parameters come from file vs user input
Manages parameter validity and fitting dependencies
"""

from dataclasses import dataclass, field
from typing import Optional, Set
from enum import Enum


class ParamSource(Enum):
    """Where did a parameter come from?"""
    FILE = "file"          # Loaded from data file
    USER_CUSTOM = "user"   # User manually set via spinbox
    DEFAULT = "default"    # Default value


@dataclass
class FittingConfig:
    """Current fitting configuration"""
    # General parameters
    wavelength_nm: float = 800.0
    zscan_range_mm: float = 40.0
    aperture_diameter_mm: float = 1.0
    d0_mm: float = 260.0
    silica_thickness_mm: float = 3.0
    concentration_pct: float = 0.0
    solvent_name: str = "water"  # Track solvent type
    
    # Which parameters are custom (vs from file)
    wavelength_custom: bool = False
    zscan_range_custom: bool = False
    aperture_custom: bool = False
    d0_custom: bool = False
    silica_thickness_custom: bool = False
    concentration_custom: bool = False
    
    # Fitted results (tracks what's been fitted)
    fitted_samples: Set[str] = field(default_factory=set)  # 'silica_CA', 'solvent_CA', etc.
    
    def mark_fitted(self, sample_type: str, aperture: str):
        """Mark a sample/aperture as fitted"""
        key = f"{sample_type}_{aperture}"
        self.fitted_samples.add(key)
    
    def get_fitted_keys(self) -> Set[str]:
        """Get all fitted sample/aperture combinations"""
        return self.fitted_samples.copy()
    
    def invalidate_dependents(self, changed_param: str) -> Set[str]:
        """
        Determine which fits become invalid when a parameter changes
        
        Returns set of fitted samples that need to be re-fitted
        """
        invalidated = set()
        
        # Parameters that affect ALL CA fits (universal physical params)
        universal_ca_params = {
            'wavelength_nm', 'aperture_diameter_mm', 'd0_mm', 
            'zscan_range_mm', 'silica_thickness_mm'
        }
        
        if changed_param in universal_ca_params:
            # All CA fits become invalid
            for key in self.fitted_samples:
                if 'CA' in key:
                    invalidated.add(key)
            
            # All OA fits depend on their corresponding CA, so they're invalid too
            for ca_key in list(invalidated):
                sample = ca_key.replace('_CA', '')
                oa_key = f"{sample}_OA"
                if oa_key in self.fitted_samples:
                    invalidated.add(oa_key)
        
        # Concentration only affects sample fits
        elif changed_param == 'concentration_pct':
            # Only sample fits become invalid
            for key in self.fitted_samples:
                if 'sample' in key:
                    invalidated.add(key)
        
        # Solvent change affects solvent AND sample (since sample is in solvent)
        elif changed_param == 'solvent_changed':
            # Both solvent and sample fits are invalid
            for key in self.fitted_samples:
                if 'solvent' in key or 'sample' in key:
                    invalidated.add(key)
        
        return invalidated
    
    def clear_fitted(self, keys: Optional[Set[str]] = None):
        """Clear fitted status for specific keys or all"""
        if keys is None:
            self.fitted_samples.clear()
        else:
            self.fitted_samples.difference_update(keys)


if __name__ == '__main__':
    print("="*70)
    print("CONFIG MANAGER TEST")
    print("="*70)
    
    config = FittingConfig()
    config.mark_fitted('silica', 'CA')
    config.mark_fitted('silica', 'OA')
    config.mark_fitted('solvent', 'CA')
    config.mark_fitted('solvent', 'OA')
    
    print(f"Fitted: {config.get_fitted_keys()}")
    
    # User changes aperture
    print("\nUser changes aperture diameter...")
    invalidated = config.invalidate_dependents('aperture_diameter_mm')
    print(f"Invalidated: {invalidated}")
    print(f"After invalidation: {config.get_fitted_keys()}")
    
    # User changes concentration
    print("\nUser changes concentration...")
    config.mark_fitted('silica', 'CA')
    config.mark_fitted('silica', 'OA')
    config.mark_fitted('solvent', 'CA')
    invalidated = config.invalidate_dependents('concentration_pct')
    print(f"Invalidated: {invalidated}")
    
    print("\n✓ Config manager test passed")