# tab_manager.py
from typing import Dict, List, Optional
from PyQt5.QtWidgets import QTabWidget, QWidget


class TabManager:
    """Manages nested tab navigation for material type and aperture selection"""
    
    def __init__(self, main_tab_widget: QTabWidget):
        """
        Initialize TabManager with the main tab widget
        
        Args:
            main_tab_widget: The main QTabWidget containing material tabs (silica, solvent, sample)
        """
        self.main_tabs: QTabWidget = main_tab_widget
        self.tab_mapping: Dict[str, Dict] = self._build_mapping()
    
    def _build_mapping(self) -> Dict[str, Dict]:
        """Build nested mapping of material -> apertures -> indices"""
        mapping: Dict[str, Dict] = {}
        
        for i in range(self.main_tabs.count()):
            material = self.main_tabs.tabText(i).lower()
            mapping[material] = {"main_index": i, "apertures": {}}
            
            # Find the nested QTabWidget for apertures
            main_tab: Optional[QWidget] = self.main_tabs.widget(i)
            if main_tab is None:
                continue
            
            aperture_tabs: Optional[QTabWidget] = main_tab.findChild(QTabWidget)
            
            if aperture_tabs is not None:
                for j in range(aperture_tabs.count()):
                    aperture_text = aperture_tabs.tabText(j).lower()
                    mapping[material]["apertures"][aperture_text] = j
                    
                    # Add shorthand keys (ca/oa) for easy access
                    if "closed" in aperture_text:
                        mapping[material]["apertures"]["ca"] = j
                    elif "open" in aperture_text:
                        mapping[material]["apertures"]["oa"] = j
        
        return mapping
    
    def switch_to(self, material: str, aperture: str = "CA") -> bool:
        """
        Switch to specified material and aperture tabs
        
        Args:
            material: Material type ("silica", "solvent", or "sample")
            aperture: Aperture type ("CA", "OA", or full text name) - default: CA
        
        Returns:
            bool: True if successful, False if tabs not found
        
        Example:
            manager.switch_to("silica")           # silica, Closed Aperture
            manager.switch_to("solvent", "OA")    # solvent, Open Aperture
            manager.switch_to("sample", "closed aperture")
        """
        material = material.lower()
        
        # Validate material exists
        if material not in self.tab_mapping:
            print(f"Material '{material}' not found. Available: {list(self.tab_mapping.keys())}")
            return False
        
        # Switch main tab
        main_index: int = self.tab_mapping[material]["main_index"]
        self.main_tabs.setCurrentIndex(main_index)
        
        # Switch aperture tab
        aperture_lower = aperture.lower()
        aperture_index: Optional[int] = self.tab_mapping[material]["apertures"].get(aperture_lower)
        
        if aperture_index is not None:
            main_tab: Optional[QWidget] = self.main_tabs.widget(main_index)
            if main_tab is None:
                return False
            
            aperture_tabs: Optional[QTabWidget] = main_tab.findChild(QTabWidget)
            if aperture_tabs is not None:
                aperture_tabs.setCurrentIndex(aperture_index)
                return True
        
        print(f"Aperture '{aperture}' not found for '{material}'")
        return False
    
    def get_current_material(self) -> str:
        """Get the currently active material tab name"""
        return self.main_tabs.tabText(self.main_tabs.currentIndex()).lower()
    
    def get_current_aperture(self) -> Optional[str]:
        """Get the currently active aperture tab name"""
        main_tab: Optional[QWidget] = self.main_tabs.currentWidget()
        if main_tab is None:
            return None
        
        aperture_tabs: Optional[QTabWidget] = main_tab.findChild(QTabWidget)
        
        if aperture_tabs is not None:
            return aperture_tabs.tabText(aperture_tabs.currentIndex()).lower()
        
        return None
    
    def get_available_materials(self) -> List[str]:
        """Get list of available material types"""
        return list(self.tab_mapping.keys())
    
    def get_available_apertures(self, material: str) -> List[str]:
        """Get list of available apertures for a material"""
        material = material.lower()
        if material in self.tab_mapping:
            # Filter out shorthand keys (ca, oa)
            apertures = [k for k in self.tab_mapping[material]["apertures"].keys() if len(k) > 2]
            return apertures
        return []
    
    def refresh(self) -> None:
        """Rebuild mapping (use if tabs are dynamically added/removed)"""
        self.tab_mapping = self._build_mapping()