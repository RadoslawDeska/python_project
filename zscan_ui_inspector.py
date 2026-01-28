"""
Inspect window.ui structure by parsing XML directly
"""

import xml.etree.ElementTree as ET
from pathlib import Path
from collections import defaultdict


def inspect_ui_xml(ui_path: str):
    """Parse UI file XML and extract widget names and types"""
    
    tree = ET.parse(ui_path)
    root = tree.getroot()
    
    # Extract namespace if present
    ns = {'ui': 'http://www.trolltech.com/Qt'}
    if '}' in root.tag:
        ns = {'ui': root.tag.split('}')[0][1:]}
    
    print("=" * 100)
    print("UI WIDGET STRUCTURE - PARSED FROM XML")
    print("=" * 100)
    print()
    
    # Find all widgets
    widgets_by_type = defaultdict(list)
    widgets_by_name = {}
    
    for widget in root.findall('.//widget', ns):
        name = widget.get('name', 'NO_NAME')
        wtype = widget.get('class', 'UNKNOWN')
        
        widgets_by_type[wtype].append(name)
        widgets_by_name[name] = wtype
        
    # Print by type
    print("WIDGETS BY TYPE:")
    print("-" * 100)
    for wtype in sorted(widgets_by_type.keys()):
        names = sorted(widgets_by_type[wtype])
        print(f"\n{wtype}:")
        for name in names:
            print(f"  - {name}")
    
    print()
    print("=" * 100)
    print("RELEVANT CONTROLS FOR Z-SCAN FITTING:")
    print("=" * 100)
    print()
    
    # Look for CA/OA specific widgets
    print("SLIDERS (sorted by sample type and aperture):")
    print("-" * 100)
    sliders = {name: wtype for name, wtype in widgets_by_name.items() if 'Slider' in wtype}
    for slider_name in sorted(sliders.keys()):
        print(f"  {slider_name}")
    
    print()
    print("FIT BUTTONS (sorted by sample type and aperture):")
    print("-" * 100)
    buttons = {name: wtype for name, wtype in widgets_by_name.items() if 'PushButton' in wtype and 'fit' in name.lower()}
    for button_name in sorted(buttons.keys()):
        print(f"  {button_name}")
    
    print()
    print("SPINBOXES (sorted):")
    print("-" * 100)
    spinboxes = {name: wtype for name, wtype in widgets_by_name.items() if 'SpinBox' in wtype}
    for spinbox_name in sorted(spinboxes.keys()):
        print(f"  {spinbox_name}")
    
    print()
    print("=" * 100)
    print("CA/OA SLIDER MAPPING (Ready for connection):")
    print("=" * 100)
    print()
    
    # Extract CA/OA sliders pattern
    ca_oa_sliders = {}
    for name in sliders.keys():
        if 'CA' in name.upper() or 'OA' in name.upper():
            sample_type = None
            aperture = None
            param = None
            
            if 'silica' in name.lower():
                sample_type = 'silica'
            elif 'solvent' in name.lower():
                sample_type = 'solvent'
            elif 'sample' in name.lower():
                sample_type = 'sample'
            
            if 'CA' in name.upper():
                aperture = 'CA'
            elif 'OA' in name.upper():
                aperture = 'OA'
            
            # Extract parameter
            if 'zeroLevel' in name or 'zero' in name.lower():
                param = 'zero_level'
            elif 'DPhi0' in name or 'deltaphi' in name.lower():
                param = 'amplitude_dphi0'
            elif 'centerPoint' in name or 'center' in name.lower():
                param = 'centerpoint'
            elif 'Rayleigh' in name or 'rayleigh' in name.lower():
                param = 'rayleigh'
            elif 'filterSize' in name or 'filter' in name.lower():
                param = 'filter_size'
            else:
                param = 'other'
            
            key = f"{sample_type}_{aperture}_{param}" if sample_type and aperture else name
            ca_oa_sliders[key] = name
    
    for key, name in sorted(ca_oa_sliders.items()):
        print(f"  {key:<40} -> {name}")
    
    print()
    print("=" * 100)


if __name__ == '__main__':
    ui_path = Path(__file__).parent / "window.ui"
    
    if not ui_path.exists():
        print(f"ERROR: UI file not found at {ui_path}")
        import sys
        sys.exit(1)
    
    inspect_ui_xml(str(ui_path))