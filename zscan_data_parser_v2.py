"""
Clean Z-scan data parsing - tested with silica and mBP files
No dependencies on UI or physics - pure data processing
Added: Track missing header fields for corrupted files
"""

from dataclasses import dataclass, field
from pathlib import Path
import numpy as np
from typing import Optional, List


@dataclass
class RawZScanData:
    """Raw data from file"""
    position_mm: np.ndarray
    ca_raw: np.ndarray  # Channel 1: Closed aperture
    ref: np.ndarray      # Channel 2: Reference
    oa_raw: np.ndarray   # Channel 3: Open aperture
    start_pos: float
    end_pos: float
    wavelength_nm: float
    sample_code: str
    missing_fields: List[str] = field(default_factory=list)  # Track missing header fields


@dataclass
class ProcessedZScanData:
    """Normalized and processed data"""
    position_mm: np.ndarray
    position_centered: np.ndarray  # Centered around focal point
    ca_raw: np.ndarray
    ca: np.ndarray
    oa: np.ndarray
    ref: np.ndarray
    wavelength_nm: float
    sample_code: str
    z_range_mm: float


class ZScanFileParser:
    """Parse Z-scan measurement files"""
    
    @staticmethod
    def parse(filepath: str) -> Optional[RawZScanData]:
        """
        Parse Z-scan data file.
        Expected format:
        - Header with metadata
        - "SNo." line marks start of data
        - 4 columns: CA, Reference, OA, Empty
        
        Returns RawZScanData or None if parsing fails.
        Tracks missing header fields for corrupted files.
        """
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            # Initialize metadata with defaults - use separate typed variables
            wavelength: float = 800.0
            start_pos: float = 40.0
            end_pos: float = 80.0
            code: str = Path(filepath).stem
            missing_fields: List[str] = []
            
            found_wavelength = False
            found_start_pos = False
            found_end_pos = False
            found_code = False
            
            # Extract metadata from header
            for line in lines[:30]:
                if 'Wavelength:' in line:
                    try:
                        header_parts = line.split(':')  # RENAMED from 'parts'
                        if len(header_parts) > 1:
                            wl_str = header_parts[1].strip().split()[0]
                            wavelength = float(wl_str)
                            found_wavelength = True
                    except (ValueError, IndexError):
                        pass
                
                if 'Code:' in line or 'Sample:' in line:
                    try:
                        code_str = line.split(':')[1].strip()
                        if code_str:
                            code = code_str
                            found_code = True
                    except IndexError:
                        pass
                
                if 'Starting pos:' in line or 'Start pos:' in line:
                    try:
                        start_pos = float(line.split(':')[1].strip())
                        found_start_pos = True
                    except (ValueError, IndexError):
                        pass
                
                if 'Ending pos:' in line or 'End pos:' in line:
                    try:
                        end_pos = float(line.split(':')[1].strip())
                        found_end_pos = True
                    except (ValueError, IndexError):
                        pass
            
            # Track missing fields
            if not found_wavelength:
                missing_fields.append('wavelength')
            if not found_start_pos:
                missing_fields.append('start_position')
            if not found_end_pos:
                missing_fields.append('end_position')
            if not found_code:
                missing_fields.append('sample_code')
            
            # Find data start
            data_start = 0
            for i, line in enumerate(lines):
                if 'SNo.' in line:
                    data_start = i + 2
                    break
            
            # Parse data rows
            data_rows: List[List[float]] = []
            for line in lines[data_start:]:
                line = line.strip()
                if not line or line.startswith('-'):
                    continue
                
                try:
                    # RENAMED from 'parts' to 'row_values' to avoid type confusion
                    row_values: List[float] = [float(x) for x in line.split() if x]
                    if len(row_values) >= 4:
                        # Skip the serial number, take the 4 voltage columns
                        data_rows.append(row_values[1:5])
                except (ValueError, IndexError):
                    continue
            
            if not data_rows:
                print(f"No data rows found in {filepath}")
                return None
            
            data_array = np.array(data_rows)
            n_points = len(data_array)
            
            # Generate position array (linear interpolation between start and end)
            position_mm = np.linspace(start_pos, end_pos, n_points)
            
            return RawZScanData(
                position_mm=position_mm,
                ca_raw=data_array[:, 0],
                ref=data_array[:, 1],
                oa_raw=data_array[:, 2],
                start_pos=start_pos,
                end_pos=end_pos,
                wavelength_nm=wavelength,
                sample_code=code,
                missing_fields=missing_fields,
            )
        
        except Exception as e:
            print(f"Error parsing {filepath}: {e}")
            return None


class ZScanProcessor:
    """Process raw Z-scan data into usable format"""
    
    @staticmethod
    def normalize(raw: RawZScanData) -> ProcessedZScanData:
        """
        Normalize data:
        0. Ensure all data is in ascending position order
        1. Divide CA and OA by reference
        2. Shift to zero level = 1.0
        3. Center positions around focal point
        """
        z_range = abs(raw.end_pos - raw.start_pos)
        
        # ================================================================
        # CHECK IF DATA IS REVERSED (backward scan)
        # ================================================================
        is_descending = raw.start_pos > raw.end_pos

        if is_descending:
            print(f"[normalize] Backward scan detected: {raw.start_pos} → {raw.end_pos}")
            print("[normalize] Reversing all arrays to ascending order")
            
            # Reverse everything to make it ascending
            ca_raw = raw.ca_raw[::-1]
            oa_raw = raw.oa_raw[::-1]
            ref = raw.ref[::-1]
            position_mm = raw.position_mm[::-1]
            
            # Update metadata
            actual_start = raw.end_pos
            actual_end = raw.start_pos
        else:
            print(f"[normalize] Forward scan: {raw.start_pos} → {raw.end_pos}")
            
            ca_raw = raw.ca_raw
            oa_raw = raw.oa_raw
            ref = raw.ref
            position_mm = raw.position_mm
            actual_start = raw.start_pos
            actual_end = raw.end_pos
        
        # ================================================================
        # Normalize by reference
        # ================================================================
        
        ca = ca_raw / ref
        oa = oa_raw / ref
        
        # Shift to zero level 1.0
        ca = 1.0 + (ca - np.mean(ca))
        oa = 1.0 + (oa - np.mean(oa))
        
        # ================================================================
        # Center positions
        # ================================================================
        center_pos = (actual_start + actual_end) / 2.0
        position_centered = position_mm - center_pos
        
        # Verify ascending order
        assert position_centered[0] < position_centered[-1], \
            f"ERROR: Data not ascending! {position_centered[0]} → {position_centered[-1]}"
        print(f"[normalize] Data verified ascending: {position_centered[0]:.1f} → {position_centered[-1]:.1f} mm")

        return ProcessedZScanData(
            position_mm=position_mm,
            position_centered=position_centered,
            ca_raw=ca_raw,
            ca=ca,
            oa=oa,
            ref=ref,
            wavelength_nm=raw.wavelength_nm,
            sample_code=raw.sample_code,
            z_range_mm=z_range,
        )


if __name__ == '__main__':
    # Test with your data files
    test_files = [
        'data/2021_07_22__10_10__silica_0-0_1600-0_2.txt',
        'data/2021_07_22__10_47__mBP_0-40_1600-0_2.txt',
    ]
    
    for filepath in test_files:
        if Path(filepath).exists():
            print(f"\n{'='*60}")
            print(f"Testing: {filepath}")
            print('='*60)
            
            raw = ZScanFileParser.parse(filepath)
            if raw:
                print(f"✓ Parsed: {raw.sample_code}")
                print(f"  Points: {len(raw.ca_raw)}")
                print(f"  Range: {raw.start_pos} - {raw.end_pos} mm")
                print(f"  λ: {raw.wavelength_nm} nm")
                
                if raw.missing_fields:
                    print(f"  ⚠ Missing fields: {', '.join(raw.missing_fields)}")
                
                processed = ZScanProcessor.normalize(raw)
                print("✓ Processed:")
                print(f"  CA mean: {np.mean(processed.ca):.4f} (should be ~1.0)")
                print(f"  OA mean: {np.mean(processed.oa):.4f} (should be ~1.0)")
                print(f"  CA range: [{np.min(processed.ca):.4f}, {np.max(processed.ca):.4f}]")
                print(f"  Position range: [{np.min(processed.position_centered):.2f}, {np.max(processed.position_centered):.2f}] mm")
                
                # Check 2: Position should be centered
                assert abs(float(np.mean(processed.position_centered))) < 0.1, "Position should be centered"
                print(f"  ✓ Position centering correct (mean={np.mean(processed.position_centered):.6f})")
                
                # Check 3: Data variance should be reasonable
                ca_std = float(np.std(processed.ca))
                assert 0.01 < ca_std < 0.2, "CA std_dev should be reasonable"
                print(f"  ✓ CA std_dev is reasonable ({ca_std:.4f})")
                
                print("\n✓ ALL CHECKS PASSED - Data is valid!")
            else:
                print("✗ Failed to parse")
        else:
            print(f"✗ File not found: {filepath}")