"""
Generate reversed Z-scan data file for testing backward-scan handling

This creates an artificial silica file with:
- Same data as original silica file
- Rows reversed (backward scan order)
- Row numbering intact (1-201)
- Position range reversed (80-40 instead of 40-80)
"""

import re
from pathlib import Path
from typing import List


def reverse_zscan_file(input_file, output_file, keep_row_numbers=True):
    """
    Reverse the order of data rows in a Z-scan file

    Args:
        input_file: Path to original Z-scan file
        output_file: Path to output reversed file
        keep_row_numbers: If True, keep row numbering 1-N (with reversed data)
                         If False, reverse both data and row numbers

    This preserves the header but reverses the data rows, which simulates
    a backward scan of the same sample.
    """

    print(f"Opening: {input_file}")
    with open(input_file, "r") as f:
        lines = f.readlines()

    # Find header/data boundary
    data_start = 0
    for i, line in enumerate(lines):
        if "SNo." in line:
            data_start = i + 4  # Skip "SNo." line and the dashes
            break

    print(f"Data starts at line {data_start}")

    # Split into header and data
    header = lines[:data_start]
    data_lines = lines[data_start:]

    def swap_positions_in_lines(lines: List[str]) -> List[str]:
        start_idx = None
        end_idx = None
        start_val = None
        end_val = None

        # znajdź linie i wartości
        for i, line in enumerate(lines):
            m_start = re.match(r"(Starting pos:\s*)(\d+)", line)
            m_end = re.match(r"(Ending pos:\s*)(\d+)", line)

            if m_start:
                start_idx = i
                start_val = m_start.group(2)

            if m_end:
                end_idx = i
                end_val = m_end.group(2)

        # jeśli nie znaleziono obu → nic nie zmieniamy
        if start_idx is None or end_idx is None or start_val is None or end_val is None:
            return lines

        # zamiana wartości z użyciem funkcji jako replacera
        def replace_start(match: re.Match[str]) -> str:
            prefix = match.group(1)
            return f"{prefix}{end_val}"

        def replace_end(match: re.Match[str]) -> str:
            prefix = match.group(1)
            return f"{prefix}{start_val}"

        lines[start_idx] = re.sub(r"(Starting pos:\s*)(\d+)", replace_start, lines[start_idx])
        lines[end_idx] = re.sub(r"(Ending pos:\s*)(\d+)", replace_end, lines[end_idx])

        return lines

    header = swap_positions_in_lines(header)

    # Parse data rows
    data_rows = []
    for line in data_lines:
        line = line.strip()
        if not line or line.startswith("-"):
            continue

        try:
            # Parse: SNo | CH1 | CH2 | CH3 | CH4
            values = [float(x) for x in line.split() if x]
            if len(values) >= 5:
                row_num = int(values[0])
                ch1, ch2, ch3, ch4 = values[1:5]
                data_rows.append(
                    {
                        "num": row_num,
                        "ch1": ch1,
                        "ch2": ch2,
                        "ch3": ch3,
                        "ch4": ch4,
                    }
                )
        except (ValueError, IndexError):
            continue

    print(f"Parsed {len(data_rows)} data rows")

    # Reverse data
    reversed_data = data_rows[::-1]

    # Re-number if requested
    if keep_row_numbers:
        # Keep original row numbers but with reversed data
        for i, row in enumerate(reversed_data):
            row["num"] = i + 1

    # Reconstruct file
    output_lines = []

    # Add header
    output_lines.extend(header)

    # Add reversed data
    for row in reversed_data:
        line = f"{int(row['num']):5d}     {row['ch1']:.8f}             {row['ch2']:.8f}             {row['ch3']:.8f}             {row['ch4']:.8f}\n"
        output_lines.append(line)

    # Write output
    print(f"Writing: {output_file}")
    with open(output_file, "w") as f:
        f.writelines(output_lines)

    print(f"✓ Reversed file created!")
    print(f"  Original rows: 1-{len(data_rows)}")
    print(
        f"  Reversed rows: {len(reversed_data)}-1 (with renumbered as 1-{len(reversed_data)})"
    )

    return len(reversed_data)


def create_reversed_silica_file():
    """
    Create reversed silica test file for testing backward scans

    Usage:
        python zscan_data_generator.py
    """

    input_path = Path("data/2021_07_22__10_16__dcm_0-0_1600-0_2.txt")
    output_path = Path(
        "data/2021_07_22__10_16__dcm_0-0_1600-0_2_REVERSED.txt"
    )

    if not input_path.exists():
        print(f"ERROR: Input file not found: {input_path}")
        return False

    print("=" * 80)
    print("Z-SCAN DATA REVERSAL TOOL")
    print("=" * 80)
    print()
    print(f"Creating backward-scan test file...")
    print(f"Input:  {input_path}")
    print(f"Output: {output_path}")
    print()

    num_rows = reverse_zscan_file(
        input_file=str(input_path),
        output_file=str(output_path),
        keep_row_numbers=True,
    )

    print()
    print("=" * 80)
    print("TESTING NOTES")
    print("=" * 80)
    print(f"""
The reversed file simulates a backward scan of the same sample.

Original scan (forward, 40→80):
  - Position increases from 40 to 80 mm
  - Row 1 at position 40, Row {num_rows} at position 80
  
Reversed scan (backward, 80→40):
  - Position decreases from 80 to 40 mm
  - Row 1 at position 80, Row {num_rows} at position 40
  - Same data values, but in reverse order
  
The Z-scan physics:
  - Z-position is reversed relative to focal point
  - Peak/valley positions should appear opposite
  - DPhi0 should have OPPOSITE sign if algorithm doesn't account for scan direction

HOW TO TEST:
1. Load: {output_path} as solvent data
2. Compare fit results with original silica:
   - If curve shape is INVERTED → Backward scan detected, need sign correction
   - If curve shape is SAME → Algorithm already handles it
   
3. Expected behavior (if algorithm is correct):
   - Forward silica: DPhi0 ≈ +0.48 rad
   - Backward solvent: DPhi0 ≈ -0.48 rad (opposite sign)
     OR
   - Forward silica: DPhi0 ≈ +0.48 rad  
   - Backward solvent: DPhi0 ≈ +0.48 rad (same sign, auto-corrected)
   
   Your current result: DPhi0 ≈ -0.43 rad
   This matches the backward scan expectation (opposite sign)!
""")

    return True


def analyze_both_files():
    """
    Analyze both original and reversed files to show the difference
    """

    original = Path("data/2021_07_22__10_10__silica_0-0_1600-0_2.txt")
    reversed_file = Path(
        "data/2021_07_22__10_10__silica_0-0_1600-0_2_REVERSED.txt"
    )

    if not original.exists() or not reversed_file.exists():
        print("ERROR: Files not found")
        return

    def get_data(filepath):
        data = []
        with open(filepath, "r") as f:
            for line in f:
                if "SNo." in line:
                    continue
                line = line.strip()
                if not line or line.startswith("-"):
                    continue
                try:
                    values = [float(x) for x in line.split() if x]
                    if len(values) >= 5:
                        data.append(values[1])  # CH1 (CA)
                except:
                    pass
        return data

    orig_data = get_data(original)
    rev_data = get_data(reversed_file)

    print("\n" + "=" * 80)
    print("DATA COMPARISON")
    print("=" * 80)

    print(f"\nOriginal file (forward scan):")
    print(f"  Rows: {len(orig_data)}")
    print(f"  Row 1 CA value: {orig_data[0]:.8f}")
    print(f"  Row {len(orig_data)} CA value: {orig_data[-1]:.8f}")

    print(f"\nReversed file (backward scan):")
    print(f"  Rows: {len(rev_data)}")
    print(f"  Row 1 CA value: {rev_data[0]:.8f}")
    print(f"  Row {len(rev_data)} CA value: {rev_data[-1]:.8f}")

    print(f"\nData order check:")
    print(
        f"  Original row 1 = Reversed row {len(rev_data)}: {abs(orig_data[0] - rev_data[-1]) < 1e-6}"
    )
    print(
        f"  Original row {len(orig_data)} = Reversed row 1: {abs(orig_data[-1] - rev_data[0]) < 1e-6}"
    )

    # Check Z-position mapping
    print(f"\nZ-position mapping:")
    print(f"  Original: position 40→80 (row 1 to {len(orig_data)})")
    print(f"  Reversed: position 80→40 (row 1 to {len(rev_data)})")
    print(f"  This reversal affects peak/valley interpretation!")


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "analyze":
        analyze_both_files()
    else:
        success = create_reversed_silica_file()
        if success:
            analyze_both_files()
