# Z-Scan v2 Bug Fix Implementation Guide

## Critical Issues (Fix First - Core Functionality)

### 1. **Slider Updates Not Reflected in UI** 🔴 HIGH
**File:** `zscan_main_v2.py` - `on_slider_changed()` method
- **Issue:** Moving sliders doesn't update fit curve or summary displays
- **Root Cause:** Slider value mapping not connected to physics calculation + no redraw
- **Fix:**
  ```python
  # Map slider value to physical quantity properly
  # Update fit summary spinboxes in real-time
  # Redraw plot after parameter change
  ```
- **Impact:** Users can't adjust fits interactively
- **Testing:** Move any slider and verify curve updates

---
<strike>
### 2. **ROI Reset When Unticking Fix ROI** 🔴 HIGH
**File:** `zscan_main_v2.py` - `_toggle_roi()` method
- **Issue:** After fitting in ROI, unticking "Fix ROI" doesn't clear selection
- **Root Cause:** `roi_limits` not cleared when checkbox unchecked
- **Fix:** 
  ```python
  def _toggle_roi_FIXED(self, canvas_key: str, enabled: bool):
      if not enabled:
          roi_selector.roi_limits = None  # CRITICAL
          # Redraw plot without ROI highlighting
  ```
- **Impact:** Users can't reset to full range without reloading data
- **Testing:** Load data → Set ROI → Uncheck Fix ROI → Verify full range shows
</strike>
---
<strike>
### 3. **Custom Checkbox Disables Value Updates** 🔴 HIGH
**File:** `zscan_main_v2.py` - `_on_custom_checkbox_changed()` method
- **Issue:** Unchecking "Set custom" greys out spinbox, prevents value updates from new files
- **Root Cause:** `setEnabled(False)` blocks value updates
- **Fix:**
  ```python
  def _on_custom_checkbox_changed_FIXED(self, spinbox_name, enabled: bool):
      widget.setReadOnly(not enabled)  # Prevent user edit, allow updates
      widget.setEnabled(True)           # Always enabled for display
  ```
- **Impact:** Loading new data doesn't update parameters if checkbox unchecked
- **Testing:** Uncheck "Set custom wavelength" → Load new file → Verify wavelength updates
</strike>
---

## High Priority Issues (Core Fitting Features)

### 4. **OA Fitting Doesn't Apply Model** 🟠 HIGH
**File:** `zscan_main_v2.py` - `on_fit_clicked()` for OA aperture
- **Issue:** OA Fit button clicked but curve doesn't change
- **Root Cause:** May need CA reference; possible missing param mappings
- **Fix:** Ensure CA fit exists first; validate beamwaist parameter
- **Testing:** Fit silica CA → Fit solvent OA → Verify curve updates

---

### 5. **Absorption Toggle Logic Reversed (OA)** 🟠 HIGH
**File:** `zscan_main_v2.py` - `_setup_absorption_initial_state()`
- **Issue:** "Assume no absorption" checkbox shows models when ON (should be opposite)
- **Root Cause:** Logic inverted in visibility check
- **Fix:**
  ```python
  # Checked = True should HIDE absorption controls
  # Unchecked = False should SHOW absorption controls
  ```
- **Testing:** Check "Assume no absorption" → Verify absorption dropdown disappears

---

### 6. **Fix ROI Toggle Doesn't Show Crosshairs** 🟠 MEDIUM
**File:** `zscan_main_v2.py` - `_toggle_roi()` method  
**OA Tab Specific**
- **Issue:** Checking "Fix ROI" on OA tab doesn't activate crosshairs
- **Root Cause:** ROI selector not initialized for OA tabs
- **Fix:** Initialize ROI selector for all CA tabs on load
- **Testing:** Load data → Go to OA tab → Check Fix ROI → Verify crosshairs appear

---

## Medium Priority Issues (Display & Calculation)

### 7. **Laser Intensity Missing/Wrong Scale** 🟡 MEDIUM
**File:** `zscan_main_v2.py` - `_update_result_display()`
- **Issue:** Laser intensity doesn't show or scale is wrong in fit summary
- **Root Cause:** Not calculated from beamwaist; calculation missing
- **Fix:**
  ```python
  intensity = power_W / (π * w0²)  # W/m²
  # Display in TW/m² scale
  ```
- **Testing:** Fit silica CA → Check laser intensity field updates

---

### 8. **n2 Not Calculated from DeltaPhi0** 🟡 MEDIUM
**File:** `zscan_main_v2.py` (Solvent/Sample CA tabs)
- **Issue:** Moving DeltaPhi0 slider doesn't update n2 in summary
- **Root Cause:** n2 calculation not connected to slider callback
- **Fix:** Add n2 calculation in `on_slider_changed()` when amplitude changes
- **Testing:** Load solvent CA → Move DeltaPhi0 slider → Verify n2 updates

---

### 9. **Z_Rayleigh Custom Toggle Not Enabled** 🟡 MEDIUM
**File:** `zscan_main_v2.py` (Solvent/Sample CA)
- **Issue:** "Custom Z_Rayleigh" toggle disabled on data load
- **Root Cause:** No initialization code for toggle state
- **Fix:**
  ```python
  def _setup_custom_rayleigh_toggles(self):
      # Enable toggle, set default unchecked
      # Connect slider enable/disable logic
  ```
- **Testing:** Load solvent/sample → Check custom Z_Rayleigh toggle is enabled

---

### 10. **Noise Reduction Sliders Don't Work** 🟡 MEDIUM
**File:** `zscan_main_v2.py` - `on_slider_changed()` (filter size parameter)
- **Issue:** Moving filter size slider doesn't reduce noise in data
- **Root Cause:** No filtering applied to data on slider change
- **Fix:** Apply Savitzky-Golay filter when filter_size slider changes
- **Testing:** Load data → Move filter size slider → Verify noise reduction in plot

---

## Low Priority Issues (UX & Menu)

### 11. **Menu Actions Not Connected** 🟢 LOW
**File:** `zscan_main_v2.py` - `__init__()` method
- **Issues:**
  - File → Load Solvents not connected
  - File → Exit doesn't confirm running fits
  - View → Theme light/dark not implemented
- **Fix:** Add `_connect_menu_actions()` method
- **Testing:** Click menu items and verify handlers execute

---

### 12. **Missing Solvent/Sample Load Warnings** 🟢 LOW
**File:** `zscan_main_v2.py` - `load_file()` method
- **Issue:** No warning if user loads solvent/sample before silica
- **Root Cause:** No validation check in load_file
- **Fix:**
  ```python
  if sample_type in ["solvent", "sample"]:
      if "silica_ca" not in self.data:
          QMessageBox.warning(...)
          return
  ```
- **Testing:** Try loading solvent before silica → Verify warning appears

---

## Implementation Checklist

### Phase 1: Critical Fixes (Do First)
- [ ] Fix slider value mapping and update callbacks
- [ ] Fix ROI reset logic when unchecking Fix ROI
- [ ] Fix custom checkbox to use setReadOnly instead of setEnabled
- [ ] Verify OA fitting returns valid results
- [ ] Fix absorption toggle logic reversal

### Phase 2: High Priority Features
- [ ] Add laser intensity calculation
- [ ] Connect n2 calculation to slider changes
- [ ] Initialize Z_Rayleigh toggle state
- [ ] Implement noise reduction filter application
- [ ] Fix ROI crosshairs for OA tabs

### Phase 3: Polish & UX
- [ ] Connect menu actions (Load Solvents, Exit, Theme)
- [ ] Add load validation warnings for solvent/sample
- [ ] Add error handling for edge cases
- [ ] Test all parameter combinations

---

## Testing Strategy

1. **Data Load Flow:**
   - Load silica CA/OA → Parameters populate ✓
   - Load solvent without silica → Warning appears ✓
   - Load sample without silica/solvent → Warning appears ✓

2. **Slider Interaction:**
   - Move any slider → Summary updates in real-time ✓
   - Move filter size → Plot shows filtered data ✓
   - Move DeltaPhi0 → n2 updates ✓

3. **ROI Selection:**
   - Check Fix ROI → Crosshairs appear ✓
   - Click twice → ROI selected ✓
   - Uncheck Fix ROI → ROI clears, full range shows ✓

4. **Fitting:**
   - Fit CA → Curve updates, R² shows ✓
   - Fit OA → Uses CA beamwaist, absorption model selected ✓
   - Change params → Invalidation warning shows ✓

5. **UI State:**
   - Check "Assume no absorption" → Absorption controls hide ✓
   - Uncheck "Assume no absorption" → Absorption controls show ✓
   - Check "Custom Z_Rayleigh" → Slider enabled ✓

---

## Priority Summary

**Do Today:**
1. Slider updates (affects every fit)
2. ROI reset (breaks workflow)
3. Custom checkbox readonly (blocks data loading)

**Do This Week:**
4. OA fitting validation
5. Absorption logic fix
6. Laser intensity + n2 calcs

**Nice to Have:**
7. Menu connections
8. Load warnings
9. Filter implementation