# Missing/faulty connections

## Menu actions
1. File->Load solvents is not connected to loading solvents file
2. File->Exit should call exiting the program with confirmation about any pending measurements or fitting process.
3. View->Theme should have working connections to switch between light and dark mode.

## Data fitting
### General problems and common to all fitting tabs
1. Data directory **Custom dir**
2. <strike>Unticking **Set custom** results in disabling (greying-out) of respective general parameters and sample properties -> This further results in values not updating when loading new data</strike>
3. Loading from file doesn't result in first manual fit (based on data geometry) neither for CA nor OA data. The line shows up only after **Fit** button.
4. Noise reduction sliders don't result in noise reduction in data or the data is not updated on the canvas.
5. Moving sliders doesn't modify the fitting curve on neither of the canvases -> Check whether it is the problem of scaling or connections.
6. <strike>After fitting the data in a ROI: if *Fix ROI* is unticked and **Fit** button is clicked, the data doesn't fit over the whole data range (ROI doesn't reset to full range).</strike>


### Silica fitting tabs
1. Moving *Z_r* slider doesn't result in changing the value of beamwaist nor the numerical aperture in the Fit summary
2. Laser intensity doesn't show up at all or check the scale of the value to be displayed on the Fit summary.
3. Moving *DeltaPhi0* slider doesn't update Laser intensity in the Fit summary.


### Solvent fitting tabs
#### General
1. After loading solvent file (prior to loading silica file) there is no warning for the user to fit silica first.

#### CA fitting tab
1. *DeltaPhi0* slider has different range than the display range in the fitting summary.
2. *n_2* doesn't get calculated or the scale is wrong so no effect is visible in the Fit summary. It should update when moving *DeltaPhi0* slider.
3. Custom *Z_Rayleigh* toggle is not enabled on data load. Make sure that when toggled to **ON** state results in enabling the *Z_Rayleigh* slider and **OFF** state disables the slider (default state on data load).

#### OA fitting tab
1. Moving sliders doesn't update the fit summary.
2. *Assume no absorption* toggle is reversed, models show up when it is in the **ON** state (should be the opposite).
3. Clicking **Fit** doesn't result in any model application to fit - no change in the curve.
    - File: zscan_main_v2.py - on_fit_clicked() for OA aperture
    - Issue: OA Fit button clicked but curve doesn't change
    - Root Cause: May need CA reference; possible missing param mappings
    - Fix: Ensure CA fit exists first; validate beamwaist parameter
    - Testing: Fit silica CA → Fit solvent OA → Verify curve updates

4. *Fix ROI* toggle doesn't activate crosshairs on the canvas.
5. Custom *Center* toggle is not required and should not block fitting of the center parameter.


### Sample fitting tabs
#### General
1. After loading sample file (prior to loading silica and/or solvent file) there is no warning for the user to fit silica/solvent first.

#### CA fitting tab
1. Moving sliders doesn't update the fit summary.
2. *n_2* doesn't get calculated or the scale is wrong so no effect is visible in the Fit summary. It should update when moving *DeltaPhi0* slider.
3. Custom *Z_Rayleigh* toggle is not enabled on data load. Make sure that when toggled to **ON** state results in enabling the *Z_Rayleigh* slider and **OFF** state disables the slider (default state on data load).

#### OA fitting tab
1. Moving sliders doesn't update the fit summary.
2. *n_2* should read the value from CA fit.
3. *T* should update with moving *DeltaT* slider. *beta* and/or *gamma* should update as well (depending on the model chosen; the saturation intensity should appear in the fit summary when saturation model is selected)
