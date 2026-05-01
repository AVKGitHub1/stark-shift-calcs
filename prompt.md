# General Stark Shift Calculator

@calcs.py  @pyqt_gui.py @marimo_gui.py

Using calculations from the arc library, create a gui with:

## Inputs

Atom Species: Drop down with all the alkali atom species in arc
n, l, J, F, mF: atomic level quantum numbers, each should have its own input box. Allow for mF to be empty
Wavelength: wavelength of light this will be an input box with the 'nm' unit specified on the side of the box
Polarization: have three input boxes each for how much part of the power is in sigma-, pi, sigma+
Power: have an input text box and a unit selection drop down box nect to it with {uW mW, W}
Waist: have an input text box and a unit selection drop down box nect to it with {um, mm}
Distance of atom from waist: have an input text box and a unit selection drop down box nect to it with {nm, um, mm, m}

## Calculation

- Verify if the given combination of quantum numbers is a correct state, if not tell the user so. If the mF state is not specified then the calculations will be done for all possible mF states given the F.
- Normalize the polarization power split
- Using ARC determine what the scalar, vector, and tensor stark shift are.

## Outputs

- If the mF state is specified just output as text the 3 shifts.
- If the mF state is not specified then draw on plot a level diagram with all the mF states for the given F and indicate with arrows the direction and size of each of the three types of shifts on all the mF states. Also have the value next to the arrows. Use units of (Hz, kHz, MHz, GHz). Use a different color for different type of shift and have a common legend.
- In the level structure plot, allow the user to save the plot in a user specified location

## General Tips

- Isolate calculations and gui code in the respective files.
- Make the calculation functions GUI agnostic.
- Make the gui with pyqt and marimo in the separate files as indicated by their names.
