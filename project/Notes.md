# Software engineering
- A class for nuclear warheads.
- A class for scenarios, which uses instances of the class of nuclear warheads.
- Keep track of the different gases used to model the chemical composition of the atmosphere to make them easy to change (e.g. an appropriate class for gases, use specific library for computation of spectral absorption).
- Try to set the number of atmospheric layers in the RCM as a constant, so that the software can be easily adapted for a multi layer analysis.
- Structures or classes? In C++ there is no difference in efficiency, just in the default visibility of members.
- Atmospheric layers and their thickness can be stored in an array each, with size fixed a priori by a constant N_LAYERS. The program should adapt automatically to the chosen number of layers.
- A class for atmospheric layers, where the values are functions depending on the spectral propertries (e.g. albedo depending on the wavelength).
- As flavour, do a routine in the main module where the spectral bands are evaluated automatically based on the studied gases, instead of specifying them manually in an array.

# General
- I need probably to create skew-T plots in Gnuplot, find how. I can resort to drwaing many curves with some kind of transparency to simulate the skew-plot lines at last. Some kind of change of variables is probably needed because in this case axes are still cartesian. Also search for the implementation of class metpy.plots.SkewT for help.
- Use the sigma adimensional pressure variable to represent heights.
- Introduce spectral dependence in albedo values for each layer surface.
- TTAPS-I works on visible wavelength of 550 nm.
- Gravitational acceleration is not constant and depends also on the altitude. Discuss why I take a constant value, or discuss the parametrization I choose if I take a non-constant value.
- Distinguish thorougly simulation time and computation time.
- Data file format for output values is set as shown in the [official Gnuplot documentation](http://gnuplot.info/docs_5.5/loc17750.html) but it is different from what suggested by Professor's (i.e. y variable fixed in each block), so I might change the format.
- Before launching the simulation, normalize quantities to reduce their magnitudes and avoid numerical divergences (e.g. time normalized in units of 1 day).

# Questions
- What is temperature profile instability.
- What is lapse rate? Is the technical term of temperature in pressure-temperature plot?
- How do I choose which formulae are suitable for specific wavelengths? Do I simply set a threshold on wavelenght value? I may also set the function as method of the class used to represent gases.
- Why there is a pi multiplying the spectrally averaged plank function? Probably for the parametrizations I use it is not needed.

# Keywords
- Tracer continuity equation (cfr. https://clima.github.io/OceananigansDocumentation/v0.47.1/physics/navier_stokes_and_tracer_conservation/, https://getm.eu/files/GETM/doc/html/node78.html and https://climate.ucdavis.edu/ATM121/AtmosphericDynamics-Chapter01-Part03-Continuity.pdf).

# Release
- Check verbal tenses and conjugations.
- Remove all personal comments.
- Remove unused appendices, maybe keep them as separate LaTeX source files and include only the ones needed.
