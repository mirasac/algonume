# Software engineering
- A class for nuclear warheads.
- A class for scenarios, which uses instances of the class of nuclear warheads.
- Keep track of the different gases used to model the chemical composition of the atmosphere to make them easy to change (e.g. an appropriate class for gases, use specific library for computation of spectral absorption).
- Try to set the number of atmospheric layers in the RCM as a constant, so that the software can be easily adapted for a multi layer analysis.

# General
- I need probably to create skew-T plots in Gnuplot, find how. I can resort to drwaing many curves with some kind of transparency to simulate the skew-plot lines at last. Some kind of change of variables is probably needed beacuase in this case axes are still cartesian. Also search for the implementation of class metpy.plots.SkewT for help.

# Questions
- What is temperature profile instability.
- What is lapse rate? Is the technical term of temperature in pressure-temperature plot?

# Keywords
- Tracer continuity equation (cfr. https://clima.github.io/OceananigansDocumentation/v0.47.1/physics/navier_stokes_and_tracer_conservation/, https://getm.eu/files/GETM/doc/html/node78.html and https://climate.ucdavis.edu/ATM121/AtmosphericDynamics-Chapter01-Part03-Continuity.pdf).

# TODO
- Visit http://climlab.readthedocs.io/en/latest/api/climlab.radiation.Radiation.html
