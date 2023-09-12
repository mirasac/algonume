# General
- I need probably to create skew-T plots in Gnuplot, find how. I can resort to drwaing many curves with some kind of transparency to simulate the skew-plot lines at last. Some kind of change of variables is probably needed because in this case axes are still cartesian. Also search for the implementation of class metpy.plots.SkewT for help.
- Use the sigma adimensional pressure variable to represent heights.
- Introduce spectral dependence in albedo values for each layer surface.
- TTAPS-I works on visible wavelength of 550 nm.
- Gravitational acceleration is not constant and depends also on the altitude. Discuss why I take a constant value, or discuss the parametrization I choose if I take a non-constant value.
- Distinguish thorougly simulation time and computation time.
- Data file format for output values is set as shown in the [official Gnuplot documentation](http://gnuplot.info/docs_5.5/loc17750.html) but it is different from what suggested by Professor's (i.e. y variable fixed in each block), so I might change the format.
- Before launching the simulation, normalize quantities to reduce their magnitudes and avoid numerical divergences (e.g. time normalized in units of 1 day).
- Divide total solar irradiance by 4 due to combination of day-night cycle and distribution of total energy over a disk (i.e. S_0 * r * pi^2 / (4 * r * pi^2) = S_0 / 4) and again by 2 due to the choice to focus on the Northern Hemisphere and not just on the illuminated face of the Earth.
- Use optical depths and formulae for soot from TTAPS-II because reproducing the aerosol formation computational model presented in Turco 1979 goes beyond the purpose of this work.
- Spectral radiance is also called specific intensity, mainly in the old literature.
- Use non-dimensional optical coordinates, cfr. \cite[290]{Modest}.
- Energy is not conserved in the atmosphere, but radiative equilibrium is forced at TOA.
- What is lapse rate? Is the technical term of temperature in pressure-temperature plot? Sort of: it is the variation of temperature with altitude and since altitude and pressure can be related quantitatively, it can be expressed also as variation with pressure.
- Absorption (or more generally, exctintion) coefficients depend on temperature and pressure. How can I use them to approximate absorption in uniform layers of gases? I must assume that their variation with temperature and pressure is negligible with respect to the varying quantities I am evaluating.
- Why there is a pi multiplying the spectrally averaged plank function? Probably for the parametrizations I use it is not needed. It is needed, since it is the factor resulting from the integration in the solid angle of an hemisphere of Planck law for the azimuthal component of the solar beam.

## Questions
- What is temperature profile instability.
- Do I add evaluation of quantities at ground level in a separate variable or as element at index N + 1 of arrays?

## Keywords
- Tracer continuity equation (cfr. https://clima.github.io/OceananigansDocumentation/v0.47.1/physics/navier_stokes_and_tracer_conservation/, https://getm.eu/files/GETM/doc/html/node78.html and https://climate.ucdavis.edu/ATM121/AtmosphericDynamics-Chapter01-Part03-Continuity.pdf).

# Code

## Software engineering
- A class for nuclear warheads.
- A class for scenarios, which uses instances of the class of nuclear warheads.
- Keep track of the different gases used to model the chemical composition of the atmosphere to make them easy to change (e.g. an appropriate class for gases, use specific library for computation of spectral absorption).
- Try to set the number of atmospheric layers in the RCM as a constant, so that the software can be easily adapted for a multi layer analysis.
- Structures or classes? In C++ there is no difference in efficiency, just in the default visibility of members.
- Atmospheric layers and their thickness can be stored in an array each, with size fixed a priori by a constant N_LAYERS. The program should adapt automatically to the chosen number of layers.
- A class for atmospheric layers, where the values are functions depending on the spectral propertries (e.g. albedo depending on the wavelength).
- As flavour, do a routine in the main module where the spectral bands are evaluated automatically based on the studied gases, instead of specifying them manually in an array. It can be also a class method of the atmopsheric layer class.
- Where possible, use integer values.
- Flush file buffer only when needed because it could become a huge bottleneck in program execution.
- How do I choose which formulae are suitable for specific wavelengths? Do I simply set a threshold on wavelenght value? I may also set the function as method of the class used to represent gases.

## Release
- Remove all personal comments.
- Remove all useless functions (e.g. in `main.cpp` array dz).

# Report
- For section label markers I put the name of the section, case sensitive, after the string `sec:`. This way the marker could contain spaces, but this should not be a drawback since nowadays most packages managing labels support spaces in them. Another possible drawback is the length that some labels can reach. Instead the advantages are many: it is a consistent way to write meaningful labels for sections, labels get immediately fixed when section titles are changed (i.e. by using a find and replace tool in any editor).
- Use nominal values then associate relative error to them.
- The chosen font size is 10pt everywhere. Set it automatically everywhere based on a single configuration file takes too much effort compared by setting it manually.
- I do not use the sigma coordinate as vertical coordinate because in the model I am considering it is simply a rescaled pressure. It could be useful if this model is then used in a model where the horizontal profile of surface pressure changes.
- Show equation dependencies only when is needed for clarity.

## Release
- Check verbal tenses and conjugations.
- Remove all personal comments.
- Remove unused appendices, maybe keep them as separate LaTeX source files and include only the ones needed.
- Correct format issues related to composition of bibliography.
- Align fields of tables in source code.
- Use cf. instead of cfr. and use see when needed.
- Check content against what prof. Palazzi taught us during lesson.
- In the hypotheses stay general by using the term medium instead of layer. Then introduce layers and use them in the transition to the numerical resolution of the model.
- Check the correctness in the use of terms equilibrium and steady state.
- Check and integrate captions.
- Move figure titles to the captions.
