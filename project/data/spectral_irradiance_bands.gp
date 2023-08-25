reset

load 'spectral_irradiance.gp'

# Show EM spectrum subdivision from ISO 21348:2007.
nu_IR_VIS=1e5 / 7.6  # / (1 / cm)
nu_VIS=1e5 / 3.8  # / (1 / cm)
nu_VIS_UV=2.5e4  # / (1 / cm)
set object 1 rectangle from graph 0, 0 to nu_IR_VIS, graph 1 behind fc rgb "0xCCFF0000" fs noborder
set object 2 rectangle from nu_IR_VIS, 0 to nu_VIS, graph 1 behind fc rgb "0xCC00FF00" fs noborder
set object 3 rectangle from nu_VIS_UV, 0 to graph 1, 1 behind fc rgb "0xCC0000FF" fs noborder

replot
