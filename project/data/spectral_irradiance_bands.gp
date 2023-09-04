reset

load 'spectral_irradiance.gp'

# Show EM spectral bands from CIE S 017:2020.
nu_IR_VIS=1e5 / 7.8 # 1 / cm
nu_VIS_UV=1e5 / 4.0 # 1 / cm
set object 1 rectangle from graph 0, 0 to nu_IR_VIS, graph 1 behind fc rgb "0xCCFF0000" fs noborder
set object 2 rectangle from nu_IR_VIS, graph 0 to nu_VIS_UV, graph 1 behind fc rgb "0xCC00FF00" fs noborder
set object 3 rectangle from nu_VIS_UV, graph 0 to graph 1, 1 behind fc rgb "0xCC0000FF" fs noborder

set label "IR" at nu_IR_VIS, graph 0.95 right offset -1, 0
set label "VIS" at nu_VIS_UV, graph 0.95 right offset -1, 0
set label "UV" at graph 1, 0.95 right offset -1, 0

replot
