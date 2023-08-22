reset

load 'configuration.gp'
nu_div=2.156e3

set title "Spectral irradiance at Earth's surface" font _font
set key top right font _font
set logscale x
set grid mxtics mytics

set xlabel '$\nu / (\unit{\per\centi\metre})$' font _font
set xrange[9e1:1.5e5]
set mxtics 10

set ylabel '$I / (\unit{\watt\centi\metre\per\square\metre})$' font _font
set yrange[0:4.5e-1]
set mytics 5

set arrow from nu_div, graph 0 to nu_div, graph 1 nohead lc rgb "black"
set xtics add ('$\nu_\text{div}$' nu_div)

# Show EM spectrum subdivision from ISO 21348:2007.
#nu_IR_VIS=1e5 / 7.6
#nu_VIS=1e5 / 3.8
#nu_VIS_UV=2.5e4
#set object 1 rectangle from graph 0, 0 to nu_IR_VIS, graph 1 behind fc rgb "0xCCFF0000" fs noborder
#set object 2 rectangle from nu_IR_VIS, 0 to nu_VIS, graph 1 behind fc rgb "0xCC00FF00" fs noborder
#set object 3 rectangle from nu_VIS_UV, 0 to graph 1, 1 behind fc rgb "0xCC0000FF" fs noborder

plot "spectral_irradiance.dat" using 1:2 notitle with lines lc rgb "blue", \
	"" using (NaN) title "Sun" lc rgb "blue" pt _pt ps 0.5, \
	"spectral_irradiance.dat" using 1:3 notitle with lines lc rgb "red", \
	"" using (NaN) title "Earth" lc rgb "red" pt _pt ps 0.5
