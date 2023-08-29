reset

load 'configuration.gp'
nu_div=2.154e3 # 1 / cm

set title "Spectral irradiance at Earth's surface" font _font
set key top left font _font
set logscale x
set grid xtics mxtics ytics mytics

set xlabel '$\nu / (\unit{\per\centi\metre})$' font _font
set xrange [9e1:1.5e5] # 1 / cm
set mxtics 10

set ylabel '$E / (\unit{\watt\centi\metre\per\square\metre})$' font _font
set yrange [0:4.5e-1]
set mytics 2

set arrow from nu_div, graph 0 to nu_div, graph 1 nohead lc rgb "black"
set xtics add ('$\nu_\text{div}$' nu_div)

plot "spectral_irradiance.dat" using 1:2 notitle with lines lc rgb "blue", \
	"" using (NaN) title "Sun" lc rgb "blue" pt _pt ps 0.5, \
	"spectral_irradiance.dat" using 1:3 notitle with lines lc rgb "red", \
	"" using (NaN) title "Earth" lc rgb "red" pt _pt ps 0.5
