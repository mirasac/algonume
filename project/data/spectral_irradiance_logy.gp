reset

load 'configuration.gp'
nu_threshold=2.156e3

set title "Blackbody spectral irradiance" font _font
set key top left font _font
set logscale xy
set grid mxtics mytics

set xlabel '$\nu / (\unit{\per\centi\metre})$' font _font
set xrange[9e1:1.5e5]
set mxtics 10

set ylabel '$I / (\unit{\watt\centi\metre\per\square\metre})$' font _font
set yrange[1e-8:1]
set mytics 10

set arrow from nu_threshold, graph 0 to nu_threshold, graph 1 nohead lc rgb "black"
set xtics add ('$\nu_\text{threshold}$' nu_threshold)

plot "spectral_irradiance.dat" using 1:2 notitle with lines lc rgb "blue", \
	"" using (NaN) title "Sun" lc rgb "blue" pt _pt ps 0.5, \
	"spectral_irradiance.dat" using 1:3 notitle with lines lc rgb "red", \
	"" using (NaN) title "Earth" lc rgb "red" pt _pt ps 0.5
