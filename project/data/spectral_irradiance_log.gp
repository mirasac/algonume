reset

_font=",10"
_pt=7

set title "Blackbody spectral irradiance" font _font
set key top left font _font
set logscale xy
set grid mxtics mytics

set xlabel '$\nu / (\unit{\per\centi\metre})$' font _font
set xrange[9e1:1.5e5]
set mxtics 10

set ylabel '$I / (\unit{\watt\metre\per\square\metre})$' font _font
set yrange[1e-10:1e-2]
set mytics 10

plot "spectral_irradiance.dat" using 1:2 notitle with lines lc rgb "blue", \
	"" using (NaN) title "Sun" lc rgb "blue" pt _pt ps 0.5, \
	"spectral_irradiance.dat" using 1:3 notitle with lines lc rgb "red", \
	"" using (NaN) title "Earth" lc rgb "red" pt _pt ps 0.5
