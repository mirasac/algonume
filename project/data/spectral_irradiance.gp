reset

v_font=",12"

set title "Blackbody spectral irradiance" font v_font
set key top right font v_font
set logscale x
set grid mxtics mytics

set xlabel '$\nu / (\unit{\per\centi\metre})$' font v_font
set xrange[9e1:1.5e5]
set mxtics 10

set ylabel '$I / (\unit{\watt\metre\per\square\metre})$' font v_font
set yrange[-5e5:4e7]
set mytics 10

plot "spectral_irradiance.dat" using 1:2 notitle lc rgb "blue" pt 7 ps 0.1, \
	"" using (NaN) title "Sun" lc rgb "blue" pt 7 ps 1, \
	"spectral_irradiance.dat" using 1:3 notitle lc rgb "red" pt 7 ps 0.1, \
	"" using (NaN) title "Earth" lc rgb "red" pt 7 ps 1
