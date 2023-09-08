reset

load 'configuration.gp'

set title "Irradiances in radiative equilibrium" font _font
set key top right font _font
set grid xtics mxtics ytics mytics

set xlabel '$E / S_\text{t}$' font _font
set xrange [-0.05:1.65]
set mxtics 4

set ylabel '$z / z_0$' font _font
set yrange [0:11]
set mytics 4

plot "irradiance_analytical.dat" using 2:1 title 'Upward' with lines lc rgb "blue", \
	"irradiance_analytical.dat" using 3:1 title 'Downward' with lines lc rgb "red"
