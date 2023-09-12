reset

load 'configuration.gp'

set title "Temperature profile in radiative equilibrium" font _font
set key top right font _font
set grid xtics mxtics ytics mytics

set xlabel '$T / \unit{\kelvin}$' font _font
set xrange [210:290] # K
set mxtics 5

set ylabel '$z / z_0$' font _font
set yrange [-0.5:11]
set mytics 4

plot "temperature_analytical.dat" using 2:1 title '$T$' with lines lc rgb "red", \
	"temperature_steady.dat" using 2:1 title ' ' pt _pt ps (2 * _ps) lc rgb "red", \
	"+" using (T_g):(z_g / z_0) title '$T_\text{g}$' pt 2 ps (4 * _ps) lc rgb "red", \
