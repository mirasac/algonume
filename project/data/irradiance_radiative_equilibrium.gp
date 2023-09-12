reset

load 'configuration.gp'

set title "Irradiances in radiative equilibrium" font _font
set key top right font _font
set grid xtics mxtics ytics mytics

set xlabel '$E / S_\text{t}$' font _font
set xrange [-0.05:1.7]
set mxtics 4

set ylabel '$z / z_0$' font _font
set yrange [-0.5:11]
set mytics 4

plot "irradiance_analytical.dat" using 2:1 title '$E_\text{U}$' with lines lc rgb "blue", \
	"irradiance_steady.dat" using 2:1 title ' ' pt _pt ps (2 * _ps) lc rgb "blue", \
	"+" using (sigma * T_g**4 / S_t):(z_g / z_0) title '$E_\text{U}(\delta_\text{g})$' pt 2 ps (4 * _ps) lc rgb "blue", \
	"irradiance_analytical.dat" using 3:1 title '$E_\text{D}$' with lines lc rgb "red", \
	"irradiance_steady.dat" using 3:1 title ' ' pt _pt ps (2 * _ps) lc rgb "red", \
	"+" using (sigma * T_g**4 / S_t - 1.0):(z_g / z_0) title '$E_\text{D}(\delta_\text{g})$' pt 2 ps (4 * _ps) lc rgb "red"
