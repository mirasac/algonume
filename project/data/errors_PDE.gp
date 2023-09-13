reset

load 'configuration.gp'

set key top right font _font
set logscale x
set offsets graph 0.05, graph 0.05, 0, 0
set grid xtics mxtics ytics mytics

set xlabel '$|x_\text{num} - x_\text{ana}|$' font _font
set xrange [7e-12:3e-2]
set xtics add ('' 1e-11, '' 1e-9, '' 1e-7, '' 1e-5, '' 1e-3, '' 1e-1)
set mxtics 10

set ylabel '$z / z_0$' font _font
set yrange [-0.5:11]
set mytics 4

plot "temperature_numerical.dat" using 5:1 every ::1 title '$T / T_0$' pt _pt ps _ps, \
	"irradiance_numerical.dat" using 7:1 every ::1 title '$E_\text{U} / S_\text{t}$' pt _pt ps _ps, \
	"irradiance_numerical.dat" using 9:1 every ::1 title '$E_\text{D} / S_\text{t}$' pt _pt ps _ps
