reset

load 'configuration.gp'

set key top right font _font
set grid xtics mxtics ytics mytics

set xlabel '$T / T_0$' font _font
set xrange [0.8:1.15]
set mxtics 5

set ylabel '$z / z_0$' font _font
set yrange [-0.5:11]
set mytics 4

plot "temperature_analytical.dat" using 2:1 title '$T$' with lines lc rgb "red", \
	"temperature_RCM.dat" using 2:1 title '$T$ RCM' with lines dt _dt lc rgb "red", \
	"+" using (T_g / T_0):(z_g / z_0) title '$T_\text{g}$' pt 2 ps (4 * _ps) lc rgb "red", \
