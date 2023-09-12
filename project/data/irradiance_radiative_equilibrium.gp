reset

load 'configuration.gp'

set key top right font _font
set grid xtics mxtics ytics mytics

set xlabel '$E / S_\text{t}$' font _font
set xrange [-0.05:1.7]
set mxtics 4

set ylabel '$z / z_0$' font _font
set yrange [-0.5:11]
set mytics 4

plot "irradiance_analytical.dat" using 2:1 title '$E_\text{U}$' with lines lc rgb "blue", \
	"irradiance_RCM.dat" using 2:1 title '$E_\text{U}$ RCM' with lines dt _dt lc rgb "blue", \
	"irradiance_analytical.dat" using 3:1 title '$E_\text{D}$' with lines lc rgb "red", \
	"irradiance_RCM.dat" using 3:1 title '$E_\text{D}$ RCM' with lines dt _dt lc rgb "red"
