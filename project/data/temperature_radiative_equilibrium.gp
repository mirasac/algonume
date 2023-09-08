reset

load 'configuration.gp'

set title "Temperature profile in radiative equilibrium" font _font
set key top right font _font
set grid xtics mxtics ytics mytics

set xlabel '$T / \unit{\kelvin}$' font _font
set xrange [210:264] # K
set mxtics 5

set ylabel '$z / z_0$' font _font
set yrange [0:11]
set mytics 4

plot "temperature_analytical.dat" using 2:1 notitle with lines lc rgb "blue"
