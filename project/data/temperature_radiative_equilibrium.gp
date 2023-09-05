reset

load 'configuration.gp'

set title "Temperature profile in radiative equilibrium" font _font
set key top left font _font
set grid xtics mxtics ytics mytics

set xlabel '$T / \unit{\kelvin}$' font _font
set xrange [210:270] # K
set mxtics 5

set ylabel '$z / z_0$' font _font
set yrange [0:12]
set mytics 5

plot "temperature_radiative_equilibrium.dat" using 2:1 notitle with lines lc rgb "blue"
