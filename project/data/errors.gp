reset

load 'configuration.gp'

set title "Normalised errors between numerical and analytical values" font _font
set key bottom left font _font
set logscale x
set grid xtics mxtics ytics mytics

set xlabel 'Normalised error' font _font
set mxtics 10

set ylabel '$z / z_0$' font _font
set yrange [0:11]
set mytics 4

plot "errors.dat" using 2:1 title '$T / T_0$' pt _pt ps _ps, \
	"errors.dat" using 4:1 title '$E_\text{U} / S_\text{t}$' pt _pt ps _ps, \
	"errors.dat" using 5:1 title '$E_\text{D} / S_\text{t}$' pt _pt ps _ps
