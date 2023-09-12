reset

load 'configuration.gp'

set title "Errors between numerical and analytical values" font _font
set key top right font _font
set logscale x
set offsets graph 0.05, graph 0.05, 0, 0
set grid xtics mxtics ytics mytics

set xlabel '$|x_\text{num} - x_\text{ana}|$' font _font
set mxtics 10

set ylabel '$z / z_0$' font _font
set yrange [-0.5:11]
set mytics 4

plot "errors_steady.dat" using 2:1 title '$T$ 2' pt _pt ps _ps, \
	"errors_steady.dat" using 4:1 title '$E_\text{U} / S_\text{t}$ 2' pt _pt ps _ps, \
	"errors_steady.dat" using 5:1 title '$E_\text{D} / S_\text{t}$ 2' pt _pt ps _ps
