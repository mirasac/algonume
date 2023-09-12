reset

load 'configuration.gp'
_logscale=2

set key top left font _font
set logscale y
set grid xtics mxtics ytics mytics

set xlabel '$\log_2(N)$' font _font
set xrange [-1:25]
set mxtics 5

set ylabel '$|x_\text{num} - x_\text{ana}|$' font _font
set mytics 5

plot "stability.dat" using (log($1)/log(_logscale)):2 title '$T(\delta_\text{g}) / T_0$' pt _pt ps _ps, \
	"stability.dat" using (log($1)/log(_logscale)):3 title '$E_\text{U}(\delta_\text{g}) / S_\text{t}$' pt _pt ps _ps, \
	"stability.dat" using (log($1)/log(_logscale)):4 title '$E_\text{D}(\delta_\text{g}) / S_\text{t}$' pt _pt ps _ps
