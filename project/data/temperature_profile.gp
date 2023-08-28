# Useful references:
# - https://gnuplot.sourceforge.net/demo/multiaxis.html;
# - https://stackoverflow.com/questions/2827650/plotting-two-axes-in-gnuplot;
# - https://stackoverflow.com/questions/63771600/is-there-a-way-to-have-3-different-y-axes-on-one-graph-using-gnuplot.

reset

load 'configuration.gp'

set title "Temperature vertical profile" font _font
set key top right font _font
set logscale y2
#set grid mxtics mytics my2tics

set xlabel '$T / \unit{\kelvin}$' font _font
#set xrange[-10.0:100.0]  # / K
set mxtics 5

set ylabel '$z / \unit{\metre}$' font _font
#set yrange[0.0:45000.0]  # / m
set ytics nomirror
set mytics 5

set y2label '$P / \unit{\pascal}$' font _font
#set y2range[1e5:1e2]  # / Pa
set y2tics nomirror
set my2tics 10

stats "temperature.dat" using 1 nooutput
index_last=STATS_blocks - 1

plot "temperature.dat" using ($3+$1):2 index 0 axes x1y1 title 'Initial condition' with lines lc rgb "0x000000FF", \
	"temperature.dat" using ($3+$1):4 index 0 axes x1y2 notitle with lines lc rgb "0xFF0000FF", \
	"temperature.dat" using ($3+$1):2 index index_last axes x1y1 title 'Steady state' with lines lc rgb "0x00FF0000", \
	"temperature.dat" using ($3+$1):4 index index_last axes x1y2 notitle with lines lc rgb "0xFFFF0000"
