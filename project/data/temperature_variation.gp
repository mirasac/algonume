reset

load 'configuration.gp'

set title "Temperature variation" font _font
unset key
set logscale y2

set view map
unset surface
set contour base
set cntrlabel start 2 font _font

set xlabel '$t / \unit{\hour}$' font _font
#set xrange[0.0:85441.0]  # / h
set mxtics 5

set ylabel '$z / \unit{\metre}$' font _font
#set yrange[0.0:45000.0]  # / m
set mytics 5

set y2label '$P / \unit{\pascal}$' font _font
#set y2range[1e5:1e2]  # / Pa
set my2tics 10

splot "temperature.dat" using 1:2:3, \
	"temperature.dat" using 1:4:3 lc rgb "0xFF000000"
