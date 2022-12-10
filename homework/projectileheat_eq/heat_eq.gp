reset

set title "Temperature distribution" font ",18"
set key top right font ",18"
set grid mxtics mytics

set xlabel "x" font ",18"
set xrange[-1.05:1.05]
set xtic font ",9"
set mxtics 10

set ylabel "T(x)" font ",18"
set yrange[0.99:1.5]
set ytic font ",9"
set mytics 10

plot "heat_eq.dat" index 0 notitle lc rgb "blue" pt 7 ps 0.1, \
	"" using (NaN) title "t = 0.05" lc rgb "blue" pt 7 ps 1, \
	"heat_eq.dat" index 1 notitle lc rgb "red" pt 7 ps 0.1, \
	"" using (NaN) title "t = 0.1" lc rgb "red" pt 7 ps 1
