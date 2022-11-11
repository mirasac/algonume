reset

set title "Plot phi(t)" font ",18"
set key top left font ",18"
set grid mxtics mytics

set xlabel "t [s]" font ",18"
set xrange[-1:31]
set xtic font ",9"
set mxtics 10

set ylabel "phi(t)" font ",18"
set yrange[-1:28]
set ytic font ",9"
set mytics 10

plot "elliptic_integrals.dat" using 1:2 title "Bisection" lc rgb "blue" pt 1 ps 1, \
	"elliptic_integrals.dat" using 1:3 title "Newton-Raphson" lc rgb "red" pt 1 ps 1

