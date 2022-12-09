reset

set title "Trajectory" font ",18"
set key top left font ",18"
set grid mxtics mytics

set xlabel "x [m]" font ",18"
set xrange[-0.1:10.1]
set xtic font ",9"
#set mxtics 10

set ylabel "y(x) [m]" font ",18"
set yrange[-0.1:5]
set ytic font ",9"
#set mytics 10

plot "projectile.dat" using 1:2 title "Initial angle phi_0" lc rgb "blue" pt 1 ps 1, \
	"projectile.dat" using 1:3 title "Initial angle phi_1" lc rgb "red" pt 1 ps 1

