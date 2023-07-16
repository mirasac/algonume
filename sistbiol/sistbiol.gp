reset

plot_file = "sistbiol.dat"

set grid mxtics mytics

set xlabel "x" font ",18"
set xrange[-0.05:1.05]
set xtic font ",9"
set mxtics 10

set ylabel "y" font ",18"
set yrange[-0.05:1.05]
set ytic font ",9"
set mytics 10

stats plot_file using 0 nooutput
plot for [i=0:(STATS_blocks - 1)] plot_file using 2:3 index i with lines notitle
