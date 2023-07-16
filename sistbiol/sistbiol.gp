reset

plot_file = 'sistbiol.dat'

stats plot_file using 0 nooutput
plot for [i=0:(STATS_blocks - 1)] plot_file using 2:3 index i with lines notitle
