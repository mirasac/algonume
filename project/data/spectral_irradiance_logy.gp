reset

load 'spectral_irradiance.gp'

set key top right font _font
set logscale y

set yrange[1e-6:1]
set mytics 10

replot
