reset

load 'spectral_irradiance_bands.gp'

set logscale x2

set xtics nomirror

set x2label '$\lambda / (\unit{\nano\metre})$' font _font
set x2range [GPVAL_X_MIN:GPVAL_X_MAX]
set x2tics add ("100000" 1e2, "10000" 1e3, "1000" 1e4, "100" 1e5) nomirror
set mx2tics 10

replot
