reset

load 'spectral_irradiance_bands.gp'

set logscale x2

set xrange [nu_min:nu_max] # 1 / cm
set xtics nomirror

set x2label '$\lambda / (\unit{\nano\metre})$' font _font
set x2range [1e7 / nu_min:1e7 / nu_max] # / nm
set x2tics nomirror
set mx2tics 10

replot
