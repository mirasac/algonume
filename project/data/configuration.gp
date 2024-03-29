_font=",10"
_pt=7
_ps=0.2
_dt=2

A=0.3
Gamma=6.5e-3 # K / m
sigma=5.670374419e-8 # W / (m^2 K^4)
S_0=1361.0 # W / m^2
S_t=0.25 * (1.0 - A) * S_0 # W / m^2
T_g=288.15 # K
z_0=2000.0 # m
z_g=0.0 # m

T_0=exp(0.25 * log(S_t / sigma))
