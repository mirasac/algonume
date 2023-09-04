#ifndef CONVECTION_H
#define CONVECTION_H

/*
Apply convective adjustment to the temperatures in K contained in T using value lapse_rate / (K / m) and altitudes in array z expressed in unit m.

Both arrays T and z have length n_layer + 1 where the element in position n_layer is referred to the ground level.
*/
void convective_adjustment(double const lapse_rate, double const n_layer, double T[], double const z[]);

#endif /* CONVECTION_H */
