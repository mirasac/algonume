#include "convection.h"

void convective_adjustment(double const lapse_rate, double const n_layer, double T[], double const z[]) {
	for (int i = n_layer - 1; i > 0; i--) {
		if ((T[i] - T[i + 1]) / (z[i] - z[i + 1]) < - lapse_rate) {
			T[i] = T[i + 1] - (z[i] - z[i + 1]) * lapse_rate;
		}
	}
}
