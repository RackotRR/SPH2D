#include "CommonIncl.h"


// calculate the smoothing kernel wij and its derivatives dwdxij
void kernel(
	const double r, // distance between particles i and j
	const heap_array<double, Params::dim>& dx, // x- y- z- distance between i and j 
	double& w, // out, kernel for all interaction pairs
	heap_array<double, Params::dim>& dwdx) // out, derivation of kernel with respect to x, y, z
{
	double hsml{ Params::hsml };

	double factor;
	double q{ r / hsml };

	// zero out params
	w = 0;
	for (size_t d{}; d < Params::dim; d++) {
		dwdx(d) = 0;
	}

	if (Params::skf == 1) { // Cubic spline  
		factor = 15.0 / (7.0 * Params::pi * sqr(hsml)); 

		if (q >= 0 && q <= 1) { 
			w = factor * (2.0 / 3.0 - q * q + q * q * q * 0.5);
			for (size_t d{}; d < Params::dim; d++) {
				dwdx(d) = factor * (-2.0 + 3.0 / 2.0 * q) / sqr(hsml) * dx(d);
			}
		}
		else if (q > 1 && q <= 2) {
			w = factor * (1.0 / 6.0 * pow(2.0 - q, 3));
			for (size_t d{}; d < Params::dim; d++) {
				dwdx(d) = -factor * 1.0 / 6.0 * 3.0 * pow(2.0 - q, 2) / hsml * (dx(d) / r);
			}
		}
		else {
			w = 0;
			for (size_t d{}; d < Params::dim; d++) {
				dwdx(d) = 0;
			}
		}
	}
	else if (Params::skf == 2) { // Gauss kernel
		factor = 1.0 / (pow(hsml, Params::dim) * pow(Params::pi, Params::dim / 2.0));

		if (q >= 0 && q <= 3) {
			w = factor * exp(-q * q);
			for (size_t d{}; d < Params::dim; d++) {
				dwdx(d) = w * (-2.0) * dx(d) / sqr(hsml);
			}
		}
		else {
			w = 0;
			for (size_t d{}; d < Params::dim; d++) {
				dwdx(d) = 0;
			}
		}
	}
	else if (Params::skf == 3) { // Quintic spline 
		factor = 7.0 / (478.0 * Params::pi * sqr(hsml)); 

		if (q >= 0 && q <= 1) {
			w = factor * (pow(3.0 - q, 5) - 6.0 * pow(2 - q, 5) + 15 * pow(1 - q, 5));
			for (size_t d{}; d < Params::dim; d++) {
				dwdx(d) = factor * ((-120.0 + 120 * q - 50 * q * q) / sqr(hsml) * dx(d));
			}
		}
		else if (q > 1 && q <= 2) {
			w = factor * (pow(3 - q, 5) - 6 * pow(2 - q, 5));
			for (size_t d{}; d < Params::dim; d++) {
				dwdx(d) = factor * (-5 * pow(3 - q, 4) + 30 * pow(2 - q, 4)) / hsml * dx(d) / r;
			}
		}
		else if (q > 2 && q <= 3) {
			w = factor * pow(3 - q, 5);
			for (size_t d{}; d < Params::dim; d++) {
				dwdx(d) = factor * (-5 * pow(3 - q, 4)) / hsml * dx(d) / r;
			}
		}
		else {
			w = 0;
			for (size_t d{}; d < Params::dim; d++) {
				dwdx(d) = 0;
			}
		}

	}
}