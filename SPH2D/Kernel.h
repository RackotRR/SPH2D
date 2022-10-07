#pragma once
#include "CommonIncl.h"

// calculate the smoothing kernel wij and its derivatives dwdxij
void kernel(
	const double r, // distance between particles i and j
	const stack_array<double, Params::dim>& dx, // x- y- z- distance between i and j 
	double& w, // out, kernel for all interaction pairs
	stack_array<double, Params::dim>& dwdx); // out, derivation of kernel with respect to x, y, z

