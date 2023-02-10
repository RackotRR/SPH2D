#pragma once
#include "CommonIncl.h"

// calculate the smoothing kernel wij and its derivatives dwdxij
void kernel(
	const rr_float2& dx, // x- y- z- distance between i and j 
	rr_float& w, // out, kernel for all interaction pairs
	rr_float2& dwdx); // out, derivation of kernel with respect to x, y, z

