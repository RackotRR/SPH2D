#pragma once
#include "CommonIncl.h"

// calculate the smoothing kernel wij and its derivatives dwdxij
void kernel(
	const rr_float2& diff, // x- y- z- distance between i and j 
	rr_float& w, // out, kernel for all interaction pairs
	rr_float2& dwdr); // out, derivation of kernel with respect to x, y, z


// calculate the smoothing kernel wij and its derivatives dwdxij
void kernel(
	const rr_float2& ri, 
	const rr_float2& rj, 
	rr_float& w, // out, kernel for all interaction pairs
	rr_float2& dwdr); // out, derivation of kernel with respect to x, y, z

// kernel for self
void kernel_self(
	rr_float& w, // out, kernel for all interaction pairs
	rr_float2& dwdr); // out, derivation of kernel with respect to x, y, z
