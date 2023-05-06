#pragma once
#include "CommonIncl.h"
#include "SmoothingKernel.h"

// calculate the smoothing kernel wij and its derivatives dwdxij
void kernel(
	const rr_float dist,
	const rr_float2& diff, // x- y- distance between i and j 
	rr_float& w, // out, kernel for all interaction pairs
	rr_float2& dwdr, // out, derivation of kernel with respect to x, y
	rr_uint skf = params.skf);


// calculate the smoothing kernel wij and its derivatives dwdxij
void kernel(
	const rr_float2& ri,
	const rr_float2& rj,
	rr_float& w, // out, kernel for all interaction pairs
	rr_float2& dwdr, // out, derivation of kernel with respect to x, y
	rr_uint skf = params.skf);

// kernel for particle itself
void kernel_self(
	rr_float& w, // out, kernel for all interaction pairs
	rr_float2& dwdr, // out, derivation of kernel with respect to x, y
	rr_uint skf = params.skf);

rr_float kernel_w(rr_float dist, rr_uint skf = params.skf);
rr_float2 kernel_dwdr(rr_float dist, const rr_float2& diff, rr_uint skf = params.skf);

// specific kernel functions:

void cubic_kernel(rr_float dist, const rr_float2& diff, rr_float& w, rr_float2& dwdr);

void gauss_kernel(rr_float dist, const rr_float2& diff, rr_float& w, rr_float2& dwdr);

void quintic_kernel(rr_float dist, const rr_float2& diff, rr_float& w, rr_float2& dwdr);

void desbrun_kernel(rr_float dist, const rr_float2& diff, rr_float& w, rr_float2& dwdr);