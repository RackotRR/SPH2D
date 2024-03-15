#pragma once
#include "CommonIncl.h"
#include "SmoothingKernel.h"

void calculate_kernels_w(
	const rr_uint ntotal,
	const heap_darray<rr_float2>& r,
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	heap_darray_md<rr_float>& w, // precomputed kernel
	rr_uint skf);
void calculate_kernels_dwdr(
	const rr_uint ntotal,
	const heap_darray<rr_float2>& r,
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	heap_darray_md<rr_float2>& dwdr, // precomputed kernel
	rr_uint skf);

// calculate the smoothing kernel wij and its derivatives dwdxij
void kernel(
	const rr_float dist,
	const rr_float2& diff, // x- y- distance between i and j 
	rr_float& w, // out, kernel for all interaction pairs
	rr_float2& dwdr, // out, derivation of kernel with respect to x, y
	rr_uint skf);


// calculate the smoothing kernel wij and its derivatives dwdxij
void kernel(
	const rr_float2& ri,
	const rr_float2& rj,
	rr_float& w, // out, kernel for all interaction pairs
	rr_float2& dwdr, // out, derivation of kernel with respect to x, y
	rr_uint skf);

// kernel for particle itself
void kernel_self(
	rr_float& w, // out, kernel for all interaction pairs
	rr_float2& dwdr, // out, derivation of kernel with respect to x, y
	rr_uint skf);

rr_float kernel_w(rr_float dist, rr_uint skf);
rr_float2 kernel_dwdr(rr_float dist, const rr_float2& diff, rr_uint skf);

// specific kernel functions:

void cubic_kernel(rr_float dist, const rr_float2& diff, rr_float& w, rr_float2& dwdr);

void gauss_kernel(rr_float dist, const rr_float2& diff, rr_float& w, rr_float2& dwdr);

void wendland_kernel(rr_float dist, const rr_float2& diff, rr_float& w, rr_float2& dwdr);

void desbrun_kernel(rr_float dist, const rr_float2& diff, rr_float& w, rr_float2& dwdr);