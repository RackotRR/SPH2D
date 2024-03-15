#include "CommonIncl.h"
#include "Kernel.h"

void cubic_kernel(rr_float dist, const rr_float2& diff, rr_float& w, rr_float2& dwdr) {
	w = cubic_kernel_w(dist);
	dwdr = cubic_kernel_dwdr(dist, diff);
}

void gauss_kernel(rr_float dist, const rr_float2& diff, rr_float& w, rr_float2& dwdr) {
	w = gauss_kernel_w(dist);
	dwdr = gauss_kernel_dwdr(dist, diff);
}

void wendland_kernel(rr_float dist, const rr_float2& diff, rr_float& w, rr_float2& dwdr) {
	w = wendland_kernel_w(dist);
	dwdr = wendland_kernel_dwdr(dist, diff);
}

void desbrun_kernel(rr_float dist, const rr_float2& diff, rr_float& w, rr_float2& dwdr) {
	w = desbrun_kernel_w(dist);
	dwdr = desbrun_kernel_dwdr(dist, diff);
}

// kernel for particle itself
void kernel_self(
	rr_float& w_ii, // out, kernel for all interaction pairs
	rr_float2& dwdr_ii, // out, derivation of kernel with respect to x, y, z
	rr_uint skf)
{
	kernel(0.f, rr_float2{ 0.f }, w_ii, dwdr_ii, skf);
}

void kernel(
	const rr_float2& ri,
	const rr_float2& rj,
	rr_float& w, // out, kernel for all interaction pairs
	rr_float2& dwdr, // out, derivation of kernel with respect to x, y, z
	rr_uint skf)
{
	rr_float2 diff = ri - rj;
	rr_float dist = length(diff);
	kernel(dist, diff, w, dwdr, skf);
}

// calculate the smoothing kernel wij and its derivatives dwdxij
void kernel(
	const rr_float dist,
	const rr_float2& diff, // x- y- distance between i and j 
	rr_float& w, // out, kernel for all interaction pairs
	rr_float2& dwdr, // out, derivation of kernel
	rr_uint skf)
{
	switch (skf) {
	case 1: cubic_kernel(dist, diff, w, dwdr); break;
	case 2: gauss_kernel(dist, diff, w, dwdr); break;
	case 3: wendland_kernel(dist, diff, w, dwdr); break;
	case 4: desbrun_kernel(dist, diff, w, dwdr); break;
	default: cubic_kernel(dist, diff, w, dwdr); break;
	}
}

rr_float kernel_w(rr_float dist, rr_uint skf) {
	switch (skf) {
	case 1: return cubic_kernel_w(dist);
	case 2: return gauss_kernel_w(dist);
	case 3: return wendland_kernel_w(dist);
	case 4: return desbrun_kernel_w(dist);
	default: return cubic_kernel_w(dist); 
	}
}

rr_float2 kernel_dwdr(rr_float dist, const rr_float2& diff, rr_uint skf) {
	switch (skf) {
	case 1: return cubic_kernel_dwdr(dist, diff);
	case 2: return gauss_kernel_dwdr(dist, diff);
	case 3: return wendland_kernel_dwdr(dist, diff);
	case 4: return desbrun_kernel_dwdr(dist, diff);
	default: return cubic_kernel_dwdr(dist, diff);
	}
}

void calculate_kernels(
	const rr_uint ntotal,
	const heap_darray<rr_float2>& r,
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	heap_darray_md<rr_float>& w, // precomputed kernel
	heap_darray_md<rr_float2>& dwdr, // precomputed kernel derivative
	rr_uint skf)
{
#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; j++) { // run through all particles
		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != ntotal; // particle near
			++n)
		{
			rr_float2 diff = r(i) - r(j);
			rr_float dist = length(diff);

			kernel(dist, diff, w(n, j), dwdr(n, j), skf);
		}
	}
}
void calculate_kernels_w(
	const rr_uint ntotal,
	const heap_darray<rr_float2>& r,
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	heap_darray_md<rr_float>& w, // precomputed kernel
	rr_uint skf)
{
	printlog_debug("calculate_kernels_w: ")(skf)();

#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; j++) { // run through all particles
		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != ntotal; // particle near
			++n)
		{
			rr_float2 diff = r(i) - r(j);
			rr_float dist = length(diff);

			w(n, j) = kernel_w(dist, skf);
		}
	}
}
void calculate_kernels_dwdr(
	const rr_uint ntotal,
	const heap_darray<rr_float2>& r,
	const heap_darray_md<rr_uint>& neighbours, // neighbours indices
	heap_darray_md<rr_float2>& dwdr, // precomputed kernel
	rr_uint skf)
{
	printlog_debug("calculate_kernels_dwdr: ")(skf)();
	
#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; j++) { // run through all particles
		rr_uint i;
		for (rr_iter n = 0;
			i = neighbours(n, j), i != ntotal; // particle near
			++n)
		{
			rr_float2 diff = r(i) - r(j);
			rr_float dist = length(diff);

			dwdr(n, j) = kernel_dwdr(dist, diff, skf);
		}
	}
}