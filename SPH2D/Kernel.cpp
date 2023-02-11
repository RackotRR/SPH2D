#include "CommonIncl.h"
#include "Kernel.h"

// kernel for particle itself
void kernel_self(
	rr_float& w_ii, // out, kernel for all interaction pairs
	rr_float2& dwdr_ii) // out, derivation of kernel with respect to x, y, z
{
	kernel(0.f, rr_float2{ 0.f }, w_ii, dwdr_ii);
}

void kernel(
	const rr_float2& ri,
	const rr_float2& rj,
	rr_float& w, // out, kernel for all interaction pairs
	rr_float2& dwdr) // out, derivation of kernel with respect to x, y, z
{
	rr_float2 diff = ri - rj;
	rr_float dist = length(diff);
	kernel(dist, diff, w, dwdr);
}

// calculate the smoothing kernel wij and its derivatives dwdxij
void kernel(
	const rr_float dist,
	const rr_float2& diff, // x- y- z- distance between i and j 
	rr_float& w, // out, kernel for all interaction pairs
	rr_float2& dwdx) // out, derivation of kernel with respect to x, y, z
{
	rr_float hsml{ Params::hsml };

	rr_float factor;
	rr_float q = dist / hsml;

	// zero out params
	w = 0.f;
	dwdx = { 0.f };

	if constexpr (Params::skf == 1) { // Cubic spline  
		factor = 15.f / (7.f * Params::pi * sqr(hsml)); 

		if (q <= 1) { 
			w = factor * (2.f / 3.f - q * q + q * q * q * 0.5f);
			dwdx = diff * factor * (-2.f + 3.f / 2.f * q) / sqr(hsml);
		}
		else if (q > 1 && q <= 2) {
			w = factor * (1.f / 6.f * pow(2.f - q, 3));
			dwdx = -diff / dist * factor * 1.f / 6.f * 3.f * pow(2.f - q, 2) / hsml;
		}
		else {
			w = 0.f;
			dwdx = { 0.f };
		}
	}
	else if constexpr (Params::skf == 2) { // Gauss kernel
		factor = 1.f / (pow(hsml, Params::dim) * pow(Params::pi, Params::dim / 2.f));

		if (q <= 3) {
			w = factor * exp(-q * q);
			dwdx = diff * w * (-2.f) / sqr(hsml);
		}
		else {
			w = 0.f;
			dwdx = { 0.f };
		}
	}
	else if constexpr (Params::skf == 3) { // Quintic spline 
		factor = 7.f / (478.f * Params::pi * sqr(hsml)); 

		if (q <= 1) {
			w = factor * (pow(3.f - q, 5) - 6.f * pow(2 - q, 5) + 15.f * pow(1.f - q, 5));
			dwdx = diff * factor * ((-120.f + 120.f * q - 50.f * q * q) / sqr(hsml));
		}
		else if (q > 1 && q <= 2) {
			w = factor * (pow(3.f - q, 5) - 6.f * pow(2.f - q, 5));
			dwdx = diff * factor * (-5.f * pow(3.f - q, 4) + 30.f * pow(2.f - q, 4)) / hsml / dist;
		}
		else if (q > 2 && q <= 3) {
			w = factor * pow(3.f - q, 5);
			dwdx = diff * factor * (-5.f * pow(3.f - q, 4)) / hsml / dist;
		}
		else {
			w = 0.f;
			dwdx = { 0.f };
		}

	}
}