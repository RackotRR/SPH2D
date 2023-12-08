#include <iostream>
#include <stdexcept>
#include <fmt/format.h>

#include "ConsistencyCheck.h"

bool check_finite(
	shared_darray<rr_float2> r,
	shared_darray<rr_int> itype,
	shared_darray<rr_float2> v,
	shared_darray<rr_float> rho,
	shared_darray<rr_float> p)
{
	printlog_debug(__func__)();
	assert(r.get());
	assert(itype.get());
	assert(v.get());
	assert(rho.get());
	assert(p.get());

	rr_int infinite_count = 0;
	rr_uint i = 0;

#pragma omp parallel for reduction(+: infinite_count)
	for (rr_iter k = 0; k < params.ntotal; ++k) {
		if (!std::isfinite(r->at(k).x) || !std::isfinite(r->at(k).y) ||
			!std::isfinite(v->at(k).x) || !std::isfinite(v->at(k).y)) 
		{
			infinite_count++;
			i = k;
		}
		else if (!std::isfinite(rho->at(k)) ||
			!std::isfinite(p->at(k)))
		{
			infinite_count++;
			i = k;
		}
	}

	if (infinite_count > 0) {
		std::string message = fmt::format("not finite particle: \n\t type: {}\n\t r: ({}; {})\n\t v: ({}; {})\n\t rho: {}\n\t p: {}\n\t k: {}",
			itype->at(i), 
			r->at(i).x, 
			r->at(i).y, 
			v->at(i).x, 
			v->at(i).y, 
			rho->at(i), 
			p->at(i), 
			i);

		if (params.consistency_treatment == CONSISTENCY_STOP) {
			throw std::runtime_error{ message };
		}
		else {
			printlog(message)();
			std::cerr << "encounter not finite value!" << std::endl;
		}
	}

	return infinite_count == 0;
}


bool check_particles_are_within_boundaries(
	shared_darray<rr_float2> r,
	shared_darray<rr_int> itype)
{
	printlog_debug(__func__)();
	assert(r.get());
	assert(itype.get());

	rr_int count_outside_boundaries = 0;
	rr_iter outside_k = 0;

#pragma omp parallel for reduction(+: count_outside_boundaries)
	for (rr_iter k = 0; k < params.ntotal; ++k) {
		rr_float x = r->at(k).x;
		rr_float y = r->at(k).y;

		if (x > params.x_maxgeom ||
			x < params.x_mingeom ||
			y > params.y_maxgeom ||
			y < params.y_mingeom ||
			!std::isfinite(x) ||
			!std::isfinite(y))
		{
			count_outside_boundaries++;
			outside_k = k;
		}
	}

	if (count_outside_boundaries > 0) {
		std::string message = fmt::format("encounter particle outside boundaries: \n\t type: {}\n\t r: ({}; {}) // ({} .. {}; {} .. {})\n\t k: {}",
			itype->at(outside_k),
			r->at(outside_k).x,
			r->at(outside_k).y,
			params.x_mingeom,
			params.x_maxgeom,
			params.y_mingeom,
			params.y_maxgeom,
			outside_k);

		if (params.consistency_treatment == CONSISTENCY_STOP) {
			throw std::runtime_error{ message };
		}
		else {
			printlog(message)();
			std::cerr << fmt::format("encounter {} particles outside boundaries!", count_outside_boundaries) << std::endl;
		}
	}
	return count_outside_boundaries == 0;
}