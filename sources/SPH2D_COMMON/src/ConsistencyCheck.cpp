#include <iostream>
#include <stdexcept>
#include <fmt/format.h>

#include "ConsistencyCheck.h"

bool check_finite(
	shared_darray<rr_float2> r,
	shared_darray<rr_int> itype,
	shared_darray<rr_float2> v,
	shared_darray<rr_float> rho,
	shared_darray<rr_float> p,
	rr_uint consistency_treatment)
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

		if (consistency_treatment == CONSISTENCY_STOP) {
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
	shared_darray<rr_int> itype,
	rr_uint consistency_treatment)
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

		if (consistency_treatment == CONSISTENCY_STOP) {
			throw std::runtime_error{ message };
		}
		else {
			printlog(message)();
			std::cerr << fmt::format("encounter {} particles outside boundaries!", count_outside_boundaries) << std::endl;
		}
	}
	return count_outside_boundaries == 0;
}


bool check_particles_have_same_position(
	shared_darray<rr_float2> r,
	shared_darray<rr_int> itype,
	rr_uint consistency_treatment)
{
	printlog_debug(__func__)();
	assert(r.get());

	rr_int count_same_position = 0;
	rr_iter any_same_position_k1 = 0;
	rr_iter any_same_position_k2 = 0;
	
	constexpr rr_float min_length_coef = 0.1;
	rr_float min_length_sqr = sqr(min_length_coef * params.hsml);

#pragma omp parallel for reduction(+: count_same_position)
	for (rr_iter k = 0; k < params.ntotal; ++k) {
		for (rr_iter other = k + 1; other < params.ntotal; ++other) {
			if (itype->at(k) == params.TYPE_NON_EXISTENT ||
				itype->at(other) == params.TYPE_NON_EXISTENT)
			{
				continue;
			}

			if (length_sqr(r->at(k) - r->at(other)) < min_length_sqr) {
				++count_same_position;

#pragma omp critical
				{
					any_same_position_k1 = k;
					any_same_position_k2 = other;
				}
			}
		}
	}


	if (count_same_position > 0) {
		std::string message = fmt::format(
			"encounter particles with same position: \n"
			"\t r1: ({}; {})\n"
			"\t r2: ({}; {})\n"
			"\t type1: {}\n"
			"\t type2: {}\n"
			"\t k1: {}\n"
			"\t k2: {}\n"
			"total pairs of particles with same positions: {}\n",
			r->at(any_same_position_k1).x,
			r->at(any_same_position_k1).y,
			r->at(any_same_position_k2).x,
			r->at(any_same_position_k2).y,
			itype->at(any_same_position_k1),
			itype->at(any_same_position_k2),
			any_same_position_k1,
			any_same_position_k2,
			count_same_position);

		if (consistency_treatment == CONSISTENCY_STOP) {
			throw std::runtime_error{ message };
		}
		else {
			printlog(message)();
		}
	}
	return count_same_position == 0;
}