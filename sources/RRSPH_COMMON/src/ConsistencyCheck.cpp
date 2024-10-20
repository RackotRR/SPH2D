#include <iostream>
#include <stdexcept>
#include <fmt/format.h>
#include <optional>

#include "ConsistencyCheck.h"

template<typename rr_floatn>
std::string format_floatn(rr_floatn v) {
	if constexpr (is_using_float3<rr_floatn>()) {
		return fmt::format("({}; {}; {})", v.x, v.y, v.z);
	}
	else {
		return fmt::format("({}; {})", v.x, v.y);
	}
}

template<typename rr_floatn>
void consistency_check_error(
	rr_uint consistency_treatment,
	std::string what,
	rr_int err_particles_count,
	rr_floatn r,
	rr_int itype,
	rr_iter k,
	std::optional<rr_floatn> v,
	std::optional<rr_float> rho,
	std::optional<rr_float> p)
{
	std::string message = what + "\n";
	message += fmt::format("\t type: {}\n", itype);

	std::string m_r = fmt::format("\t r: {}", format_floatn(r));
	std::string m_bounds;
	std::string m_xbounds = fmt::format("{} .. {}", params.x_mingeom, params.x_maxgeom);
	std::string m_ybounds = fmt::format("{} .. {}", params.y_mingeom, params.y_maxgeom);
	if constexpr (is_using_float3<rr_floatn>()) {
		std::string m_zbounds = fmt::format("{} .. {}", params.z_mingeom, params.z_maxgeom);
		m_bounds = fmt::format("({}; {}; {})\n", m_xbounds, m_ybounds, m_zbounds);
	}
	else {
		m_bounds = fmt::format("({}; {})\n", m_xbounds, m_ybounds);
	}
	message += m_r + " // " + m_bounds;

	if (v.has_value()) {
		message += fmt::format("\t v: {}\n", format_floatn(v.value()));
	}

	if (rho.has_value()) {
		message += fmt::format("\t rho: {}\n", rho.value());
	}

	if (p.has_value()) {
		message += fmt::format("\t p: {}\n", p.value());
	}

	message += fmt::format("\t k: {}\n", k);
	message += fmt::format("\t total_error_particles: {}", err_particles_count);

	printlog(message)();
	if (consistency_treatment == CONSISTENCY_STOP) {
		throw std::runtime_error{ message };
	}
	else {
		std::cerr << what << std::endl;
	}
}

static bool check_finite(rr_float2 val) {
	return std::isfinite(val.x) && std::isfinite(val.y);
}
static bool check_finite(rr_float3 val) {
	return check_finite(rr_float2{ val.x, val.y }) && std::isfinite(val.z);
}

template<typename rr_floatn>
bool check_finite(
	heap_darray<rr_floatn>& r,
	shared_darray<rr_int> itype,
	heap_darray<rr_floatn>& v,
	shared_darray<rr_float> rho,
	shared_darray<rr_float> p,
	rr_uint consistency_treatment)
{
	printlog_debug(__func__)();

	rr_int infinite_count = 0;
	rr_uint i = 0;

#pragma omp parallel for reduction(+: infinite_count)
	for (rr_iter k = 0; k < params.ntotal; ++k) {
		if (!check_finite(r.at(k)) || !check_finite(v.at(k))) {
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
		consistency_check_error<rr_floatn>(
			consistency_treatment,
			"Encounter not finite value",
			infinite_count,
			r.at(i),
			itype->at(i),
			i,
			v.at(i),
			rho->at(i),
			p->at(i));
	}

	return infinite_count == 0;
}

bool check_finite(
	shared_vheap_darray_floatn r_var,
	shared_darray<rr_int> itype,
	shared_vheap_darray_floatn v_var,
	shared_darray<rr_float> rho,
	shared_darray<rr_float> p,
	rr_uint consistency_treatment)
{
	assert(r_var.get());
	assert(itype.get());
	assert(v_var.get());
	assert(rho.get());
	assert(p.get());

	if (params.dim == 3) {
		auto& r = r_var->get_flt3();
		auto& v = v_var->get_flt3();
		return check_finite(r, itype, v, rho, p, consistency_treatment);
	}
	else {
		auto& r = r_var->get_flt2();
		auto& v = v_var->get_flt2();
		return check_finite(r, itype, v, rho, p, consistency_treatment);
	}
}

static bool check_particle_is_within_boundaries(rr_float2 r) {
	return r.x < params.x_maxgeom &&
		r.x > params.x_mingeom &&
		r.y < params.y_maxgeom &&
		r.y > params.y_mingeom &&
		std::isfinite(r.x) &&
		std::isfinite(r.y);
}
static bool check_particle_is_within_boundaries(rr_float3 r) {
	return check_particle_is_within_boundaries(rr_float2{ r.x, r.y }) &&
		r.z < params.z_maxgeom &&
		r.z > params.z_mingeom &&
		std::isfinite(r.z);
}

template<typename rr_floatn>
bool check_particles_are_within_boundaries(
	const heap_darray<rr_floatn>& r,
	const heap_darray<rr_int>& itype,
	rr_uint consistency_treatment)
{
	printlog_debug(__func__)();

	rr_int count_outside_boundaries = 0;
	rr_iter outside_k = 0;

#pragma omp parallel for reduction(+: count_outside_boundaries)
	for (rr_iter k = 0; k < params.ntotal; ++k) 
	{
		if (!check_particle_is_within_boundaries(r(k))) {
			count_outside_boundaries++;
			outside_k = k;
		}
	}

	if (count_outside_boundaries > 0) {
		consistency_check_error<rr_floatn>(
			consistency_treatment,
			"Encounter particle outside boundaries",
			count_outside_boundaries,
			r(outside_k),
			itype(outside_k),
			outside_k,
			std::nullopt,
			std::nullopt,
			std::nullopt);
	}
	return count_outside_boundaries == 0;
}

bool check_particles_are_within_boundaries(
	const vheap_darray_floatn& r_var,
	const heap_darray<rr_int>& itype,
	rr_uint consistency_treatment)
{
	if (params.dim == 3) {
		auto& r = r_var.get_flt3();
		return check_particles_are_within_boundaries(r, itype, consistency_treatment);
	}
	else {
		auto& r = r_var.get_flt2();
		return check_particles_are_within_boundaries(r, itype, consistency_treatment);
	}
}

bool check_particles_are_within_boundaries(
	shared_vheap_darray_floatn r_var,
	shared_darray<rr_int> itype,
	rr_uint consistency_treatment)
{
	return check_particles_are_within_boundaries(*r_var, *itype, consistency_treatment);
}

template<typename rr_floatn>
bool check_particles_have_same_position(
	const heap_darray<rr_floatn>& r,
	const heap_darray<rr_int>& itype,
	rr_uint consistency_treatment)
{
	printlog_debug(__func__)();

	rr_int count_same_position = 0;
	rr_iter any_same_position_k1 = 0;
	rr_iter any_same_position_k2 = 0;

	constexpr rr_float min_length_coef = 0.1;
	rr_float min_length_sqr = sqr(min_length_coef * params.hsml);

#pragma omp parallel for reduction(+: count_same_position)
	for (rr_iter k = 0; k < params.ntotal; ++k) {
		for (rr_iter other = k + 1; other < params.ntotal; ++other) {
			if (itype.at(k) == params.TYPE_NON_EXISTENT ||
				itype.at(other) == params.TYPE_NON_EXISTENT)
			{
				continue;
			}

			if (length_sqr(r.at(k) - r.at(other)) < min_length_sqr) {
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
			"\t r1: {}\n"
			"\t r2: {}\n"
			"\t type1: {}\n"
			"\t type2: {}\n"
			"\t k1: {}\n"
			"\t k2: {}\n"
			"total pairs of particles with same positions: {}\n",
			format_floatn(r.at(any_same_position_k1)),
			format_floatn(r.at(any_same_position_k2)),
			itype.at(any_same_position_k1),
			itype.at(any_same_position_k2),
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

bool check_particles_have_same_position(
	shared_vheap_darray_floatn r_var,
	shared_darray<rr_int> itype,
	rr_uint consistency_treatment)
{
	if (params.dim == 3) {
		auto& r = r_var->get_flt3();
		return check_particles_have_same_position(r, *itype, consistency_treatment);
	}
	else {
		auto& r = r_var->get_flt2();
		return check_particles_have_same_position(r, *itype, consistency_treatment);
	}
}