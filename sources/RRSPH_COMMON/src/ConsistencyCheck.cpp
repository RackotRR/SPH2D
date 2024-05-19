#include <iostream>
#include <stdexcept>
#include <fmt/format.h>
#include <optional>

#include "ConsistencyCheck.h"

template<typename rr_floatn>
void consistency_check_error(
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

	std::string m_r;
	std::string m_bounds;
	std::string m_xbounds = fmt::format("{} .. {}", params.x_mingeom, params.x_maxgeom);
	std::string m_ybounds = fmt::format("{} .. {}", params.y_mingeom, params.y_maxgeom);
	if constexpr (is_using_float3<rr_floatn>()) {
		m_r = fmt::format("\t r: ({}; {}; {})", r.x, r.y, r.z);
		std::string m_zbounds = fmt::format("{} .. {}", params.z_mingeom, params.z_maxgeom);
		m_bounds = fmt::format("({}; {}; {})\n", m_xbounds, m_ybounds, m_zbounds);
	}
	else {
		m_r = fmt::format("\t r: ({}; {})", r.x, r.y);
		m_bounds = fmt::format("({}; {})\n", m_xbounds, m_ybounds);
	}
	message += m_r + " // " + m_bounds;

	if (v.has_value()) {
		std::string m_v;
		if constexpr (is_using_float3<rr_floatn>()) {
			auto& value_v = v.value();
			m_v = fmt::format("({}; {}; {})", value_v.x, value_v.y, value_v.z);
		}
		else {
			auto& value_v = v.value();
			m_v = fmt::format("({}; {})", value_v.x, value_v.y);
		}

		message += fmt::format("\t v: {}\n", m_v);
	}

	if (rho.has_value()) {
		message += fmt::format("\t rho: {}\n", rho.value());
	}

	if (p.has_value()) {
		message += fmt::format("\t p: {}\n", p.value());
	}

	message += fmt::format("\t k: {}\n", k);
	message += fmt::format("\t total_error_particles: {}", err_particles_count);

	if (params.consistency_treatment == CONSISTENCY_STOP) {
		throw std::runtime_error{ message };
	}
	else {
		printlog(message)();
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
	shared_darray<rr_float> p)
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
	shared_darray<rr_float> p)
{
	assert(r.get());
	assert(itype.get());
	assert(v.get());
	assert(rho.get());
	assert(p.get());

	if (params.dim == 3) {
		auto& r = r_var->get_flt3();
		auto& v = v_var->get_flt3();
		return check_finite(r, itype, v, rho, p);
	}
	else {
		auto& r = r_var->get_flt2();
		auto& v = v_var->get_flt2();
		return check_finite(r, itype, v, rho, p);
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
	const heap_darray<rr_int>& itype)
{
	printlog_debug(__func__)();
	assert(r.get());
	assert(itype.get());

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
	const heap_darray<rr_int>& itype)
{
	if (params.dim == 3) {
		auto& r = r_var.get_flt3();
		return check_particles_are_within_boundaries(r, itype);
	}
	else {
		auto& r = r_var.get_flt2();
		return check_particles_are_within_boundaries(r, itype);
	}
}

bool check_particles_are_within_boundaries(
	shared_vheap_darray_floatn r_var,
	shared_darray<rr_int> itype)
{
	return check_particles_are_within_boundaries(*r_var, *itype);
}