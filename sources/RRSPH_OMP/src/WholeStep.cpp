#include "WholeStep.h"
#include "Density.h"

inline bool is_point_within_geometry(const rr_float2& point) {
	return point.x < params.x_maxgeom &&
		point.x > params.x_mingeom &&
		point.y < params.y_maxgeom &&
		point.y > params.y_mingeom;
}
inline bool is_point_within_geometry(const rr_float3& point) {
	return point.x < params.x_maxgeom &&
		point.x > params.x_mingeom &&
		point.y < params.y_maxgeom &&
		point.y > params.y_mingeom &&
		point.z < params.z_maxgeom &&
		point.z > params.z_mingeom;
}

template<typename rr_floatn>
void whole_step(
	const rr_uint ntotal,
	const rr_uint timestep,
	const heap_darray<rr_float>& drho,	// density change
	const heap_darray<rr_floatn>& a,	// acceleration
	const heap_darray<rr_floatn>& av,	// average velocity
	heap_darray<rr_int>& itype, // material type 
	heap_darray<rr_float>& rho, // density
	heap_darray<rr_floatn>& v,	// velocities
	heap_darray<rr_floatn>& r)	// coordinates of all particles
{
	printlog()(__func__)();

	rr_float r_dt = params.dt;
	rr_float v_dt = timestep ? params.dt : params.dt * 0.5f;

	for (rr_uint i = 0; i < ntotal; i++) {
		if (density_is_using_continuity()) {
			rho(i) = rho(i) + drho(i) * v_dt;
		}

		if (itype(i) > 0) {
			v(i) += a(i) * v_dt;

			if (params.average_velocity) {
				v(i) += av(i);
			}

			if (params.consistency_treatment == CONSISTENCY_FIX) {
				rr_float2 new_r = r(i) + v(i) * r_dt;
				if (is_point_within_geometry(new_r)) {
					r(i) = new_r;
				}
				else {
					itype(i) = params.TYPE_NON_EXISTENT;
				}
			}
			else {
				r(i) += v(i) * r_dt;
			}
		}
	}
}

void whole_step(
	const rr_uint ntotal,
	const rr_uint timestep,
	const heap_darray<rr_float>& drho,	// density change
	const vheap_darray_floatn& a_var,	// acceleration
	const vheap_darray_floatn& av_var,	// average velocity
	heap_darray<rr_int>& itype, // material type 
	heap_darray<rr_float>& rho, // density
	vheap_darray_floatn& v_var,	// velocities
	vheap_darray_floatn& r_var)	// coordinates of all particles
{

	if (params.dim == 2) {
		auto& r = r_var.get_flt2();
		auto& v = v_var.get_flt2();
		const auto& a = a_var.get_flt2();
		const auto& av = av_var.get_flt2();
		whole_step(ntotal, timestep, drho, a, av, itype, rho, v, r);
	}
	else if (params.dim == 3) {
		auto& r = r_var.get_flt3();
		auto& v = v_var.get_flt3();
		const auto& a = a_var.get_flt3();
		const auto& av = av_var.get_flt3();
		whole_step(ntotal, timestep, drho, a, av, itype, rho, v, r);
	}
	else {
		assert(0);
	}
}