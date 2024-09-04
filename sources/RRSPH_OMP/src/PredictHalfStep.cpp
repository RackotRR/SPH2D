#include "PredictHalfStep.h"
#include "Density.h"

template<typename rr_floatn>
void predict_half_step(
	const heap_darray<rr_int>& itype, // material type 
	const heap_darray<rr_float>& rho, // density
	const heap_darray<rr_float>& drho,	// density change
	const heap_darray<rr_floatn>& v,	// velocities
	const heap_darray<rr_floatn>& a,	// acceleration
	heap_darray<rr_float>& rho_predict, // half step for density
	heap_darray<rr_floatn>& v_predict)	// half step for velocities
{
	printlog()(__func__)();

	for (rr_uint i = 0; i < params.ntotal; i++) {
		if (density_is_using_continuity()) {
			rho_predict(i) = rho(i) + drho(i) * params.dt * 0.5f;
		}

		if (itype(i) > 0) {
			v_predict(i) = v(i) + a(i) * params.dt * 0.5f;
		}
	}
}

void predict_half_step(
	const heap_darray<rr_int>& itype, // material type 
	const heap_darray<rr_float>& rho, // density
	const heap_darray<rr_float>& drho,	// density change
	const vheap_darray_floatn& v_var,	// velocities of all particles
	const vheap_darray_floatn& a_var,	// acceleration
	heap_darray<rr_float>& rho_predict, // half step for density
	vheap_darray_floatn& v_predict_var) // half step for velocities
{
	if (params.dim == 2) {
		const auto& v = v_var.get_flt2();
		const auto& a = a_var.get_flt2();
		auto& v_predict = v_predict_var.get_flt2();
		predict_half_step(itype, rho, drho, v, a, rho_predict, v_predict);
	}
	else if (params.dim == 3) {
		const auto& v = v_var.get_flt3();
		const auto& a = a_var.get_flt3();
		auto& v_predict = v_predict_var.get_flt3();
		predict_half_step(itype, rho, drho, v, a, rho_predict, v_predict);
	}
	else {
		assert(0);
	}
}