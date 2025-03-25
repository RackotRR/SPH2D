#include "testArtificialViscosity.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "GridFind.h"
#include "SmoothingKernel.h"
#include "ArtificialViscosity.h"

class TestArtificialViscosityOMP : public ::testing::Test {};

namespace {
	TestExperimentDirectory test_directory{ TestExperimentDirectory::Case::DamBreak };

	heap_darray<rr_int> itype;
	vheap_darray_floatn r_var;
	vheap_darray_floatn v_var;
	heap_darray<rr_float> rho;
	heap_darray<rr_float> p;

	heap_darray_md<rr_uint> neighbours;

	void init() {
		// prepare
		fileInput(r_var, v_var, rho, p, itype, test_directory.initial_dump_path, test_directory.experiment_dir);

		params.dt_correction_method = DT_CORRECTION_DYNAMIC;
		params.artificial_viscosity = true;

		neighbours = heap_darray_md<rr_uint>(params.max_neighbours, params.maxn);
		grid_find(
			r_var,
			itype,
			neighbours);
	}

	auto calc_artificial_viscosity() {
		init();

		const auto& r = r_var.get_flt2();
		const auto& v = v_var.get_flt2();

		heap_darray<rr_float> art_visc_mu(params.maxn);
		heap_darray<rr_float2> art_visc_dvdt(params.maxn);

		for_neighbour_particles(neighbours, 
			[&](rr_uint j, rr_uint i) {
				art_visc_dvdt(j) += artificial_viscosity_part(
					r(j),
					r(i),
					v(j),
					v(i),
					rho(j),
					rho(i),
					&art_visc_mu(j)
				);
			});
		return std::make_pair(std::move(art_visc_dvdt), std::move(art_visc_mu));
	}
}


TEST_F(TestArtificialViscosityOMP, test_artificial_viscosity_omp) {
	auto [art_visc_dvdt, art_visc_mu] = calc_artificial_viscosity();
	testArtificialViscosityCommon::check_artificial_viscosity(
		testArtificialViscosityCommon::ART_VISC_FILENAME,
		art_visc_dvdt,
		art_visc_mu
	);
}

#ifdef GEN_TEST_DATA
TEST_F(TestArtificialViscosityOMP, prepare_artificial_viscosity_omp) {
	auto [art_visc_dvdt, art_visc_mu] = calc_artificial_viscosity();
	testArtificialViscosityCommon::prepare_artificial_viscosity(
		testArtificialViscosityCommon::ART_VISC_FILENAME,
		art_visc_dvdt,
		art_visc_mu
	);
}
#endif