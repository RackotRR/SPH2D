#include "testArtificialViscosity.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "GridFind.h"
#include "Kernel.h"
#include "ArtificialViscosity.h"

class TestArtificialViscosityOMP : public ::testing::Test {};

namespace {
	TestExperimentDirectory test_directory{ TestExperimentDirectory::Case::DamBreak };

	rr_uint ntotal;
	rr_uint nfluid;
	heap_darray<rr_int> itype;
	heap_darray<rr_float2> r;
	heap_darray<rr_float2> v;
	heap_darray<rr_float> rho;
	heap_darray<rr_float> p;

	heap_darray_md<rr_uint> neighbours;
	heap_darray_md<rr_float2> dwdr;

	void init() {
		// prepare
		fileInput(r, v, rho, p, itype, ntotal, nfluid, test_directory.initial_dump_path, test_directory.experiment_dir);

		params.dt_correction_method = DT_CORRECTION_DYNAMIC;
		params.artificial_viscosity = true;

		neighbours = heap_darray_md<rr_uint>(params.max_neighbours, params.maxn);
		grid_find(ntotal,
			r,
			itype,
			neighbours);

		dwdr = heap_darray_md<rr_float2>(params.max_neighbours, params.maxn);
		calculate_kernels_dwdr(ntotal,
			r, neighbours, dwdr, params.artificial_viscosity_skf);
	}
}

static void test_artificial_viscosity(std::string filename) {
	init();

	heap_darray<rr_float> art_visc_mu(params.maxn);
	heap_darray<rr_float2> art_visc_dvdt(params.maxn);
	artificial_viscosity(ntotal,
		r, v, rho,
		neighbours, dwdr,
		art_visc_dvdt, art_visc_mu);

	testArtificialViscosityCommon::check_artificial_viscosity(
		filename,
		art_visc_dvdt,
		art_visc_mu
	);
}

TEST_F(TestArtificialViscosityOMP, test_artificial_viscosity_omp) {
	test_artificial_viscosity(testArtificialViscosityCommon::ART_VISC_FILENAME);
}

#if 0
static void prepare_artificial_viscosity(std::string filename) {
	init();

	heap_darray<rr_float> art_visc_mu(params.maxn);
	heap_darray<rr_float2> art_visc_dvdt(params.maxn);
	artificial_viscosity(ntotal,
		r, v, rho,
		neighbours, dwdr,
		art_visc_dvdt, art_visc_mu);

	testArtificialViscosityCommon::prepare_artificial_viscosity(
		filename,
		art_visc_dvdt,
		art_visc_mu
	);
}
TEST_F(TestArtificialViscosityOMP, prepare_artificial_viscosity_omp) {
	prepare_artificial_viscosity(testArtificialViscosityCommon::ART_VISC_FILENAME);
}
#endif