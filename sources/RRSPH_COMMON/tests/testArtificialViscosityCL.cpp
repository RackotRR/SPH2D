#include "testArtificialViscosity.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "CLAdapter.h"

#include "ParamsIO.h"

class TestArtificialViscosityCL : public ::testing::Test {};

namespace {
	TestExperimentDirectory test_directory{ TestExperimentDirectory::Case::DamBreak };

	heap_darray<rr_int> itype;
	vheap_darray_floatn r_var;
	vheap_darray_floatn v_var;
	heap_darray<rr_float> rho;
	heap_darray<rr_float> p;

	heap_darray_md<rr_uint> neighbours;

	cl::Buffer r_;
	cl::Buffer rho_;
	cl::Buffer v_;
	cl::Buffer itype_;
	cl::Buffer neighbours_;


	void init() {
		fileInput(r_var, v_var, rho, p, itype, test_directory.initial_dump_path, test_directory.experiment_dir);

		// rewrite params
		params.dt_correction_method = DT_CORRECTION_DYNAMIC;
		params.artificial_viscosity = true;
		params_make_header(std::filesystem::current_path() / "cl" / "clparams.h");

		neighbours = heap_darray_md<rr_uint>(params.max_neighbours, params.maxn);
		r_ = makeBufferCopyHost(r_var.get_flt2());
		v_ = makeBufferCopyHost(v_var.get_flt2());
		rho_ = makeBufferCopyHost(rho);
		itype_ = makeBufferCopyHost(itype);
		neighbours_ = makeBufferCopyHost(neighbours);

		clProgramAdapter grid_find_adapter{ makeProgram("GridFind.cl"), cl_grid_find };
		grid_find_adapter(
			r_, itype_,
			neighbours_);
	}
}

static void test_artificial_viscosity(std::string filename) {
	init();

	// init test array
	heap_darray<rr_float> art_visc_mu(params.maxn);
	heap_darray<rr_float2> art_visc_dvdt(params.maxn);
	cl::Buffer art_visc_mu_ = makeBufferCopyHost(art_visc_mu);
	cl::Buffer art_visc_dvdt_ = makeBufferCopyHost(art_visc_dvdt);

	// run test
	auto art_visc_adapter = clProgramAdapter{ makeProgram("ArtificialViscosity.cl"), cl_artificial_viscosity };
	art_visc_adapter(
		r_, v_, rho_,
		neighbours_,
		art_visc_dvdt_, art_visc_mu_);

	cl::copy(art_visc_dvdt_, std::begin(art_visc_dvdt), std::end(art_visc_dvdt));
	cl::copy(art_visc_mu_, std::begin(art_visc_mu), std::end(art_visc_mu));

	// check result
	testArtificialViscosityCommon::check_artificial_viscosity(
		filename,
		art_visc_dvdt,
		art_visc_mu
	);
}

TEST_F(TestArtificialViscosityCL, test_artificial_viscosity_cl) {
	test_artificial_viscosity(testArtificialViscosityCommon::ART_VISC_FILENAME);
}