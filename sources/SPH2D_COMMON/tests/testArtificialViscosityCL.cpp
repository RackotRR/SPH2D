#include "testArtificialViscosity.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "CLAdapter.h"

#include "ParamsIO.h"

class TestArtificialViscosityCL : public ::testing::Test {};

namespace {
	TestExperimentDirectory test_directory{ TestExperimentDirectory::Case::DamBreak };

	rr_uint ntotal;
	rr_uint nfluid;
	heap_darray<rr_int> itype;
	heap_darray<rr_float2> r;
	heap_darray<rr_float2> v;
	heap_darray<rr_float> rho;
	heap_darray<rr_float> p;

	cl::Buffer r_;
	cl::Buffer v_;
	cl::Buffer rho_;
	cl::Buffer itype_;
	cl::Buffer neighbours_;

	std::unordered_map<rr_uint, cl::Buffer> smoothing_kernels_w;
	std::unordered_map<rr_uint, cl::Buffer> smoothing_kernels_dwdr;


	void init() {
		fileInput(r, v, rho, p, itype, ntotal, nfluid, test_directory.initial_dump_path, test_directory.experiment_dir);

		// rewrite params with config
		params.dt_correction_method = DT_CORRECTION_DYNAMIC;
		params.artificial_viscosity = true;
		params_make_header(std::filesystem::current_path() / "cl" / "clparams.h");

		heap_darray_md<rr_uint> neighbours(params.max_neighbours, params.maxn);
		r_ = makeBufferCopyHost(r);
		v_ = makeBufferCopyHost(v);
		rho_ = makeBufferCopyHost(rho);
		itype_ = makeBufferCopyHost(itype);
		neighbours_ = makeBufferCopyHost(neighbours);

		cl::Buffer dwdr_ = makeBuffer<rr_float2>(CL_MEM_READ_WRITE, params.max_neighbours * params.maxn);
		smoothing_kernels_dwdr[params.artificial_viscosity_skf] = std::move(dwdr_);

		clProgramAdapter grid_find_adapter{ makeProgram("GridFind.cl"), cl_grid_find };
		clProgramAdapter skf_adapter{ makeProgram("SmoothingKernel.cl"), cl_calculate_kernels };

		// calc grid and skf

		grid_find_adapter(
			r_, itype_,
			neighbours_);

		skf_adapter(
			r_, itype_, neighbours_,
			smoothing_kernels_w, smoothing_kernels_dwdr);
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
		neighbours_, smoothing_kernels_dwdr[params.artificial_viscosity_skf],
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