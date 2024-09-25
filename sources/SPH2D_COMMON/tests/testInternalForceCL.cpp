#include "testInternalForce.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "CLAdapter.h"

#include "ParamsIO.h"

class TestInternalForceCL : public ::testing::Test {};

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
	cl::Buffer w_;
	cl::Buffer dwdr_;

	std::unordered_map<rr_uint, cl::Buffer> smoothing_kernels_w;
	std::unordered_map<rr_uint, cl::Buffer> smoothing_kernels_dwdr;

	struct TestConfig {
		rr_uint intf_sph_approximation;
		bool artificial_pressure;
	};

	void init(TestConfig config) {
		fileInput(r, v, rho, p, itype, ntotal, nfluid, test_directory.initial_dump_path, test_directory.experiment_dir);

		// rewrite params with config
		params.intf_sph_approximation = config.intf_sph_approximation;
		params.artificial_pressure = config.artificial_pressure;
		params_make_header(std::filesystem::current_path() / "cl" / "clparams.h");

		heap_darray_md<rr_uint> neighbours(params.max_neighbours, params.maxn);
		r_ = makeBufferCopyHost(r);
		v_ = makeBufferCopyHost(v);
		rho_ = makeBufferCopyHost(rho);
		itype_ = makeBufferCopyHost(itype);
		neighbours_ = makeBufferCopyHost(neighbours);

		cl::Buffer w_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE, params.max_neighbours * params.maxn);
		cl::Buffer dwdr_ = makeBuffer<rr_float2>(CL_MEM_READ_WRITE, params.max_neighbours * params.maxn);
		smoothing_kernels_w[params.artificial_pressure_skf] = std::move(w_);
		smoothing_kernels_dwdr[params.intf_skf] = std::move(dwdr_);

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

static void test_intf(TestConfig config, std::string filename) {
	init(config);

	// init test array
	heap_darray<rr_float> test_p(params.maxn);
	heap_darray<rr_float2> test_indvdt(params.maxn);
	cl::Buffer test_p_ = makeBufferCopyHost(test_p);
	cl::Buffer test_indvdt_ = makeBufferCopyHost(test_indvdt);

	// run test
	clProgramAdapter intf_adapter{ makeProgram("InternalForce.cl"), cl_internal_force };
	intf_adapter(
		v_, rho_, neighbours_,
		smoothing_kernels_w[params.artificial_pressure_skf],
		smoothing_kernels_dwdr[params.intf_skf],
		test_p_, test_indvdt_);

	cl::copy(test_indvdt_, std::begin(test_indvdt), std::end(test_indvdt));
	cl::copy(test_p_, std::begin(test_p), std::end(test_p));

	// check result
	testInternalForceCommon::check_internal_force(filename, test_indvdt, test_p);
}

TEST_F(TestInternalForceCL, intf1_no_art_pressure_cl) {
	test_intf({ INTF_SPH_APPROXIMATION_1, false }, testInternalForceCommon::INTF1_NO_ARTP_FILENAME);
}
TEST_F(TestInternalForceCL, intf2_no_art_pressure_cl) {
	test_intf({ INTF_SPH_APPROXIMATION_2, false }, testInternalForceCommon::INTF2_NO_ARTP_FILENAME);
}
TEST_F(TestInternalForceCL, intf1_art_pressure_cl) {
	test_intf({ INTF_SPH_APPROXIMATION_1, true }, testInternalForceCommon::INTF1_ARTP_FILENAME);
}
TEST_F(TestInternalForceCL, intf2_art_pressure_cl) {
	test_intf({ INTF_SPH_APPROXIMATION_2, true }, testInternalForceCommon::INTF2_ARTP_FILENAME);
}