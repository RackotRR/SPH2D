#include "testExternalForce.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "CLAdapter.h"

#include "ParamsIO.h"

class TestExternalForceCL : public ::testing::Test {};

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

	struct TestConfig {
		rr_uint boundary_treatment;
	};

	void init(TestConfig config) {
		fileInput(r, v, rho, p, itype, ntotal, nfluid, test_directory.initial_dump_path, test_directory.experiment_dir);

		// rewrite params with config
		params.boundary_treatment = config.boundary_treatment;
		params_make_header(std::filesystem::current_path() / "cl" / "clparams.h");

		heap_darray_md<rr_uint> neighbours(params.max_neighbours, params.maxn);
		r_ = makeBufferCopyHost(r);
		v_ = makeBufferCopyHost(v);
		rho_ = makeBufferCopyHost(rho);
		itype_ = makeBufferCopyHost(itype);
		neighbours_ = makeBufferCopyHost(neighbours);

		clProgramAdapter grid_find_adapter{ makeProgram("GridFind.cl"), cl_grid_find };

		// calc grid 

		grid_find_adapter(
			r_, itype_,
			neighbours_);
	}
}

static void check_external_forces(TestConfig config, std::string filename) {
	init(config);

	// init test array
	heap_darray<rr_float2> test_exdvdt(params.maxn);
	cl::Buffer test_exdvdt_ = makeBufferCopyHost(test_exdvdt);

	// run test
	clProgramAdapter external_force_adapter{ makeProgram("ExternalForce.cl"), cl_external_force };
	external_force_adapter(
		r_, neighbours_, itype_,
		test_exdvdt_);

	cl::copy(test_exdvdt_, std::begin(test_exdvdt), std::end(test_exdvdt));

	// check result
	testExternalForceCommon::check_external_forces(filename, test_exdvdt);
}

TEST_F(TestExternalForceCL, test_sbt_dynamic_cl) {
	check_external_forces({ SBT_DYNAMIC }, testExternalForceCommon::FILENAME_DYNAMIC);
}
TEST_F(TestExternalForceCL, test_sbt_repulsive_cl) {
	check_external_forces({ SBT_REPULSIVE }, testExternalForceCommon::FILENAME_REPULSIVE);
}