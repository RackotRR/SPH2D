#include "testInternalForce.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "CLAdapter.h"

#include "ParamsIO.h"

class TestInternalForceCL : public ::testing::Test {};

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
	cl::Buffer p_;
	cl::Buffer v_;
	cl::Buffer itype_;
	cl::Buffer neighbours_;
	
	struct TestConfig {
		rr_uint intf_sph_approximation;
		bool artificial_pressure;
	};

	void init(TestConfig config) {
		fileInput(r_var, v_var, rho, p, itype, test_directory.initial_dump_path, test_directory.experiment_dir);

		// rewrite params with config
		params.intf_sph_approximation = config.intf_sph_approximation;
		params.artificial_pressure = config.artificial_pressure;
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

		p_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE, params.maxn);
		cl::Buffer temp_drho_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE, params.maxn);

		clProgramAdapter con_density_adapter{ makeProgram("Density.cl"), cl_con_density };
		con_density_adapter(
			r_,
			v_,
			neighbours_,
			rho_,
			temp_drho_,
			p_
		);
	}
}

static void test_intf(TestConfig config, std::string filename) {
	init(config);

	// init test array
	heap_darray<rr_float2> test_indvdt(params.maxn);
	auto test_indvdt_ = makeBufferCopyHost(test_indvdt);
	
	// run test
	clProgramAdapter intf_adapter{ makeProgram("InternalForce.cl"), cl_internal_force };
	intf_adapter(
		r_, v_, rho_, p_,
		neighbours_,
		test_indvdt_);

	cl::copy(test_indvdt_, std::begin(test_indvdt), std::end(test_indvdt));

	// check result
	testInternalForceCommon::check_internal_force(filename, test_indvdt);
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