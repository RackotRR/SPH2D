#include "testDensity.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "CLAdapter.h"

#include "ParamsIO.h"

class TestDensityCL : public ::testing::Test {};

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

	void init(rr_uint density_treatment) {
		fileInput(r_var, v_var, rho, p, itype, test_directory.initial_dump_path, test_directory.experiment_dir);

		// rewrite params with target density_treatment
		params.density_treatment = density_treatment;
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

TEST_F(TestDensityCL, con_density_cl)
{
	init(DENSITY_CONTINUITY);

	// init test array
	heap_darray<rr_float> test_drho(params.maxn);
	heap_darray<rr_float> test_p(params.maxn);
	cl::Buffer test_drho_ = makeBufferCopyHost(test_drho);
	cl::Buffer test_p_ = makeBufferCopyHost(p);

	// run test
	clProgramAdapter con_density_adapter{ makeProgram("Density.cl"), cl_con_density };
	con_density_adapter(
		r_,
		v_,
		neighbours_,
		rho_,
		test_drho_,
		test_p_
	);

	// check result
	cl::copy(test_drho_, std::begin(test_drho), std::end(test_drho));
	cl::copy(test_p_, std::begin(test_p), std::end(test_p));

	testDensityCommon::check_density(
		testDensityCommon::CON_DENSITY_FILENAME,
		test_drho,
		test_p,
		"drho"
	);
}

TEST_F(TestDensityCL, sum_density_cl)
{
	init(DENSITY_SUMMATION);

	// init test array
	heap_darray<rr_float> test_rho(params.maxn);
	heap_darray<rr_float> test_p(params.maxn);
	cl::Buffer test_rho_ = makeBufferCopyHost(test_rho);
	cl::Buffer test_p_ = makeBufferCopyHost(p);

	// run test
	clProgramAdapter sum_density_adapter{ makeProgram("Density.cl"), cl_sum_density };
	sum_density_adapter(
		r_,
		neighbours_,
		test_rho_,
		test_p_
	);

	// check result
	cl::copy(test_rho_, std::begin(test_rho), std::end(test_rho));
	cl::copy(test_p_, std::begin(test_p), std::end(test_p));

	testDensityCommon::check_density(
		testDensityCommon::SUM_DENSITY_FILENAME,
		test_rho,
		test_p,
		"rho"
	);
}