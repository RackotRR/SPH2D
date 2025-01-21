#include "testGridFind.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "CLAdapter.h"

#include "ParamsIO.h"

class TestGridFindCL : public ::testing::Test {};

namespace {
	TestExperimentDirectory test_directory{ TestExperimentDirectory::Case::DamBreak };

	heap_darray<rr_int> itype;
	vheap_darray_floatn r_var;
	vheap_darray_floatn v_var;
	heap_darray<rr_float> rho;
	heap_darray<rr_float> p;

	heap_darray_md<rr_uint> neighbours;

	cl::Buffer r_;
	cl::Buffer itype_;
	cl::Buffer neighbours_;

	void init() {
		fileInput(r_var, v_var, rho, p, itype, test_directory.initial_dump_path, test_directory.experiment_dir);

		neighbours = heap_darray_md<rr_uint>(params.max_neighbours, params.maxn);
		r_ = makeBufferCopyHost(r_var.get_flt2());
		itype_ = makeBufferCopyHost(itype);
		neighbours_ = makeBufferCopyHost(neighbours);

		// rewrite params
		params_make_header(std::filesystem::current_path() / "cl" / "clparams.h");
	}
}

TEST_F(TestGridFindCL, grid_find_cl)
{
	init();

	clProgramAdapter grid_find_adapter{ makeProgram("GridFind.cl"), cl_grid_find};
	grid_find_adapter(
		r_, itype_,
		neighbours_);

	cl::copy(neighbours_, std::begin(neighbours), std::end(neighbours));
	testGridFindCommon::check_neighbours(neighbours);
}