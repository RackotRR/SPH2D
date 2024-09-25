#include "testGridFind.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "CLAdapter.h"

class TestGridFindCL : public ::testing::Test {};

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

	cl::Buffer r_;
	cl::Buffer itype_;
	cl::Buffer neighbours_;

	void init() {
		fileInput(r, v, rho, p, itype, ntotal, nfluid, test_directory.initial_dump_path, test_directory.experiment_dir);

		neighbours = heap_darray_md<rr_uint>(params.max_neighbours, params.maxn);
		r_ = makeBufferCopyHost(r);
		itype_ = makeBufferCopyHost(itype);
		neighbours_ = makeBufferCopyHost(neighbours);
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