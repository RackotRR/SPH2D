#include "testGridFind.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "GridFind.h"

class TestGridFindOMP : public ::testing::Test {};

namespace {
	TestExperimentDirectory test_directory{ TestExperimentDirectory::Case::DamBreak };

	rr_uint ntotal;
	rr_uint nfluid;
	heap_darray<rr_int> itype;
	heap_darray<rr_float2> r;
	heap_darray<rr_float2> v;
	heap_darray<rr_float> rho;
	heap_darray<rr_float> p;

	void init() {
		fileInput(r, v, rho, p, itype, ntotal, nfluid, test_directory.initial_dump_path, test_directory.experiment_dir);
	}
}

TEST_F(TestGridFindOMP, grid_find_omp)
{
	init();

	heap_darray_md<rr_uint> neighbours(params.max_neighbours, params.maxn);
	grid_find(ntotal,
		r,
		itype,
		neighbours);

	testGridFindCommon::check_neighbours(neighbours);
}

static void prepare_test_data() {
	init();

	heap_darray_md<rr_uint> neighbours(params.max_neighbours, params.maxn);
	grid_find(ntotal,
		r,
		itype,
		neighbours);

	testGridFindCommon::neighbours_array_prepare(neighbours);
}