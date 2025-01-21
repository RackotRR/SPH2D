#include "testGridFind.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "GridFind.h"

class TestGridFindOMP : public ::testing::Test {};

namespace {
	TestExperimentDirectory test_directory{ TestExperimentDirectory::Case::DamBreak };

	heap_darray<rr_int> itype;
	vheap_darray_floatn r_var;
	vheap_darray_floatn v_var;
	heap_darray<rr_float> rho;
	heap_darray<rr_float> p;

	void init() {
		fileInput(r_var, v_var, rho, p, itype, test_directory.initial_dump_path, test_directory.experiment_dir);
	}
}

TEST_F(TestGridFindOMP, grid_find_omp)
{
	init();

	heap_darray_md<rr_uint> neighbours(params.max_neighbours, params.maxn);
	grid_find(
		r_var,
		itype,
		neighbours);

	testGridFindCommon::check_neighbours(neighbours);
}

#if 0
TEST_F(TestGridFindOMP, prepare_grid_find_omp) {
	init();

	heap_darray_md<rr_uint> neighbours(params.max_neighbours, params.maxn);
	grid_find(
		r_var,
		itype,
		neighbours);

	testGridFindCommon::neighbours_array_prepare(neighbours);
}
#endif