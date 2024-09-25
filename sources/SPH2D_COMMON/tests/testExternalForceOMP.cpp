#include "testExternalForce.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "GridFind.h"
#include "Kernel.h"
#include "ExtForce.h"

class TestExternalForceOMP : public ::testing::Test {};

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
	heap_darray_md<rr_float> w;
	heap_darray_md<rr_float2> dwdr;

	void init() {
		// prepare
		fileInput(r, v, rho, p, itype, ntotal, nfluid, test_directory.initial_dump_path, test_directory.experiment_dir);

		neighbours = heap_darray_md<rr_uint>(params.max_neighbours, params.maxn);
		grid_find(ntotal,
			r,
			itype,
			neighbours);

		w = heap_darray_md<rr_float>(params.max_neighbours, params.maxn);
		calculate_kernels_w(ntotal,
			r, neighbours, w, params.density_skf);

		dwdr = heap_darray_md<rr_float2>(params.max_neighbours, params.maxn);
		calculate_kernels_dwdr(ntotal,
			r, neighbours, dwdr, params.density_skf);
	}
}

static void check_external_forces(std::string filename, rr_uint boundary_treatment) {
	init();

	params.boundary_treatment = boundary_treatment;
	heap_darray<rr_float2> test_exdvdt(params.maxn);
	external_force(ntotal,
		r,
		neighbours, itype,
		test_exdvdt);

	testExternalForceCommon::check_external_forces(filename, test_exdvdt);
}

TEST_F(TestExternalForceOMP, test_sbt_dynamic_omp) {
	check_external_forces(testExternalForceCommon::FILENAME_DYNAMIC, SBT_DYNAMIC);
}
TEST_F(TestExternalForceOMP, test_sbt_repulsive_omp) {
	check_external_forces(testExternalForceCommon::FILENAME_REPULSIVE, SBT_REPULSIVE);
}


#if 0
static void prepare_external_forces(std::string filename, rr_uint boundary_treatment) {
	init();

	params.boundary_treatment = boundary_treatment;
	heap_darray<rr_float2> test_exdvdt(params.maxn);
	external_force(ntotal,
		r,
		neighbours, itype,
		test_exdvdt);

	testExternalForceCommon::prepare_external_forces(filename, test_exdvdt);
}


TEST_F(TestExternalForceOMP, prepare_sbt_dynamic) {
	prepare_external_forces(testExternalForceCommon::FILENAME_DYNAMIC, SBT_DYNAMIC);
}
TEST_F(TestExternalForceOMP, prepare_sbt_repulsive) {
	prepare_external_forces(testExternalForceCommon::FILENAME_REPULSIVE, SBT_REPULSIVE);
}
#endif