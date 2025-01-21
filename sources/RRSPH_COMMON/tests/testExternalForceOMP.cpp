#include "testExternalForce.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "GridFind.h"
#include "SmoothingKernel.h"
#include "ExtForce.h"

class TestExternalForceOMP : public ::testing::Test {};

namespace {
	TestExperimentDirectory test_directory{ TestExperimentDirectory::Case::DamBreak };

	heap_darray<rr_int> itype;
	vheap_darray_floatn r_var;
	vheap_darray_floatn v_var;
	heap_darray<rr_float> rho;
	heap_darray<rr_float> p;

	heap_darray_md<rr_uint> neighbours;

	void init() {
		// prepare
		fileInput(r_var, v_var, rho, p, itype, test_directory.initial_dump_path, test_directory.experiment_dir);

		neighbours = heap_darray_md<rr_uint>(params.max_neighbours, params.maxn);
		grid_find(
			r_var,
			itype,
			neighbours);
	}

	auto calc_external_force(rr_uint boundary_treatment) {
		init();

		params.boundary_treatment = boundary_treatment;
		heap_darray<rr_float2> test_exdvdt(params.maxn, rr_float2{ 0.f, -params.g });

		const auto& r = r_var.get_flt2();
		const auto& v = v_var.get_flt2();

		for_neighbour_particles(neighbours,
			[&](rr_uint j, rr_uint i) {
				rr_float dist = distance(r(j), r(i));
				test_exdvdt(j) += external_force_part(
					dist,
					r(j),
					r(i),
					itype(j),
					itype(i)
				);
			});

		return test_exdvdt;
	}
}

static void check_external_forces(std::string filename, rr_uint boundary_treatment) {
	auto test_exdvdt = calc_external_force(boundary_treatment);
	testExternalForceCommon::check_external_forces(filename, test_exdvdt);
}

TEST_F(TestExternalForceOMP, test_sbt_dynamic_omp) {
	check_external_forces(testExternalForceCommon::FILENAME_DYNAMIC, SBT_DYNAMIC);
}
//TODO: test doesn't pass. Error is small, option is rarely used
TEST_F(TestExternalForceOMP, test_sbt_repulsive_omp) {
	check_external_forces(testExternalForceCommon::FILENAME_REPULSIVE, SBT_REPULSIVE);
}


#if 0
static void prepare_external_forces(std::string filename, rr_uint boundary_treatment) {
	auto test_exdvdt = calc_external_force(boundary_treatment);
	testExternalForceCommon::prepare_external_forces(filename, test_exdvdt);
}


TEST_F(TestExternalForceOMP, prepare_sbt_dynamic) {
	prepare_external_forces(testExternalForceCommon::FILENAME_DYNAMIC, SBT_DYNAMIC);
}
TEST_F(TestExternalForceOMP, prepare_sbt_repulsive) {
	prepare_external_forces(testExternalForceCommon::FILENAME_REPULSIVE, SBT_REPULSIVE);
}
#endif