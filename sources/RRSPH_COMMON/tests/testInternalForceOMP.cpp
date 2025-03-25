#include "testInternalForce.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "GridFind.h"
#include "SmoothingKernel.h"
#include "InternalForce.h"

#include "Density.h"

class TestInternalForceOMP : public ::testing::Test {};

namespace {
	TestExperimentDirectory test_directory{ TestExperimentDirectory::Case::DamBreak };

	heap_darray<rr_int> itype;
	vheap_darray_floatn r_var;
	vheap_darray_floatn v_var;
	heap_darray<rr_float> rho;
	heap_darray<rr_float> p;

	heap_darray_md<rr_uint> neighbours;

	struct TestConfig {
		rr_uint intf_sph_approximation;
		bool artificial_pressure;
	};

	void init(TestConfig config) {
		// prepare
		fileInput(r_var, v_var, rho, p, itype, test_directory.initial_dump_path, test_directory.experiment_dir);

		params.intf_sph_approximation = config.intf_sph_approximation;
		params.artificial_pressure = config.artificial_pressure;

		neighbours = heap_darray_md<rr_uint>(params.max_neighbours, params.maxn);
		grid_find(
			r_var,
			itype,
			neighbours);
	}

	auto calc_intf(TestConfig config) {
		init(config);

		const auto& r = r_var.get_flt2();
		const auto& v = v_var.get_flt2();
		heap_darray<rr_float2> indvdt(params.maxn);

		heap_darray<rr_float> temp_drho(params.maxn);
		heap_darray<rr_float> test_p(params.maxn);

		params.density_treatment = DENSITY_CONTINUITY;
		density_con(
			r,
			v,
			neighbours,
			rho,
			temp_drho,
			test_p
		);

		for_neighbour_particles(neighbours,
			[&](rr_uint j, rr_uint i) {
				rr_float2 diff_ij = r(i) - r(j);
				rr_float dist_ij = length(diff_ij);
				indvdt(j) += find_internal_changes_part(
					diff_ij,
					dist_ij,
					test_p(j),
					test_p(i),
					rho(j),
					rho(i)
				);
			}
		);

		return indvdt;
	}

	void test_intf(TestConfig config, std::string filename) {
		auto test_indvdt = calc_intf(config);
		testInternalForceCommon::check_internal_force(filename, test_indvdt);
	}
}

TEST_F(TestInternalForceOMP, intf1_no_art_pressure_omp) {
	test_intf({ INTF_SPH_APPROXIMATION_1, false }, testInternalForceCommon::INTF1_NO_ARTP_FILENAME);
}
TEST_F(TestInternalForceOMP, intf2_no_art_pressure_omp) {
	test_intf({ INTF_SPH_APPROXIMATION_2, false }, testInternalForceCommon::INTF2_NO_ARTP_FILENAME);
}
TEST_F(TestInternalForceOMP, intf1_art_pressure_omp) {
	test_intf({ INTF_SPH_APPROXIMATION_1, true }, testInternalForceCommon::INTF1_ARTP_FILENAME);
}
TEST_F(TestInternalForceOMP, intf2_art_pressure_omp) {
	test_intf({ INTF_SPH_APPROXIMATION_2, true }, testInternalForceCommon::INTF2_ARTP_FILENAME);
}

#ifdef GEN_TEST_DATA
static void prepare_intf(TestConfig config, std::string filename) {
	auto test_indvdt = calc_intf(config);
	testInternalForceCommon::prepare_internal_force(filename, test_indvdt);
}

TEST_F(TestInternalForceOMP, prepare_intf1_no_art_pressure) {
	prepare_intf({ INTF_SPH_APPROXIMATION_1, false }, testInternalForceCommon::INTF1_NO_ARTP_FILENAME);
}
TEST_F(TestInternalForceOMP, prepare_intf2_no_art_pressure) {
	prepare_intf({ INTF_SPH_APPROXIMATION_2, false }, testInternalForceCommon::INTF2_NO_ARTP_FILENAME);
}
TEST_F(TestInternalForceOMP, prepare_intf1_art_pressure) {
	prepare_intf({ INTF_SPH_APPROXIMATION_1, true }, testInternalForceCommon::INTF1_ARTP_FILENAME);
}
TEST_F(TestInternalForceOMP, prepare_intf2_art_pressure) {
	prepare_intf({ INTF_SPH_APPROXIMATION_2, true }, testInternalForceCommon::INTF2_ARTP_FILENAME);
}
#endif