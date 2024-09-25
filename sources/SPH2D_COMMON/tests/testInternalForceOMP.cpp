#include "testInternalForce.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "GridFind.h"
#include "Kernel.h"
#include "InternalForce.h"

class TestInternalForceOMP : public ::testing::Test {};

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

	struct TestConfig {
		rr_uint intf_sph_approximation;
		bool artificial_pressure;
	};

	void init(TestConfig config) {
		// prepare
		fileInput(r, v, rho, p, itype, ntotal, nfluid, test_directory.initial_dump_path, test_directory.experiment_dir);

		params.intf_sph_approximation = config.intf_sph_approximation;
		params.artificial_pressure = config.artificial_pressure;

		neighbours = heap_darray_md<rr_uint>(params.max_neighbours, params.maxn);
		grid_find(ntotal,
			r,
			itype,
			neighbours);

		if (params.artificial_pressure) {
			w = heap_darray_md<rr_float>(params.max_neighbours, params.maxn);
			calculate_kernels_w(ntotal,
				r, neighbours, w, params.artificial_pressure_skf);
		}
		else {
			w = heap_darray_md<rr_float>{};
		}

		dwdr = heap_darray_md<rr_float2>(params.max_neighbours, params.maxn);
		calculate_kernels_dwdr(ntotal,
			r, neighbours, dwdr, params.intf_skf);
	}
}

static void test_intf(TestConfig config, std::string filename) {
	init(config);

	heap_darray<rr_float> test_p(params.maxn);
	heap_darray<rr_float2> test_indvdt(params.maxn);
	int_force(ntotal,
		r, v, rho,
		neighbours,
		w,
		dwdr,
		test_p, test_indvdt);

	testInternalForceCommon::check_internal_force(filename, test_indvdt, test_p);
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

#if 0
static void prepare_intf(TestConfig config, std::string filename) {
	init(config);

	heap_darray<rr_float> test_p(params.maxn);
	heap_darray<rr_float2> test_indvdt(params.maxn);
	int_force(ntotal,
		r, v, rho,
		neighbours,
		w,
		dwdr,
		test_p, test_indvdt);

	testInternalForceCommon::prepare_internal_force(filename, test_indvdt, test_p);
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