#include "testDensity.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "GridFind.h"
#include "SmoothingKernel.h"
#include "Density.h"

class TestDensityOMP : public ::testing::Test {};

namespace {
	TestExperimentDirectory test_directory{ TestExperimentDirectory::Case::DamBreak };

	heap_darray<rr_int> itype;
	vheap_darray_floatn r_var;
	vheap_darray_floatn v_var;
	heap_darray<rr_float> rho;
	heap_darray<rr_float> p;

	heap_darray_md<rr_uint> neighbours;

	struct TestConfig {
		rr_uint density_treatment;
		std::optional<rr_float> delta_sph_coef;
		std::optional<bool> normalization;
	};

	void init(TestConfig config) {
		// prepare
		fileInput(r_var, v_var, rho, p, itype, test_directory.initial_dump_path, test_directory.experiment_dir);

		params.density_treatment = config.density_treatment;
		params.density_delta_sph_coef = config.delta_sph_coef.value_or(0);
		params.density_normalization = config.normalization.value_or(0);

		neighbours = heap_darray_md<rr_uint>(params.max_neighbours, params.maxn);
		grid_find(
			r_var,
			itype,
			neighbours);
	}
}

TEST_F(TestDensityOMP, con_density_omp)
{
	init({ DENSITY_CONTINUITY });

	heap_darray<rr_float> test_drho(params.maxn);
	heap_darray<rr_float> test_p(params.maxn);

	density_con(
		r_var, v_var,
		neighbours,
		rho,
		test_drho, test_p);

	testDensityCommon::check_density(
		testDensityCommon::CON_DENSITY_FILENAME,
		test_drho,
		test_p,
		"drho"
	);
}

TEST_F(TestDensityOMP, sum_density_omp)
{
	init({ DENSITY_SUMMATION });

	heap_darray<rr_float> test_rho(params.maxn);
	heap_darray<rr_float> test_p(params.maxn);

	density_sum(
		r_var,
		neighbours,
		test_rho,
		test_p);

	testDensityCommon::check_density(
		testDensityCommon::SUM_DENSITY_FILENAME,
		test_rho,
		p,
		"rho"
	);
}


#if 0
//TEST_F(TestDensityOMP, prepare_con_delta_density)
//{
//	init({ DENSITY_CONTINUITY_DELTA, 0.1 });
//
//	heap_darray<rr_float> test_drho(params.maxn);
//	heap_darray<rr_float> test_p(params.maxn);
//
//	density_con(
//		r_var, v_var,
//		neighbours,
//		rho,
//		test_drho, test_p);
//
//	scalar_array_prepare(
//		"con_delta_density_target.csv",
//		"drho",
//		test_drho);
//}

TEST_F(TestDensityOMP, prepare_con_density)
{
	init({ DENSITY_CONTINUITY });

	heap_darray<rr_float> test_drho(params.maxn);
	heap_darray<rr_float> test_p(params.maxn);

	density_con(
		r_var, v_var,
		neighbours,
		rho,
		test_drho, test_p);

	testDensityCommon::prepare_density(
		testDensityCommon::CON_DENSITY_FILENAME,
		test_drho,
		test_p,
		"drho"
	);
}

TEST_F(TestDensityOMP, prepare_sum_density)
{
	init({ DENSITY_SUMMATION });

	heap_darray<rr_float> test_rho(params.maxn);
	heap_darray<rr_float> test_p(params.maxn);

	density_sum(
		r_var,
		neighbours,
		test_rho,
		test_p);

	testDensityCommon::prepare_density(
		testDensityCommon::SUM_DENSITY_FILENAME,
		test_rho,
		test_p,
		"rho"
	);
}
#endif