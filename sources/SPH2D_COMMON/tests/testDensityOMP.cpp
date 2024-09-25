#include "testArraysCommon.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "GridFind.h"
#include "Kernel.h"
#include "Density.h"

class TestDensityOMP : public ::testing::Test {};

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
		rr_uint density_treatment;
		std::optional<rr_float> delta_sph_coef;
		std::optional<bool> normalization;
	};

	void init(TestConfig config) {
		// prepare
		fileInput(r, v, rho, p, itype, ntotal, nfluid, test_directory.initial_dump_path, test_directory.experiment_dir);

		params.density_treatment = config.density_treatment;
		params.density_delta_sph_coef = config.delta_sph_coef.value_or(0);
		params.density_normalization = config.normalization.value_or(0);

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

TEST_F(TestDensityOMP, con_density_omp)
{
	init({ DENSITY_CONTINUITY });

	heap_darray<rr_float> test_drho(params.maxn);
	con_density(ntotal,
		r, v,
		neighbours, dwdr,
		rho,
		test_drho);

	scalar_array_check<2>("con_density_target.csv", "drho", test_drho);
}

TEST_F(TestDensityOMP, sum_density_omp)
{
	init({ DENSITY_SUMMATION });

	heap_darray<rr_float> test_rho(params.maxn);
	sum_density(ntotal,
		neighbours, w,
		test_rho);

	scalar_array_check<2>("sum_density_target.csv", "rho", test_rho);
}


#if 0
TEST_F(TestDensityOMP, prepare_con_delta_density)
{
	init({ DENSITY_CONTINUITY_DELTA, 0.1 });

	heap_darray<rr_float> test_drho(params.maxn);
	con_density(ntotal,
		r, v,
		neighbours, dwdr,
		rho,
		test_drho);

	scalar_array_prepare(
		"con_delta_density_target.csv",
		"drho",
		test_drho);
}

TEST_F(TestDensityOMP, prepare_con_density)
{
	init({ DENSITY_CONTINUITY });

	heap_darray<rr_float> test_drho(params.maxn);
	con_density(ntotal,
		r, v,
		neighbours, dwdr,
		rho,
		test_drho);

	scalar_array_prepare(
		"con_density_target.csv",
		"drho",
		test_drho);
}

TEST_F(TestDensityOMP, prepare_sum_density)
{
	init({ DENSITY_SUMMATION });

	heap_darray<rr_float> test_rho(params.maxn);
	sum_density(ntotal,
		neighbours, w,
		test_rho);

	scalar_array_prepare(
		"sum_density_target.csv",
		"rho",
		test_rho);
}
#endif