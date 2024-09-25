#include "testArraysCommon.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "GridFind.h"
#include "Kernel.h"
#include "AverageVelocity.h"

class TestAverageVelocityOMP : public ::testing::Test {};

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

	constexpr const char* AV_VEL_FILENAME = "average_velocity_target.csv";

	void init() {
		// prepare
		fileInput(r, v, rho, p, itype, ntotal, nfluid, test_directory.initial_dump_path, test_directory.experiment_dir);

		params.average_velocity = true;
		params.average_velocity_skf = 1;
		params.average_velocity_coef = 0.05;

		neighbours = heap_darray_md<rr_uint>(params.max_neighbours, params.maxn);
		grid_find(ntotal,
			r,
			itype,
			neighbours);

		w = heap_darray_md<rr_float>(params.max_neighbours, params.maxn);
		calculate_kernels_w(ntotal,
			r, neighbours, w, params.average_velocity_skf);
	}
}

TEST_F(TestAverageVelocityOMP, test_average_velocity_omp)
{
	init();

	heap_darray<rr_float2> av_vel(params.maxn);
	average_velocity(nfluid,
		r, itype, v, rho,
		neighbours, w,
		av_vel);

	scalar_array_check<2>(AV_VEL_FILENAME, "av", av_vel);
}

#if 0
TEST_F(TestAverageVelocityOMP, prepare_average_velocity_omp)
{
	init();

	heap_darray<rr_float2> av_vel(params.maxn);
	average_velocity(nfluid,
		r, itype, v, rho,
		neighbours, w,
		av_vel);

	scalar_array_prepare(
		AV_VEL_FILENAME,
		"av",
		av_vel);
}
#endif