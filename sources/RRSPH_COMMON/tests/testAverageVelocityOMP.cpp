#include "testArraysCommon.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "GridFind.h"
#include "SmoothingKernel.h"
#include "AverageVelocity.h"

class TestAverageVelocityOMP : public ::testing::Test {};

namespace {
	TestExperimentDirectory test_directory{ TestExperimentDirectory::Case::DamBreak };

	heap_darray<rr_int> itype;
	vheap_darray_floatn r_var;
	vheap_darray_floatn v_var;
	heap_darray<rr_float> rho;
	heap_darray<rr_float> p;

	heap_darray_md<rr_uint> neighbours;

	constexpr const char* AV_VEL_FILENAME = "average_velocity_target.csv";

	void init() {
		// prepare
		fileInput(r_var, v_var, rho, p, itype, test_directory.initial_dump_path, test_directory.experiment_dir);

		params.average_velocity = true;
		params.average_velocity_skf = 1;
		params.average_velocity_coef = 0.05;

		neighbours = heap_darray_md<rr_uint>(params.max_neighbours, params.maxn);
		grid_find(
			r_var,
			itype,
			neighbours);
	}

	auto calc_average_velocity() {
		init();

		heap_darray<rr_float2> av_vel(params.maxn);
		const auto& r = r_var.get_flt2();
		const auto& v = v_var.get_flt2();

		for_neighbour_particles(neighbours, 
			[&](rr_uint j, rr_uint i) {
				if (j >= params.nfluid) return;
				rr_float dist = distance(r(j), r(i));

				av_vel(j) += average_velocity_part(
					dist,
					itype(i),
					v(j),
					v(i),
					rho(j),
					rho(i)
				);

			});
		return av_vel;
	}
}



TEST_F(TestAverageVelocityOMP, test_average_velocity_omp)
{
	auto av_vel = calc_average_velocity();
	scalar_array_check<2>(AV_VEL_FILENAME, "av", av_vel);
}

#ifdef GEN_TEST_DATA
TEST_F(TestAverageVelocityOMP, prepare_average_velocity_omp)
{
	auto av_vel = calc_average_velocity();

	scalar_array_prepare(
		AV_VEL_FILENAME,
		"av",
		av_vel);
}
#endif