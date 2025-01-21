#include "testArraysCommon.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "CLAdapter.h"

#include "ParamsIO.h"

class TestAverageVelocityCL : public ::testing::Test {};

namespace {
	TestExperimentDirectory test_directory{ TestExperimentDirectory::Case::DamBreak };

	heap_darray<rr_int> itype;
	vheap_darray_floatn r_var;
	vheap_darray_floatn v_var;
	heap_darray<rr_float> rho;
	heap_darray<rr_float> p;

	heap_darray_md<rr_uint> neighbours;

	cl::Buffer r_;
	cl::Buffer rho_;
	cl::Buffer v_;
	cl::Buffer itype_;
	cl::Buffer neighbours_;

	void init() {
		fileInput(r_var, v_var, rho, p, itype, test_directory.initial_dump_path, test_directory.experiment_dir);

		// rewrite params
		params.average_velocity = true;
		params.average_velocity_skf = 1;
		params.average_velocity_coef = 0.05;
		params_make_header(std::filesystem::current_path() / "cl" / "clparams.h");

		neighbours = heap_darray_md<rr_uint>(params.max_neighbours, params.maxn);
		r_ = makeBufferCopyHost(r_var.get_flt2());
		v_ = makeBufferCopyHost(v_var.get_flt2());
		rho_ = makeBufferCopyHost(rho);
		itype_ = makeBufferCopyHost(itype);
		neighbours_ = makeBufferCopyHost(neighbours);

		clProgramAdapter grid_find_adapter{ makeProgram("GridFind.cl"), cl_grid_find };
		grid_find_adapter(
			r_, itype_,
			neighbours_);
	}
}

TEST_F(TestAverageVelocityCL, test_average_velocity_cl)
{
	init();

	// init test array
	heap_darray<rr_float2> av_vel(params.maxn);
	cl::Buffer av_vel_ = makeBufferCopyHost(av_vel);

	// run test
	auto average_velocity_adapter = clProgramAdapter{ makeProgram("AverageVelocity.cl"), cl_average_velocity };
	average_velocity_adapter(
		r_, itype_, v_, rho_,
		neighbours_,
		av_vel_
	);

	// check result
	cl::copy(av_vel_, std::begin(av_vel), std::end(av_vel));
	scalar_array_check<2>("average_velocity_target.csv", "av", av_vel);
}