#include "testArraysCommon.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "CLAdapter.h"

#include "ParamsIO.h"

class TestAverageVelocityCL : public ::testing::Test {};

namespace {
	TestExperimentDirectory test_directory{ TestExperimentDirectory::Case::DamBreak };

	rr_uint ntotal;
	rr_uint nfluid;
	heap_darray<rr_int> itype;
	heap_darray<rr_float2> r;
	heap_darray<rr_float2> v;
	heap_darray<rr_float> rho;
	heap_darray<rr_float> p;

	cl::Buffer r_;
	cl::Buffer v_;
	cl::Buffer rho_;
	cl::Buffer itype_;
	cl::Buffer neighbours_;

	std::unordered_map<rr_uint, cl::Buffer> smoothing_kernels_w;
	std::unordered_map<rr_uint, cl::Buffer> smoothing_kernels_dwdr;

	void init() {
		fileInput(r, v, rho, p, itype, ntotal, nfluid, test_directory.initial_dump_path, test_directory.experiment_dir);

		// rewrite params
		params.average_velocity = true;
		params.average_velocity_skf = 1;
		params.average_velocity_coef = 0.05;
		params_make_header(std::filesystem::current_path() / "cl" / "clparams.h");

		heap_darray_md<rr_uint> neighbours(params.max_neighbours, params.maxn);
		r_ = makeBufferCopyHost(r);
		v_ = makeBufferCopyHost(v);
		rho_ = makeBufferCopyHost(rho);
		itype_ = makeBufferCopyHost(itype);
		neighbours_ = makeBufferCopyHost(neighbours);

		cl::Buffer w_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE, params.max_neighbours * params.maxn);
		smoothing_kernels_w[params.average_velocity_skf] = std::move(w_);

		clProgramAdapter grid_find_adapter{ makeProgram("GridFind.cl"), cl_grid_find };
		clProgramAdapter skf_adapter{ makeProgram("SmoothingKernel.cl"), cl_calculate_kernels };

		// calc grid and skf

		grid_find_adapter(
			r_, itype_,
			neighbours_);

		skf_adapter(
			r_, itype_, neighbours_,
			smoothing_kernels_w, smoothing_kernels_dwdr);
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
		neighbours_, smoothing_kernels_w[params.average_velocity_skf],
		av_vel_
	);

	// check result
	cl::copy(av_vel_, std::begin(av_vel), std::end(av_vel));
	scalar_array_check<2>("average_velocity_target.csv", "av", av_vel);
}