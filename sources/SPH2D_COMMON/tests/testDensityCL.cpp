#include "testArraysCommon.h"
#include "testExperimentDirectoryCase.h"

#include "Input.h"
#include "CLAdapter.h"

#include "ParamsIO.h"

class TestDensityCL : public ::testing::Test {};

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

	void init(rr_uint density_treatment) {
		fileInput(r, v, rho, p, itype, ntotal, nfluid, test_directory.initial_dump_path, test_directory.experiment_dir);

		// rewrite params with target density_treatment
		params.density_treatment = density_treatment;
		params_make_header(std::filesystem::current_path() / "cl" / "clparams.h");

		heap_darray_md<rr_uint> neighbours(params.max_neighbours, params.maxn);
		r_ = makeBufferCopyHost(r);
		v_ = makeBufferCopyHost(v);
		rho_ = makeBufferCopyHost(rho);
		itype_ = makeBufferCopyHost(itype);
		neighbours_ = makeBufferCopyHost(neighbours);

		cl::Buffer w_ = makeBuffer<rr_float>(CL_MEM_READ_WRITE, params.max_neighbours * params.maxn);
		cl::Buffer dwdr_ = makeBuffer<rr_float2>(CL_MEM_READ_WRITE, params.max_neighbours * params.maxn);
		smoothing_kernels_w[params.density_skf] = std::move(w_);
		smoothing_kernels_dwdr[params.density_skf] = std::move(dwdr_);

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

TEST_F(TestDensityCL, con_density_cl)
{
	init(DENSITY_CONTINUITY);

	// init test array
	heap_darray<rr_float> test_drho(params.maxn);
	cl::Buffer test_drho_ = makeBufferCopyHost(test_drho);

	// run test
	clProgramAdapter con_density_adapter{ makeProgram("Density.cl"), cl_con_density };
	con_density_adapter(
		r_,
		v_,
		neighbours_, smoothing_kernels_dwdr[params.density_skf],
		rho_,
		test_drho_
	);

	// check result
	cl::copy(test_drho_, std::begin(test_drho), std::end(test_drho));
	scalar_array_check<2>(
		"con_density_target.csv",
		"drho",
		test_drho);
}

TEST_F(TestDensityCL, sum_density_cl)
{
	init(DENSITY_SUMMATION);

	// init test array
	heap_darray<rr_float> test_rho(params.maxn);
	cl::Buffer test_rho_ = makeBufferCopyHost(test_rho);

	// run test
	clProgramAdapter sum_density_adapter{ makeProgram("Density.cl"), cl_sum_density };
	sum_density_adapter(
		neighbours_, smoothing_kernels_w[params.density_skf],
		test_rho_
	);

	// check result
	cl::copy(test_rho_, std::begin(test_rho), std::end(test_rho));
	scalar_array_check<2>(
		"sum_density_target.csv",
		"rho",
		test_rho);
}