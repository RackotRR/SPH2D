#include "CLCommon.h"
#include "Test.h"
#include "TimeIntegration.h"
#include "Input.h"
#include "GridFind.h"
#include "Density.h"
#include "ExtForce.h"
#include "InternalForce.h"
#include "ArtificialViscosity.h"
#include "SingleStep.h"

void update_change_rate_gpu(rr_uint nfluid,
	const heap_array<rr_float2, Params::maxn>& indvxdt_cl,
	const heap_array<rr_float2, Params::maxn>& exdvxdt_cl,
	const heap_array<rr_float2, Params::maxn>& arvdvxdt_cl,
	const heap_array<rr_float, Params::maxn>& arvdudt_cl,
	const heap_array<rr_float, Params::maxn>& ahdudt_cl,
	heap_array<rr_float2, Params::maxn>& a_cl,
	heap_array<rr_float, Params::maxn>& dudt_cl)
{
	printlog_debug(__func__)();

	static RRKernel kernel(makeProgram("TimeIntegration.cl"), "single_step");

	auto indvxdt_ = makeBufferCopyHost(CL_MEM_READ_ONLY, indvxdt_cl);
	auto exdvxdt_ = makeBufferCopyHost(CL_MEM_READ_ONLY, exdvxdt_cl);
	auto arvdvxdt_ = makeBufferCopyHost(CL_MEM_READ_ONLY, arvdvxdt_cl);
	auto arvdudt_ = makeBufferCopyHost(CL_MEM_READ_ONLY, arvdudt_cl);
	auto a_ = makeBufferCopyHost(CL_MEM_WRITE_ONLY, a_cl);
	auto dudt_ = makeBufferCopyHost(CL_MEM_WRITE_ONLY, dudt_cl);

	kernel(
		dudt_,
		arvdudt_,
		indvxdt_,
		exdvxdt_,
		arvdvxdt_,
		dudt_,
		a_
	).execute(Params::maxn, Params::localThreads);

	cl::copy(a_, a_cl.begin(), a_cl.end());
	cl::copy(dudt_, dudt_cl.begin(), dudt_cl.end());
}

