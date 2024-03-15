#ifdef RRSPH_OMP
#include "RRSPHOMPVersion.h"
#include "TimeIntegration.h" 
#elif defined RRSPH_CL
#include "RRSPHCLVersion.h"
#include "CLCommon.h"
#include "CLAdapter.h"
#else
static_assert(false, "undefined RRSPH Simulator");
#endif

#include <iostream>
#include <string>
#include <thread>
#include <RR/Time/Timer.h>

#include "CommonIncl.h"
#include "Input.h"
#include "TimeFormat.h"

#ifdef RRSPH_OMP
rr_uint RRSPH_GetSpecificVersionMajor() {
	return RRSPH_OMP_VERSION_MAJOR;
}
rr_uint RRSPH_GetSpecificVersionMinor() {
	return RRSPH_OMP_VERSION_MINOR;
}
rr_uint RRSPH_GetSpecificVersionPatch() {
	return RRSPH_OMP_VERSION_PATCH;
}
std::string RRSPH_GetSpecificVersionName() {
	return "RRSPH_OMP";
}
#else
rr_uint RRSPH_GetSpecificVersionMajor() {
	return RRSPH_CL_VERSION_MAJOR;
}
rr_uint RRSPH_GetSpecificVersionMinor() {
	return RRSPH_CL_VERSION_MINOR;
}
rr_uint RRSPH_GetSpecificVersionPatch() {
	return RRSPH_CL_VERSION_PATCH;
}
std::string RRSPH_GetSpecificVersionName() {
	return "RRSPH_CL";
}
#endif

void simulation() {
	rr_uint ntotal; // number of particles 
	rr_uint nfluid; 
	heap_darray<rr_int> itype; // material type of particles
	heap_darray<rr_float2> r; // coordinates of all particles
	heap_darray<rr_float2> v; // velocities of all particles
	heap_darray<rr_float> rho; // density
	heap_darray<rr_float> p; // pressure

	cli(r, v, rho, p, itype, ntotal, nfluid);

#ifdef RRSPH_CL
	logCLInfo();
#endif

	RR::Timer timer;
	timer.start();
#ifdef SPH2D_OMP
	time_integration(r, v, rho, p, itype, ntotal, nfluid);
#elif defined RRSPH_CL
	cl_time_integration(r, v, rho, p, itype, ntotal, nfluid);
#endif
	timer.finish();

	auto total_time = format_timer(timer.value<std::chrono::nanoseconds>());
	std::cout << "total time: " << total_time << std::endl;
	printlog("total time passed: ")(total_time)();
}

int main(int arc, const char* argv[]) {
	try {
		simulation();
	}
	catch (const std::exception& ex) {
		printlog("catch exception: ")(ex.what())();
	}
	return 0;
}