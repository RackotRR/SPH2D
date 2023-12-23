#ifdef SPH2D_OMP
#include "SPH2DOMPVersion.h"
#include "TimeIntegration.h" 
#elif defined SPH2D_CL
#include "SPH2DCLVersion.h"
#include "CLCommon.h"
#include "CLAdapter.h"
#else
static_assert(false, "undefined SPH2D Simulator");
#endif

#include <iostream>
#include <string>
#include <thread>
#include <RR/Time/Timer.h>

#include "CommonIncl.h"
#include "Input.h"
#include "TimeFormat.h"

#ifdef SPH2D_OMP
rr_uint SPH2D_GetSpecificVersionMajor() {
	return SPH2D_OMP_VERSION_MAJOR;
}
rr_uint SPH2D_GetSpecificVersionMinor() {
	return SPH2D_OMP_VERSION_MINOR;
}
rr_uint SPH2D_GetSpecificVersionPatch() {
	return SPH2D_OMP_VERSION_PATCH;
}
std::string SPH2D_GetSpecificVersionName() {
	return "SPH2D_OMP";
}
#else
rr_uint SPH2D_GetSpecificVersionMajor() {
	return SPH2D_CL_VERSION_MAJOR;
}
rr_uint SPH2D_GetSpecificVersionMinor() {
	return SPH2D_CL_VERSION_MINOR;
}
rr_uint SPH2D_GetSpecificVersionPatch() {
	return SPH2D_CL_VERSION_PATCH;
}
std::string SPH2D_GetSpecificVersionName() {
	return "SPH2D_CL";
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

#ifdef SPH2D_CL
	logCLInfo();
#endif

	RR::Timer timer;
	timer.start();
#ifdef SPH2D_OMP
	time_integration(r, v, rho, p, itype, ntotal, nfluid);
#elif defined SPH2D_CL
	cl_time_integration(r, v, rho, p, itype, ntotal, nfluid);
#endif
	timer.finish();

	std::cout << "total time: " << format_timer(timer.value<std::chrono::nanoseconds>()) << std::endl;
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