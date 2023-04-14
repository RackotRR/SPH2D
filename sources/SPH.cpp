#ifdef SPH2D_OMP
#include "TimeIntegration.h" 
#else
#include "CLCommon.h"
#include "CLAdapter.h"
#endif // !SPH2D_OMP

#include <iostream>
#include <string>
#include <thread>
#include <RR/Time/Timer.h>

#include "CommonIncl.h"
#include "Input.h"

void simulation() {
	rr_uint ntotal; // number of particles 
	rr_uint nfluid; 
	heap_darray<rr_float> mass(0); // particle masses
	heap_darray<rr_int> itype(0); // material type of particles
	heap_darray<rr_float2> r(0); // coordinates of all particles
	heap_darray<rr_float2> v(0); // velocities of all particles
	heap_darray<rr_float> rho(0); // density
	heap_darray<rr_float> p(0); // pressure

	cli(r, v, mass, rho, p, itype, ntotal, nfluid);

#ifndef SPH2D_OMP
	logCLInfo();
#endif // !SPH2D_OMP

	printlog("Experiment: ")(params.experiment_name)();

	RR::Timer timer;
	timer.start();
#ifdef SPH2D_OMP
	time_integration(r, v, mass, rho, p, itype, ntotal, nfluid);
#else
	cl_time_integration(r, v, mass, rho, p, itype, ntotal, nfluid);
#endif // SPH2D_OMP
	timer.finish();

	std::cout << "total time in minutes: " << timer.value<std::chrono::minutes>() << std::endl;
}

int main(int arc, const char* argv[]) {
	try {
		simulation();
	}
	catch (const std::exception& ex) {
		printlog("catch exception: ")(ex.what())();
	}

#ifdef _WIN32
	system("pause");
#else // all detached output threads have to finish their work
	std::this_thread::sleep_for(std::chrono::minutes{ 1 });
#endif // _WIN32
	return 0;
}