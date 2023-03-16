#include "CLAdapter.h"
#include <iostream>
#include <string>
#include "RR/Time/Timer.h"

#include "CommonIncl.h"
#include "Input.h"
#include "Output.h"
#include "TimeIntegration.h" 
#include "Test.h"

void sph() {
	rr_uint ntotal; // number of particles
	rr_uint nfluid;
	heap_array<rr_float, Params::maxn> mass; // particle masses
	heap_array<rr_int, Params::maxn> itype;// material type of particles
	heap_array<rr_float2, Params::maxn> r;	// coordinates of all particles
	heap_array<rr_float2, Params::maxn> v;// velocities of all particles
	heap_array<rr_float, Params::maxn> rho; // density
	heap_array<rr_float, Params::maxn> p;	// pressure
	heap_array<rr_float, Params::maxn> u;	// specific internal energy
	heap_array<rr_float, Params::maxn> c;	// sound velocity 

	input(r, v, mass, rho, p, u, itype, ntotal, nfluid);
	//time_integration(r, v, mass, rho, p, u, c, itype, ntotal, nfluid);
	cl_time_integration(r, v, mass, rho, p, u, c, itype, ntotal, nfluid);
}

void testing() {
	try {
		Params::experimentName = "test";
		setupOutput();
		initConsts();
		Test{};
	}
	catch (const std::exception& ex) {
		std::cerr << "error at " << ex.what() << std::endl;
	}
}

void simulation() {
	std::cout << "Experiment name: ";
	std::getline(std::cin, Params::experimentName);
	setupOutput();

	printlog("Experiment: ")(Params::experimentName)();

	try {
		RR::Timer timer;

		timer.start();
		sph();
		timer.finish();

		std::cout << "total time in minutes: " << timer.value<std::chrono::minutes>() << std::endl;
	}
	catch (const std::exception& ex) {
		printlog("catch exception: ")(ex.what())();
	}

	system("pause");
}

int main(int arc, const char* argv[]) {
	//testing();
	simulation();
	return 0;
}