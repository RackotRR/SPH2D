#include <iostream>
#include <string>
#include <RRTime/Timer.h>

#include "CommonIncl.h"
#include "Input.h"
#include "Output.h"
#include "TimeIntegration.h" 

rr_float sqr(rr_float value) { 
	return value * value;
}
rr_float pow(rr_float value, rr_int power) {
	rr_float result{ 1 };
	for (rr_int i{ power }; i > 0; i--) {
		result *= value;
	}
	for (rr_int i{ power }; i < 0; i++) {
		result /= value;
	}
	return result;
}
rr_float pow(rr_float value, rr_uint power) {
	rr_float result{ 1 };
	for (rr_uint i{ power }; i > 0; i--) {
		result *= value;
	}
	return result;
}

void sph() {
	heap_array<rr_float, Params::maxn> mass; // particle masses
	rr_uint ntotal; // number of particles
	rr_uint nfluid;
	heap_array<rr_int, Params::maxn> itype;// material type of particles
	heap_array<rr_float2, Params::maxn> r;	// coordinates of all particles
	heap_array<rr_float2, Params::maxn> v;// velocities of all particles
	heap_array<rr_float, Params::maxn> rho; // density
	heap_array<rr_float, Params::maxn> p;	// pressure
	heap_array<rr_float, Params::maxn> u;	// specific internal energy
	heap_array<rr_float, Params::maxn> c;	// sound velocity 
	heap_array<rr_float, Params::maxn> e;	// total energy of particles

	input(r, v, mass, rho, p, u, itype, ntotal, nfluid);
	setupOutput();

	time_integration(r, v, mass, rho, p, u, c, e, itype, ntotal, nfluid);
}
 

int main(int arc, const char* argv[]) {
	try {
		RR::Timer timer;
		timer.start();
		sph();
		timer.finish();
		std::cout << "total time in seconds: " << timer.value<std::chrono::seconds>() << std::endl;
	}
	catch (std::exception& ex) {
		std::cerr << ex.what() << std::endl;
	}
	system("pause");
	return 0;
}