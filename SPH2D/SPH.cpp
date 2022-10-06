#include <iostream>
#include <string>

#include "CommonIncl.h"
#include "Input.h"
#include "Output.h"
#include "TimeIntegration.h" 

double sqr(double value) { 
	return value * value;
}

double pow(double value, int power) {
	double result{ 1 };
	for (int i{ power }; i > 0; i--) {
		result *= value;
	}
	for (int i{ power }; i < 0; i++) {
		result /= value;
	}
	return result;
}  

void sph() {
	heap_array<double, Params::maxn> mass; // particle masses
	size_t ntotal; // number of particles
	size_t nfluid;
	heap_array<int, Params::maxn> itype;// material type of particles
	heap_array_md<double, Params::dim, Params::maxn> x;	// coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn> vx;// velocities of all particles
	heap_array<double, Params::maxn> rho; // density
	heap_array<double, Params::maxn> p;	// pressure
	heap_array<double, Params::maxn> u;	// specific internal energy
	heap_array<double, Params::maxn> c;	// sound velocity 
	heap_array<double, Params::maxn> e;	// total energy of particles

	input(x, vx, mass, rho, p, u, itype, ntotal, nfluid);
	setupOutput();

	time_integration(x, vx, mass, rho, p, u, c, e, itype, ntotal, nfluid);
}
 

int main(int arc, const char* argv[]) {
	try {
		sph();
	}
	catch (std::exception& ex) {
		std::cerr << ex.what() << std::endl;
	}
	system("pause");
	return 0;
}