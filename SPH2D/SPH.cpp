#include <iostream>

#include "CommonIncl.h"
#include "Input.h"
#include "Output.h"
#include "TimeIntegration.h" 

double sqr(double value) { 
	return value * value;
}

double pow(double value, unsigned power) {
	double result{ 1 };
	for (unsigned i{}; i < power; i++) {
		result *= value;
	}
	return result;
}  

void sph() {
	heap_array<double, Params::maxn> mass; // particle masses
	size_t ntotal; // number of particles
	double dt; // timestep
	heap_array<int, Params::maxn> itype;// material type of particles
	heap_array_md<double, Params::dim, Params::maxn> x;	// coordinates of all particles
	heap_array_md<double, Params::dim, Params::maxn> vx;// velocities of all particles
	heap_array<double, Params::maxn> rho; // density
	heap_array<double, Params::maxn> p;	// pressure
	heap_array<double, Params::maxn> u;	// specific internal energy
	heap_array<double, Params::maxn> c;	// sound velocity 
	heap_array<double, Params::maxn> e;	// total energy of particles
	
	Params::maxtimestep = 2500;
	if (Params::shocktube) {
		dt = 0.005;
	}
	if (Params::shearcavity) {
		dt = 4e-6;
	}

	input(x, vx, mass, rho, p, u, itype, ntotal);

	try {
		time_integration(x, vx, mass, rho, p, u, c, e, itype, ntotal, Params::maxtimestep, dt);
	}
	catch (std::exception& e) {
		std::cerr << e.what() << std::endl;
	}
}
 

int main() {  
	sph(); 
	system("pause");
	return 0;
}