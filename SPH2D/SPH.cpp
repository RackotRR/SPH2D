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
	
	double simulationTime = 2;
	double dt = 5e-4;
	if (simulationTime / dt < 0) {
		std::cerr << "maxtimestep error" << std::endl;
		return;
	}

	Params::maxtimestep = static_cast<size_t>(simulationTime / dt);

	input(x, vx, mass, rho, p, u, itype, ntotal, nfluid);

	try {
		time_integration(x, vx, mass, rho, p, u, c, e, itype, ntotal, nfluid, Params::maxtimestep, dt);
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