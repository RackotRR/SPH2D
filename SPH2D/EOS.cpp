#include "CommonIncl.h"
 

// artificial equation of state for the artificial compressibility
void p_art_water(
	const double rho, // density
	const double u,	// internal energy
	double& p, // pressure
	double& c) // sound velocity
{
	// artificial EOS, Form 1 (Monaghan, 1994)
	 //static constexpr double gamma{ 7.0 };
	 //static constexpr double rho0{ 1000.0 };
	 //static constexpr double b{0.1};//{ 1.013e5 };
	 //p = b * (pow(rho / rho0, gamma) - 1.0); 
	 //
	 //static constexpr double c0{0.01};//{ 1480.0 };
	 //c = c0;

	// Lennard-Jones EOS
	constexpr double rho0{ 1000 };
	double H{ Params::d };
	double g{ Params::g };
	double cSqr{ 200 * g * H };
	double alpha{ 3 };
	double beta{ 2 };
	double K{ (rho0 * cSqr) / (alpha + beta) };

	p = K;
	if (rho >= 2 * rho0) {
		p *= pow((rho - rho0) / rho0, alpha) - pow((rho - rho0) / rho0, beta);
	}
	else if (rho >= rho0) {
		p *= pow(rho / rho0, alpha) - pow(rho / rho0, beta);
	}
	else {
		p = 0;
	}

	// Artificial EOS, Form2 (Morris, 1997)
	//c = 0.05; 
	//p = sqr(c) * rho;
}