#include "CommonIncl.h"
 

// artificial equation of state for the artificial compressibility
void p_art_water(
    const double rho, // density
    const double u,	// internal energy
    double& p, // pressure
    double& c) // sound velocity
{

    constexpr double rho0 = 1000;
    double H = Params::d;
    double g = Params::g;
    double cSqr = 200 * g * H;
    if (Params::eos == 1) {

        // Lennard-Jones EOS
        int alpha = 3;
        int beta = 2;
        double K = (rho0 * cSqr) / (alpha + beta);

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
    }
    else {
        // artificial EOS, Form (Monaghan, 1994)
        static constexpr double gamma = 7.0;
        double B = cSqr * rho0 / gamma;
        if (rho > rho0) {
            p = B * (pow(rho / rho0, gamma) - 1.0);
        }
        else {
            p = 0;
        }
    }
}