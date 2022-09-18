#include "CommonIncl.h"
#include "EOS.h"

// 144 page, 4.38; 4.41; 4.59; 4.58 equations 
/*	Calculate the internal forces on the rigth hand side of the Navier-Stokes equations,
*	i.e  the pressure gradient and the gradient of the viscous stress tensor, used by the time integration.
*	Moreover the entropy production due to viscous dissipation, tds/dt,
*	and the change of internal energy per mass, de/dt, are calculated
*/
void int_force(
	const double dt, // time step
	const size_t ntotal, // number of particles,
	const heap_array<double, Params::maxn>& mass,// particle masses
	const heap_array_md<double, Params::dim, Params::maxn>& vx,	// velocities of all particles
	const size_t niac,	// number of interaction pairs
	const heap_array<double, Params::maxn>& rho,	// density
	const heap_array<double, Params::maxn>& eta,	// dynamic viscosity
	const heap_array<size_t, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<size_t, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array_md<double, Params::dim, Params::max_interaction>& dwdx,   // derivative of kernel with respect to x, y, z
	const heap_array<int, Params::maxn>& itype,	 // particle material type
	const heap_array<double, Params::maxn>& u,	// specific internal energy
	const heap_array_md<double, Params::dim, Params::maxn>& x,	// coordinates of all particles 
	heap_array<double, Params::maxn>& c,	// particle sound speed
	heap_array<double, Params::maxn>& p,	// particle pressure
	heap_array_md<double, Params::dim, Params::maxn>& dvxdt,	// acceleration with respect to x, y, z
	heap_array<double, Params::maxn>& tdsdt,	// production of viscous entropy
	heap_array<double, Params::maxn>& dedt)	// change of specific internal energy
{
	size_t i, j;
	double  hxx, hyy, hxy, h, hvcc, he, rhoij;
	heap_array<double, Params::dim> dvx;
	heap_array<double, Params::maxn> vcc;

	heap_array<double, Params::maxn> txx;
	heap_array<double, Params::maxn> tyy;
	heap_array<double, Params::maxn> txy;

	// initialization of shear tensor, velocity divergence, viscous energy, internal energy, acceleration
	for (int i = 0; i < ntotal; i++) {
		tdsdt(i) = 0;
		dedt(i) = 0;
		for (int d = 0; d < Params::dim; d++) {
			dvxdt(d, i) = 0;
		}
	}

	// calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c
	if (Params::visc) {
		for (int k = 0; k < niac; k++) {
			i = pair_i(k);
			j = pair_j(k);
			for (int d = 0; d < Params::dim; d++) {
				dvx(d) = vx(d, j) - vx(d, i);
			}

			hxx = 2.0 * dvx(0) * dwdx(0, k) - dvx(1) * dwdx(1, k);
			hxy = dvx(0) * dwdx(1, k) + dvx(1) * dwdx(0, k);
			hyy = 2.0 * dvx(1) * dwdx(1, k) - dvx(0) * dwdx(0, k);

			hxx *= 2.0 / 3.0;
			hyy *= 2.0 / 3.0;

			txx(i) += mass(j) * hxx / rho(j);
			txx(j) += mass(i) * hxx / rho(i);
			txy(i) += mass(j) * hxy / rho(j);
			txy(j) += mass(i) * hxy / rho(i);
			tyy(i) += mass(j) * hyy / rho(j);
			tyy(j) += mass(i) * hyy / rho(i);

			// calculate SPH sum for vc, c = dvx/dx + dvy/dy + dvz/dz
			hvcc = 0;
			for (int d = 0; d < Params::dim; d++) {
				hvcc += dvx(d) * dwdx(d, k);
			}
			vcc(i) += mass(j) * hvcc / rho(j);
			vcc(j) += mass(i) * hvcc / rho(i);
		}
	}

	for (int i = 0; i < ntotal; i++) {
		// viscous entropy Tds/dt = 1/2 eta/rho Tab Tab
		if (Params::visc) {
			tdsdt(i) = sqr(txx(i)) + 2 * sqr(txy(i)) + sqr(tyy(i));
			tdsdt(i) = 0.5 * eta(i) / rho(i) * tdsdt(i);
		}

		// pressure from equation of state 
		p_art_water(rho(i), u(i), p(i), c(i));
	}

	// calculate SPH sum for pressure force -p, a/rho
	// and viscous force (eta Tab), b / rho
	// and the internal energy change de/dt due to -p/rho vc, c
	for (int k = 0; k < niac; k++) {
		i = pair_i(k);
		j = pair_j(k);
		he = 0;

		if (Params::pa_sph == 1) { // for sph algorithm 1
			rhoij = 1.0 / (rho(i) * rho(j));
			for (int d = 0; d < Params::dim; d++) {
				// pressure part
				h = -(p(i) + p(j)) * dwdx(d, k);
				he += (vx(d, j) - vx(d, i)) * h;

				// viscous force
				if (Params::visc) {
					if (d == 1) { // x-coordinate of acceleration 
						h += (eta(i) * txx(i) + eta(j) * txx(j)) * dwdx(0, k);
						h += (eta(i) * txy(i) + eta(j) * txy(j)) * dwdx(1, k);
					}
					else if (d == 2) { // y-coordinate of acceleration 
						h += (eta(i) * txy(i) + eta(j) * txy(j)) * dwdx(0, k) +
							(eta(i) * tyy(i) + eta(j) * tyy(j)) * dwdx(1, k);
					}
				}
				h *= rhoij;
				dvxdt(d, i) += mass(j) * h;
				dvxdt(d, j) -= mass(i) * h;
			}
			he *= rhoij;
			dedt(i) += mass(j) * he;
			dedt(j) += mass(i) * he;
		}
		else if (Params::pa_sph == 2) { // for sph algorithm 2
			for (int d = 0; d < Params::dim; d++) {
				h = -(p(i) / sqr(rho(i)) + p(j) / sqr(rho(j))) * dwdx(d, k);
				he += (vx(d, j) - vx(d, i)) * h;

				// viscous force
				if (Params::visc) {
					if (d == 1) { // x-coordinate of acceleration
						h += (eta(i) * txx(i) / sqr(rho(i)) + eta(j) * txx(j) / sqr(rho(j))) * dwdx(0, k);
						h += (eta(i) * txy(i) / sqr(rho(i)) + eta(j) * txy(j) / sqr(rho(j))) * dwdx(1, k);
					}
					else if (d == 2) { // y-coordinate of acceleration
						h += (eta(i) * txy(i) / sqr(rho(i)) + eta(j) * txy(j) / sqr(rho(j))) * dwdx(0, k) +
							(eta(i) * tyy(i) / sqr(rho(i)) + eta(j) * tyy(j) / sqr(rho(j))) * dwdx(1, k);
					}
				}
				dvxdt(d, i) += mass(j) * h;
				dvxdt(d, j) -= mass(i) * h;
			}
			dedt(i) += mass(j) * he;
			dedt(j) += mass(i) * he;
		}
	}

	// change of specific internal energy de/dt = T ds/ds - p/rho vc, c:
	for (int i = 0; i < ntotal; i++) {
		dedt(i) = tdsdt(i) + 0.5 * dedt(i);
	}
}