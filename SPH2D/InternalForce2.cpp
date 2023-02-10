#include "CommonIncl.h"
#include "EOS.h"
#include "Kernel.h"

// 144 page, 4.38; 4.41; 4.59; 4.58 equations 
/*	Calculate the internal forces on the rigth hand side of the Navier-Stokes equations,
*	i.e  the pressure gradient and the gradient of the viscous stress tensor, used by the time integration.
*	Moreover the entropy production due to viscous dissipation, tds/dt,
*	and the change of internal energy per mass, de/dt, are calculated
*/
void int_force2(
	const rr_uint ntotal, // number of particles, 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& eta,	// dynamic viscosity
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	const heap_array<rr_uint, Params::maxn>& grid, // particles indices sorted so particles in the same cell are one after another
	const heap_array<rr_uint, Params::max_cells>& cell_starts_in_grid, // indices of first particle in cell
	heap_array<rr_float, Params::maxn>& c,	// particle sound speed
	heap_array<rr_float, Params::maxn>& p,	// particle pressure
	heap_array<rr_float2, Params::maxn>& a,	// acceleration with respect to x, y, z
	heap_array<rr_float, Params::maxn>& tdsdt,	// production of viscous entropy
	heap_array<rr_float, Params::maxn>& dedt)	// change of specific internal energy
{
	static heap_array<rr_float, Params::maxn> vcc;
	static heap_array<rr_float, Params::maxn> txx;
	static heap_array<rr_float, Params::maxn> tyy;
	static heap_array<rr_float, Params::maxn> txy;

	// initialization of shear tensor, velocity divergence, viscous energy, internal energy, acceleration
	tdsdt.fill(0);
	dedt.fill(0);
	vcc.fill(0);
	txx.fill(0);
	tyy.fill(0);
	txy.fill(0);
	a.fill({ 0.f });

	// calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c
	if constexpr (Params::visc) {
#pragma omp parallel for
		for (rr_iter j = 0; j < ntotal; j++) { // run through all particles
			rr_uint center_cell_idx = get_cell_idx(r(j));

			rr_uint neighbour_cells[9];
			get_neighbouring_cells(center_cell_idx, neighbour_cells);
			for (rr_uint cell_i = 0; cell_i < 9; ++cell_i) { // run through neighbouring cells
				rr_uint cell_idx = neighbour_cells[cell_i];
				if (cell_idx == Params::max_cells) continue; // invalid cell

				for (rr_uint grid_i = cell_starts_in_grid(cell_idx); // run through all particles in cell
					grid_i < cell_starts_in_grid(cell_idx + 1ull);
					++grid_i)
				{
					rr_uint i = grid(grid_i); // index of particle
					// j - current particle; i - particle near

					rr_float wij;
					rr_float2 dwdr;
					kernel(r(i), r(j), wij, dwdr);

					rr_float2 dvx = v(j) - v(i);

					rr_float hxx = 2.f * dvx.x * dwdr.x - dvx.y * dwdr.y;
					rr_float hxy = dvx.x * dwdr.y + dvx.y * dwdr.x;
					rr_float hyy = 2.f * dvx.y * dwdr.y - dvx.x * dwdr.x;

					hxx *= 2.f / 3.f;
					hyy *= 2.f / 3.f;

					txx(j) += mass(i) * hxx / rho(i);
					txy(j) += mass(i) * hxy / rho(i);
					tyy(j) += mass(i) * hyy / rho(i);

					// calculate SPH sum for vc, c = dvx/dx + dvy/dy + dvz/dz
					rr_float hvcc = dvx.x * dwdr.x + dvx.y * dwdr.y;
					vcc(j) += mass(i) * hvcc / rho(i);
				}
			}
		}
	}

	for (rr_uint i = 0; i < ntotal; i++) {
		// viscous entropy Tds/dt = 1/2 eta/rho Tab Tab
		if constexpr (Params::visc) {
			tdsdt(i) = sqr(txx(i)) + 2.f * sqr(txy(i)) + sqr(tyy(i));
			tdsdt(i) = 0.5f * eta(i) / rho(i) * tdsdt(i);
		}

		// pressure from equation of state 
		p_art_water(rho(i), u(i), p(i), c(i));
	}

	// calculate SPH sum for pressure force -p, a/rho
	// and viscous force (eta Tab), b / rho
	// and the internal energy change de/dt due to -p/rho vc, c
#pragma omp parallel for
	for (rr_iter j = 0; j < ntotal; j++) { // run through all particles
		rr_uint center_cell_idx = get_cell_idx(r(j));

		rr_uint neighbour_cells[9];
		get_neighbouring_cells(center_cell_idx, neighbour_cells);
		for (rr_uint cell_i = 0; cell_i < 9; ++cell_i) { // run through neighbouring cells
			rr_uint cell_idx = neighbour_cells[cell_i];
			if (cell_idx == Params::max_cells) continue; // invalid cell

			for (rr_uint grid_i = cell_starts_in_grid(cell_idx); // run through all particles in cell
				grid_i < cell_starts_in_grid(cell_idx + 1ull);
				++grid_i)
			{
				rr_uint i = grid(grid_i); // index of particle
				// j - current particle; i - particle near

				rr_float wij;
				rr_float2 dwdr;
				kernel(r(i), r(j), wij, dwdr);

				rr_float he = 0.f;
				if constexpr (Params::pa_sph == 1) { // for sph algorithm 1
					rr_float rhoij = 1.f / (rho(i) * rho(j));
					{
						rr_float h = -(p(i) + p(j)) * dwdr.x;
						he += (v(j).x - v(i).x) * h;

						// viscous force
						if (Params::visc) {
							h += (eta(i) * txx(i) + eta(j) * txx(j)) * dwdr.x;
							h += (eta(i) * txy(i) + eta(j) * txy(j)) * dwdr.y;
						}

						h *= rhoij;
						a(j).x -= mass(j) * h;
					}
					{
						// pressure part
						rr_float h = -(p(i) + p(j)) * dwdr.y;
						he += (v(j).y - v(i).y) * h;

						// viscous force
						if (Params::visc) {
							h += (eta(i) * txy(i) + eta(j) * txy(j)) * dwdr.x +
								(eta(i) * tyy(i) + eta(j) * tyy(j)) * dwdr.y;
						}

						h *= rhoij;
						a(j).y -= mass(j) * h;
					}

					he *= rhoij;
					dedt(j) += mass(i) * he;
				}
				else if constexpr (Params::pa_sph == 2) { // for sph algorithm 2
					{
						rr_float h = -(p(i) / sqr(rho(i)) + p(j) / sqr(rho(j))) * dwdr.x;
						he += (v(j).x - v(i).x) * h;

						// viscous force
						if constexpr (Params::visc) {
							h += (eta(i) * txx(i) / sqr(rho(i)) + eta(j) * txx(j) / sqr(rho(j))) * dwdr.x;
							h += (eta(i) * txy(i) / sqr(rho(i)) + eta(j) * txy(j) / sqr(rho(j))) * dwdr.y;
						}
						a(j).x -= mass(i) * h;
					}
					{
						rr_float h = -(p(i) / sqr(rho(i)) + p(j) / sqr(rho(j))) * dwdr.y;
						he += (v(j).y - v(i).y) * h;

						// viscous force
						if constexpr (Params::visc) {
							h += (eta(i) * txy(i) / sqr(rho(i)) + eta(j) * txy(j) / sqr(rho(j))) * dwdr.x +
								(eta(i) * tyy(i) / sqr(rho(i)) + eta(j) * tyy(j) / sqr(rho(j))) * dwdr.y;
						}
						a(j).y -= mass(i) * h;
					}
					dedt(j) += mass(i) * he;
				}
			}
		}
	}

	// change of specific internal energy de/dt = T ds/ds - p/rho vc, c:
	for (rr_uint i = 0; i < ntotal; i++) {
		dedt(i) = tdsdt(i) + 0.5f * dedt(i);
	}
}

// 144 page, 4.38; 4.41; 4.59; 4.58 equations 
/*	Calculate the internal forces on the rigth hand side of the Navier-Stokes equations,
*	i.e  the pressure gradient and the gradient of the viscous stress tensor, used by the time integration.
*	Moreover the entropy production due to viscous dissipation, tds/dt,
*	and the change of internal energy per mass, de/dt, are calculated
*/
void int_force(
	const rr_uint ntotal, // number of particles, 
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles 
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const rr_uint niac,	// number of interaction pairs
	const heap_array<rr_float, Params::maxn>& rho,	// density
	const heap_array<rr_float, Params::maxn>& eta,	// dynamic viscosity
	const heap_array<rr_uint, Params::max_interaction>& pair_i,  // list of first partner of interaction pair
	const heap_array<rr_uint, Params::max_interaction>& pair_j,  // list of second partner of interaction pair
	const heap_array<rr_float2, Params::max_interaction>& dwdx,   // derivative of kernel with respect to x, y, z
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	heap_array<rr_float, Params::maxn>& c,	// particle sound speed
	heap_array<rr_float, Params::maxn>& p,	// particle pressure
	heap_array<rr_float2, Params::maxn>& a,	// acceleration with respect to x, y, z
	heap_array<rr_float, Params::maxn>& tdsdt,	// production of viscous entropy
	heap_array<rr_float, Params::maxn>& dedt)	// change of specific internal energy
{
	rr_uint i, j;
	rr_float  hxx, hyy, hxy, h, hvcc, he, rhoij;
	rr_float2 dvx = { 0.f };

	static heap_array<rr_float, Params::maxn> vcc;
	static heap_array<rr_float, Params::maxn> txx;
	static heap_array<rr_float, Params::maxn> tyy;
	static heap_array<rr_float, Params::maxn> txy;

	// initialization of shear tensor, velocity divergence, viscous energy, internal energy, acceleration
	for (rr_uint i = 0; i < ntotal; i++) {
		tdsdt(i) = 0;
		dedt(i) = 0;
		vcc(i) = 0;
		txx(i) = 0;
		tyy(i) = 0;
		txy(i) = 0;
		a(i) = { 0.f };
	}

	// calculate SPH sum for shear tensor Tab = va,b + vb,a - 2/3 delta_ab vc,c
	if constexpr (Params::visc) {
		for (rr_uint k = 0; k < niac; k++) {
			i = pair_i(k);
			j = pair_j(k);

			dvx = v(j) - v(i);

			hxx = 2.f * dvx.x * dwdx(k).x - dvx.y * dwdx(k).y;
			hxy = dvx.x * dwdx(k).y + dvx.y * dwdx(k).x;
			hyy = 2.f * dvx.y * dwdx(k).y - dvx.x * dwdx(k).x;

			hxx *= 2.f / 3.f;
			hyy *= 2.f / 3.f;

			txx(i) += mass(j) * hxx / rho(j);
			txx(j) += mass(i) * hxx / rho(i);
			txy(i) += mass(j) * hxy / rho(j);
			txy(j) += mass(i) * hxy / rho(i);
			tyy(i) += mass(j) * hyy / rho(j);
			tyy(j) += mass(i) * hyy / rho(i);

			// calculate SPH sum for vc, c = dvx/dx + dvy/dy + dvz/dz
			hvcc = dvx.x * dwdx(k).x + dvx.y * dwdx(k).y;
			vcc(i) += mass(j) * hvcc / rho(j);
			vcc(j) += mass(i) * hvcc / rho(i);
		}
	}

	for (rr_uint i = 0; i < ntotal; i++) {
		// viscous entropy Tds/dt = 1/2 eta/rho Tab Tab
		if constexpr (Params::visc) {
			tdsdt(i) = sqr(txx(i)) + 2.f * sqr(txy(i)) + sqr(tyy(i));
			tdsdt(i) = 0.5f * eta(i) / rho(i) * tdsdt(i);
		}

		// pressure from equation of state 
		p_art_water(rho(i), u(i), p(i), c(i));
	}

	// calculate SPH sum for pressure force -p, a/rho
	// and viscous force (eta Tab), b / rho
	// and the internal energy change de/dt due to -p/rho vc, c
	for (rr_uint k = 0; k < niac; k++) {
		i = pair_i(k);
		j = pair_j(k);
		he = 0.f;

		if constexpr (Params::pa_sph == 1) { // for sph algorithm 1
			rhoij = 1.f / (rho(i) * rho(j));
			{
				h = -(p(i) + p(j)) * dwdx(k).x;
				he += (v(j).x - v(i).x) * h;

				// viscous force
				if (Params::visc) {
					h += (eta(i) * txx(i) + eta(j) * txx(j)) * dwdx(k).x;
					h += (eta(i) * txy(i) + eta(j) * txy(j)) * dwdx(k).y;
				}

				h *= rhoij;
				a(i).x += mass(j) * h;
				a(j).x -= mass(j) * h;
			}
			{
				// pressure part
				h = -(p(i) + p(j)) * dwdx(k).y;
				he += (v(j).y - v(i).y) * h;

				// viscous force
				if (Params::visc) {
						h += (eta(i) * txy(i) + eta(j) * txy(j)) * dwdx(k).x +
							(eta(i) * tyy(i) + eta(j) * tyy(j)) * dwdx(k).y;
				}

				h *= rhoij;
				a(i).y += mass(j) * h;
				a(j).y -= mass(j) * h;
			}

			he *= rhoij;
			dedt(i) += mass(j) * he;
			dedt(j) += mass(i) * he;
		}
		else if constexpr (Params::pa_sph == 2) { // for sph algorithm 2
			{
				h = -(p(i) / sqr(rho(i)) + p(j) / sqr(rho(j))) * dwdx(k).x;
				he += (v(j).x - v(i).x) * h;

				// viscous force
				if constexpr (Params::visc) {
					h += (eta(i) * txx(i) / sqr(rho(i)) + eta(j) * txx(j) / sqr(rho(j))) * dwdx(k).x;
					h += (eta(i) * txy(i) / sqr(rho(i)) + eta(j) * txy(j) / sqr(rho(j))) * dwdx(k).y;
				}
				a(i).x += mass(j) * h;
				a(j).x -= mass(i) * h;
			}
			{
				h = -(p(i) / sqr(rho(i)) + p(j) / sqr(rho(j))) * dwdx(k).y;
				he += (v(j).y - v(i).y) * h;

				// viscous force
				if constexpr (Params::visc) {
					h += (eta(i) * txy(i) / sqr(rho(i)) + eta(j) * txy(j) / sqr(rho(j))) * dwdx(k).x +
						(eta(i) * tyy(i) / sqr(rho(i)) + eta(j) * tyy(j) / sqr(rho(j))) * dwdx(k).y;
				}
				a(i).y += mass(j) * h;
				a(j).y -= mass(i) * h;
			}
			dedt(i) += mass(j) * he;
			dedt(j) += mass(i) * he;
		}
	}

	// change of specific internal energy de/dt = T ds/ds - p/rho vc, c:
	for (rr_uint i = 0; i < ntotal; i++) {
		dedt(i) = tdsdt(i) + 0.5f * dedt(i);
	}
}