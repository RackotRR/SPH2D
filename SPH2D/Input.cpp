#include "CommonIncl.h"
#include "Input.h"
#include "EOS.h"
#include "VirtualParticles.h"
#include "Output.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <format>
#include <stdexcept>

void initConsts() {
	printlog()(__func__)();

	constexpr rr_float H = 0.1f;
	constexpr rr_float L = 2.f;
	constexpr rr_float depth = 0.7f;
	constexpr rr_float ratio = L / depth;
	constexpr rr_float tank_length = 5.f * L;
	constexpr rr_float tank_height = 2.f * depth;
	constexpr rr_uint particlesPer_d = 50;
	constexpr rr_uint particlesPer_L = static_cast<rr_uint>(particlesPer_d * ratio);
	constexpr rr_uint fluid_particles_x = static_cast<rr_uint>(particlesPer_L * tank_length / L);
	constexpr rr_uint fluid_particles_y = static_cast<rr_uint>(particlesPer_d * tank_height / depth * 0.5f);
	constexpr rr_uint fluid_particles = fluid_particles_x * fluid_particles_y;
	constexpr rr_float delta = depth / particlesPer_d;

	Params::delta = delta;
	//Params::hsml = delta * 1.2f;
	printlog("dx=dy=delta: ")(Params::delta)();
	printlog("hsml: ")(Params::hsml)();

	Params::fluid_particles_per_d = particlesPer_d;
	Params::x_fluid_particles = fluid_particles_x;
	Params::y_fluid_particles = fluid_particles_y;
	Params::x_fluid_min = 0.1f;
	Params::y_fluid_min = 0.1f;
	Params::x_fluid_max = Params::x_fluid_min + Params::x_fluid_particles * Params::delta;
	Params::y_fluid_max = Params::y_fluid_min + Params::y_fluid_particles * Params::delta;

	Params::x_maxgeom = 11.f;
	Params::x_mingeom = 0;
	Params::y_maxgeom = 4.f;
	Params::y_mingeom = 0;
	Params::beachX = 0; // is not in use now
	printlog("x_maxgeom: ")(Params::x_maxgeom)();
	printlog("x_mingeom: ")(Params::x_mingeom)();
	printlog("y_maxgeom: ")(Params::y_maxgeom)();
	printlog("y_mingeom: ")(Params::y_mingeom)();

	constexpr rr_float k = 2.f * Params::pi / L;
	constexpr rr_float kd = k * depth;
	Params::freq = sqrt(k * Params::g * tanh(kd));
	Params::A = H * 0.5f / sqr(sinh(kd)) * (sinh(kd) * cosh(kd) + kd);
	Params::H = H;
	Params::k = k;
	Params::L = L;
	Params::d = depth;
	printlog("generation wave k: ")(Params::k)();
	printlog("generation freq: ")(Params::freq)();
	printlog("generation wave height: ")(Params::H)();
	printlog("generation wave length: ")(Params::L)();
	printlog("generation fluid level (depth): ")(Params::d)();
	printlog("generator amplitude: ")(Params::A)();

	Params::normal_check_step = 500;
	Params::save_step = 500;
	Params::simulation_time = 2.f;
	Params::dt = 1e-4f;
	rr_float steps = Params::simulation_time / Params::dt;
	if (steps < 0) {
		throw std::runtime_error{ "maxtimestep error" };
	}
	Params::maxtimestep = static_cast<size_t>(steps);
	printlog("check normal step: ")(Params::normal_check_step)();
	printlog("save step: ")(Params::save_step)();
	printlog("simulation time: ")(Params::simulation_time)();
	printlog("dt: ")(Params::dt)();
	printlog("maxtimestep: ")(Params::maxtimestep)();

	Params::generator_time_wait = 0.f;
	printlog("generator time wait: ")(Params::generator_time_wait)();
}

// loading or generating initial particle information
void input(
	heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	heap_array<rr_float, Params::maxn>& mass,// particle masses
	heap_array<rr_float, Params::maxn>& rho, // particle densities
	heap_array<rr_float, Params::maxn>& p,	 // particle pressure
	heap_array<rr_float, Params::maxn>& u,	 // particle internal energy
	heap_array<rr_int, Params::maxn>& itype,	  // particle material type
	rr_uint& ntotal, // total particle number
	rr_uint& nfluid)
{
	printlog()(__func__)();

	ntotal = 0;
	nfluid = 0;
	rr_uint nvirt = 0;

	initConsts();
	printParams();

	generateParticles(r, v, mass, rho, p, u, itype, nfluid);
	virt_part(nfluid, nvirt, mass, r, v, rho, u, p, itype);
	ntotal = nfluid + nvirt;

	Params::particles_fluid = nfluid;
	Params::particles_boundary = nvirt;
	Params::particles_total = ntotal;


	printlog("nfluid: ")(nfluid)();
	printlog("nvirt: ")(nvirt)();
	printlog("ntotal: ")(ntotal)();

	printParams();
	makeParamsHeader(ntotal, nfluid, nvirt, Params::experimentName + "\\clparams.h");
}


void generateParticles(
	heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	heap_array<rr_float, Params::maxn>& mass,// particle masses
	heap_array<rr_float, Params::maxn>& rho, // particle densities
	heap_array<rr_float, Params::maxn>& p,	 // particle pressure
	heap_array<rr_float, Params::maxn>& u,	 // particle internal energy
	heap_array<rr_int, Params::maxn>& itype,	 // particle material type
	rr_uint& nfluid) // total particle number
{
	printlog()(__func__)();

	nfluid = 0;
	for (rr_uint x_i = 0; x_i < Params::x_fluid_particles; ++x_i) {
		for (rr_uint y_i = 0; y_i < Params::y_fluid_particles; ++y_i) {

			r(nfluid).x = Params::x_fluid_min + x_i * Params::delta;
			r(nfluid).y = Params::y_fluid_min + y_i * Params::delta;
			nfluid++;
		}
	}

	for (rr_uint i = 0; i < nfluid; i++) {
		v(i) = { 0.f };

		rho(i) = 1000.f;
		mass(i) = Params::delta * Params::delta * rho(i);
		u(i) = 357.1f;
		itype(i) = Params::TYPE_WATER;

		rr_float c = 0.f;
		p_art_water(rho(i), u(i), p(i), c);
	}
}


namespace {
	class ParamsHeader {
	public:
		template<typename T>
		void set_param(const char* param_name, const T& value) {
			buffer << std::format("#define params_{} {}\n", param_name, value);
		}
		template<>
		void set_param(const char* param_name, const float& value) {
			buffer << std::format("#define params_{} {:.10f}f\n", param_name, value);
			//buffer << "#define params_" << param_name << " " << std::setprecision(10) << std::fixed << val << std::endl;
		}
		template<>
		void set_param(const char* param_name, const bool& value) {
			if (value) {
				buffer << std::format("#define params_{}\n", param_name);
			}
		}

		ParamsHeader(unsigned ntotal, unsigned nfluid, unsigned nvirt) {
			buffer << "#ifndef CL_PARAMS_H" << std::endl;
			buffer << "#define CL_PARAMS_H" << std::endl << std::endl;

#define set(param) set_param(#param, param);
			using namespace Params;
			set(dim);
			set(maxn);
			set(max_neighbours);
			set(max_cells);
			set(x_maxgeom);
			set(x_mingeom);
			set(y_maxgeom);
			set(y_mingeom);
			set(L);
			set(d);
			set(freq);
			set(A);
			set(left_wall_start);
			set(left_wall_end);
			set(generator_time_wait);
			set(dt);
			set(eos);
			set(pa_sph);
			set(skf);
			set(nwm);
			set(hsml);
			set(delta);
			set(boundary_delta);
			set(summation_density);
			set(nor_density);
			set(average_velocity);
			set(visc);
			set(ex_force);
			set(self_gravity);
			set(visc_artificial);
			set(heat_artificial);
			set_param<int>("TYPE_BOUNDARY", TYPE_BOUNDARY);
			set_param<int>("TYPE_NON_EXISTENT", TYPE_NON_EXISTENT);
			set_param<int>("TYPE_WATER", TYPE_WATER);
			set(pi);
			set(g);
			set(ntotal);
			set(nfluid);
			set(nvirt);
#undef set

			buffer << std::endl << "#endif" << std::endl;
		}

		std::string string() {
			return buffer.str();
		}
	private:
		std::stringstream buffer;
	};
}

void makeParamsHeader(unsigned ntotal, unsigned nfluid, unsigned nvirt, std::string path) {
	::ParamsHeader header(ntotal, nfluid, nvirt);
	std::string params = header.string();
	std::ofstream stream(path);
	stream << params;
}