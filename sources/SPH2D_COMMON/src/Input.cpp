#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <csv-parser/csv.hpp>

#include "CommonIncl.h"
#include "Input.h"
#include "VirtualParticles.h"
#include "Output.h"

static rr_uint setup_virt_part(
	rr_float delta,
	rr_float hsml,
	rr_float x_fluid_min,
	rr_float y_fluid_min,
	rr_float x_fluid_max,
	rr_float y_maxgeom)
{
	printlog()(__func__)();

	params.boundary_layers_num = 2;
	params.boundary_delta = delta * 2.f;

	params.x_boundary_min = x_fluid_min - 2 * hsml;
	params.y_boundary_min = y_fluid_min - 2 * hsml;
	params.x_boundary_max = x_fluid_max + 2 * hsml;
	params.y_boundary_max = y_maxgeom;

	params.beach_x = 0; // is not in use now
	params.nwm = 2;
	params.generator_time_wait = 0.5f;

	return count_virt_part_num();
}

void loadDefaultParams() {
	printlog()(__func__)();
	
	constexpr rr_float H = 0.2f;
	constexpr rr_float L = 8.975979804992676f;
	constexpr rr_float depth = 0.7f;
	constexpr rr_float ratio = L / depth;
	constexpr rr_float tank_length = 3.5f * L;
	constexpr rr_float tank_height = 2.f;
	constexpr rr_uint particlesPer_d = 50;
	constexpr rr_uint particlesPer_L = static_cast<rr_uint>(particlesPer_d * ratio);
	constexpr rr_uint fluid_particles_x = static_cast<rr_uint>(particlesPer_L / L * tank_length);
	constexpr rr_uint fluid_particles_y = static_cast<rr_uint>(particlesPer_d);
	constexpr rr_uint fluid_particles = fluid_particles_x * fluid_particles_y;
	constexpr rr_float delta = depth / particlesPer_d;

	params.delta = delta;
	params.hsml = delta * 1.2f;

	params.fluid_particles_per_d = particlesPer_d;
	params.x_fluid_particles = fluid_particles_x;
	params.y_fluid_particles = fluid_particles_y;
	params.x_fluid_min = 0.f;
	params.y_fluid_min = 0.f;
	params.x_fluid_max = params.x_fluid_min + (particlesPer_L / L * tank_length) * params.delta;
	params.y_fluid_max = params.y_fluid_min + (particlesPer_d / depth * tank_height) * params.delta;

	params.x_maxgeom = tank_length + 0.1;
	params.x_mingeom = -0.3;
	params.y_maxgeom = tank_height + 0.2;
	params.y_mingeom = -0.3;

	rr_uint nvirt = setup_virt_part(
		params.delta,
		params.hsml,
		params.x_fluid_min,
		params.y_fluid_min,
		params.x_fluid_max,
		params.y_maxgeom);

	params.nvirt = nvirt;
	params.nfluid = fluid_particles;
	params.ntotal = params.nvirt + params.nfluid;
	params.maxn = 1 << (1 + intlog2(params.ntotal));
	params.max_cells = params.max_neighbours * params.maxn;

	constexpr rr_float k = 2.f * params.pi / L;
	constexpr rr_float kd = k * depth;
	params.freq = sqrt(k * params.g * tanh(kd));
	params.piston_amp = H * 0.5f / sqr(sinh(kd)) * (sinh(kd) * cosh(kd) + kd);
	params.wave_amp = H;
	params.wave_number = k;
	params.wave_length = L;
	params.depth = depth;

	params.eos_csqr_k = 1;
	params.average_velocity = true;
	params.average_velocity_epsilon = 0.3f;
	params.int_force_kernel = false;

	params.save_step = 1000;
	params.dump_step = 10 * params.save_step;
	params.normal_check_step = params.save_step;
	params.simulation_time = 15.f;
	constexpr rr_float k_CFL = 0.125f;
	params.dt = k_CFL * params.hsml / (2.2f * sqrt(200 * params.g * params.depth * params.eos_csqr_k));
	params.maxtimestep = static_cast<rr_uint>(params.simulation_time / params.dt);
	if (params.maxtimestep % params.save_step != 0) { // fix last save step
		params.maxtimestep = params.maxtimestep + (params.save_step - params.maxtimestep % params.save_step);
	}

	params.print_time_est_step = 500;
	params.local_threads = 256;
}

static void generateParticles(
	heap_darray<rr_float2>& r,	// coordinates of all particles
	heap_darray<rr_float2>& v,	// velocities of all particles
	heap_darray<rr_float>& mass,// particle masses
	heap_darray<rr_float>& rho, // particle densities
	heap_darray<rr_float>& p,	 // particle pressure
	heap_darray<rr_int>& itype)	 // particle material type
{
	printlog()(__func__)();

	rr_uint nfluid = 0;
	for (rr_uint x_i = 0; x_i < params.x_fluid_particles; ++x_i) {
		for (rr_uint y_i = 0; y_i < params.y_fluid_particles; ++y_i) {

			r(nfluid).x = params.x_fluid_min + x_i * params.delta;
			r(nfluid).y = params.y_fluid_min + y_i * params.delta;
			nfluid++;
		}
	}
	if (nfluid != params.nfluid) {
		throw std::runtime_error{"Generated nfluid doesn't equal to estimated nfluid" };
	}

	for (rr_uint i = 0; i < nfluid; i++) {
		v(i) = { 0.f };

		rho(i) = 1000.f;
		mass(i) = params.delta * params.delta * rho(i);
		itype(i) = params.TYPE_WATER;
		p(i) = 0;
	}
}

// loading or generating initial particle information
void input(
	heap_darray<rr_float2>& r,	// coordinates of all particles
	heap_darray<rr_float2>& v,	// velocities of all particles
	heap_darray<rr_float>& mass,	// particle masses
	heap_darray<rr_float>& rho,	// particle densities
	heap_darray<rr_float>& p,	// particle pressure
	heap_darray<rr_int>& itype,	// particle material type 
	rr_uint& ntotal, // total particle number
	rr_uint& nfluid, // total fluid particles
	bool load_default_params)
{
	setupOutput();
	printlog("Experiment name: ")(params.experiment_name)();
	printlog()(__func__)();

	if (load_default_params) {
		loadDefaultParams();
	}

	r = heap_darray<rr_float2>(params.maxn); 
	v = heap_darray<rr_float2>(params.maxn);
	mass = heap_darray<rr_float>(params.maxn);
	rho = heap_darray<rr_float>(params.maxn);
	p = heap_darray<rr_float>(params.maxn);
	itype = heap_darray<rr_int>(params.maxn);

	ntotal = params.ntotal;
	nfluid = params.nfluid;

	generateParticles(r, v, mass, rho, p, itype);
	virt_part(nfluid, r, v, mass, rho, p, itype);
	printParams();
}


void fileInput(
	heap_darray<rr_float2>& r,	// coordinates of all particles
	heap_darray<rr_float2>& v,	// velocities of all particles
	heap_darray<rr_float>& mass,	// particle masses
	heap_darray<rr_float>& rho,	// particle densities
	heap_darray<rr_float>& p,	// particle pressure
	heap_darray<rr_int>& itype,// particle material type 
	rr_uint& ntotal, // total particle number
	rr_uint& nfluid, // total fluid particles
	std::string particles_path,
	std::string params_path)
{
	setupOutput();
	printlog("Experiment name: ")(params.experiment_name)();
	printlog()(__func__)();

	if (!params_path.empty()) {
		params.load(params_path);
	}

	r = heap_darray<rr_float2>(params.maxn);
	v = heap_darray<rr_float2>(params.maxn);
	mass = heap_darray<rr_float>(params.maxn, 1000.f * params.delta * params.delta);
	rho = heap_darray<rr_float>(params.maxn);
	p = heap_darray<rr_float>(params.maxn);
	itype = heap_darray<rr_int>(params.maxn);

	ntotal = params.ntotal;
	nfluid = params.nfluid;

	std::string particles_filename = std::filesystem::path(particles_path).stem().string();
	params.starttimestep = atoi(particles_filename.c_str()) + 1;

	std::cout << "read data...";

	auto read_path = std::filesystem::current_path() / params.experiment_name / "dump" / (particles_filename + ".csv");
	csv::CSVReader reader(read_path.string());

	size_t j = 0;
	for (const auto& row : reader) {
		r(j).x = row["x"].get<float>();
		r(j).y = row["y"].get<float>();
		itype(j) = row["itype"].get<int>();
		v(j).x = row["vx"].get<float>();
		v(j).y = row["vy"].get<float>();
		rho(j) = row["rho"].get<float>();
		p(j) = row["p"].get<float>();
		++j;
	}

	if (j != params.ntotal) {
		throw std::runtime_error{ "dump corrupted: dump ntotal doesn't match params.ntotal!" };
	}
	std::cout << "...success" << std::endl;

	printParams();
}

