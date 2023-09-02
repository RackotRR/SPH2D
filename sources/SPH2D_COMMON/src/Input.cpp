#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <csv-parser/csv.hpp>

#include "CommonIncl.h"
#include "Input.h"
#include "VirtualParticles.h"
#include "Output.h"
#include "GridUtils.h"


static rr_uint setup_virt_part(
	rr_float delta,
	rr_float hsml,
	rr_float x_fluid_min,
	rr_float y_fluid_min,
	rr_float x_fluid_max,
	rr_float y_maxgeom)
{
	printlog()(__func__)();

	params.boundary_layers_num = 3;
	params.boundary_delta = delta * 1.f;

	constexpr rr_float spacing = 2;

	params.x_boundary_min = x_fluid_min - spacing * hsml;
	params.y_boundary_min = y_fluid_min - spacing * hsml;
	params.x_boundary_max = x_fluid_max + spacing * hsml;
	params.y_boundary_max = y_maxgeom;

	params.beach_x = 0; // is not in use now
	params.nwm = 0;
	params.sbt = 0;
	params.generator_time_wait = 0.5f;

	return count_virt_part_num();
}

rr_uint countCells(
	rr_float hsml,
	rr_float x_mingeom,
	rr_float y_mingeom,
	rr_float x_maxgeom,
	rr_float y_maxgeom)
{
	rr_uint x_id = get_cell_x_from_coordinate(params.x_maxgeom);
	rr_uint y_id = get_cell_x_from_coordinate(params.y_maxgeom);

	if (x_id > UINT16_MAX || y_id > UINT16_MAX) {
		throw std::runtime_error{ "can't make grid with so many cells" };
	}

	rr_uint id = get_cell_idx_by_cell_xy(x_id, y_id);
	rr_uint max_cells = 1 << (intlog2(id) + 1);

	return max_cells;
}

void loadDefaultParams() {
	printlog()(__func__)();
	
	constexpr rr_float H = 0.0f;
	constexpr rr_float L = 1.2f;
	constexpr rr_float depth = 0.6f;
	constexpr rr_float ratio = L / depth;
	constexpr rr_float tank_length = 3.22f;
	constexpr rr_float tank_height = 2.f;
	constexpr rr_uint particlesPer_d = 125;
	constexpr rr_uint particlesPer_L = static_cast<rr_uint>(particlesPer_d * ratio);
	constexpr rr_uint fluid_particles_x = static_cast<rr_uint>(particlesPer_L);
	constexpr rr_uint fluid_particles_y = static_cast<rr_uint>(particlesPer_d);
	constexpr rr_uint fluid_particles = fluid_particles_x * fluid_particles_y;
	constexpr rr_float delta = depth / particlesPer_d;

	params.delta = delta;
	params.hsml = delta * 1.f;

	params.mass = 1000 * sqr(delta);

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

	params.max_cells = countCells(
		params.hsml,
		params.x_mingeom,
		params.y_mingeom,
		params.x_maxgeom,
		params.y_maxgeom);

	constexpr rr_float k = 2.f * params.pi / L;
	constexpr rr_float kd = k * depth;
	params.freq = sqrt(k * params.g * tanh(kd));
	params.piston_amp = H * 0.5f / sqr(sinh(kd)) * (sinh(kd) * cosh(kd) + kd);
	params.wave_amp = H;
	params.wave_number = k;
	params.wave_length = L;
	params.depth = depth;

	params.density_skf = 1;
	params.int_force_skf = 4;
	params.artificial_viscosity_skf = 4;
	params.average_velocity_skf = 3;
	params.cell_scale_k = get_cell_scale_k(
		params.density_skf,
		params.int_force_skf,
		params.artificial_viscosity_skf,
		params.average_velocity_skf);;

	params.eos_csqr_k = 1;
	params.average_velocity = true;
	params.average_velocity_epsilon = 0.01f;
	params.artificial_shear_visc = 0.01f;
	params.artificial_bulk_visc = 0;

	params.save_step = 1000;
	params.dump_step = 10 * params.save_step;
	params.normal_check_step = params.save_step;
	params.simulation_time = 3.f;
	params.CFL_coef = 0.125f;
	params.dt = params.CFL_coef * params.hsml / (2.2f * sqrt(200 * params.g * params.depth * params.eos_csqr_k));
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
		itype(i) = params.TYPE_WATER;
		p(i) = 0;
	}
}

// loading or generating initial particle information
void input(
	heap_darray<rr_float2>& r,	// coordinates of all particles
	heap_darray<rr_float2>& v,	// velocities of all particles
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
	rho = heap_darray<rr_float>(params.maxn);
	p = heap_darray<rr_float>(params.maxn);
	itype = heap_darray<rr_int>(params.maxn);

	ntotal = params.ntotal;
	nfluid = params.nfluid;

	generateParticles(r, v, rho, p, itype);
	virt_part(nfluid, r, v, rho, p, itype);
	printParams();
}


void fileInput(
	heap_darray<rr_float2>& r,	// coordinates of all particles
	heap_darray<rr_float2>& v,	// velocities of all particles
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

