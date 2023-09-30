#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <csv-parser/csv.hpp>
#include <fmt/format.h>
#include <nlohmann/json.hpp>

#include "CommonIncl.h"
#include "Input.h"
#include "Output.h"
#include "GridUtils.h"
#include "ParamsIO.h"

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

rr_float find_depth(
	rr_uint nfluid, 
	const heap_darray<rr_float2>& r) 
{
	if (params.depth != 0) return params.depth;

	rr_float y_fluid_min = FLT_MAX;
	rr_float y_fluid_max = -FLT_MAX;

	for (rr_uint i = 0; i < nfluid; ++i) {
		y_fluid_max = std::max(r(i).y, y_fluid_max);
		y_fluid_min = std::min(r(i).y, y_fluid_min);
	}

	return y_fluid_max - y_fluid_min;
}

void fillInSPH2DParams() {
	params.format_line = "fmt: vx vy p ";
	params.hsml = params.delta * params.intf_hsml_coef;

	params.maxn = 1 << (1 + intlog2(params.ntotal));

	params.max_cells = countCells(
		params.hsml,
		params.x_mingeom,
		params.y_mingeom,
		params.x_maxgeom,
		params.y_maxgeom);

	params.mass = params.rho0 * params.delta * params.delta;

	if (params.artificial_viscosity_skf == SKF_GAUSS ||
		params.average_velocity_skf == SKF_GAUSS ||
		params.density_skf == SKF_GAUSS ||
		params.intf_skf == SKF_GAUSS)
	{
		params.cell_scale_k = 3;
	}
	else {
		params.cell_scale_k = 2;
	}
}
void postFillInModelParams(ModelParams& model_params)
{
	if (params.eos_sound_vel_method == EOS_SOUND_VEL_DAM_BREAK) {
		model_params.eos_sound_vel = params.eos_sound_vel = sqrt(200 * params.g * params.depth * params.eos_sound_vel_coef);
	}

	if (params.dt_correction_method == DT_CORRECTION_CONST_CFL) {
		model_params.dt = params.dt = params.CFL_coef * params.hsml / (params.eos_sound_vel * (1 + 1.2 * params.artificial_shear_visc));
	}
	else if (params.dt_correction_method == DT_CORRECTION_DYNAMIC) {
		throw std::runtime_error{ "not implemented: DT_CORRECTION_DYNAMIC" };
	}

	if (params.dt_correction_method == DT_CORRECTION_DYNAMIC) {
		params.maxtimestep = 0;
		if (params.step_time_estimate == 0) {
			model_params.step_time_estimate = params.step_time_estimate = 1;
		}
	}
	else {
		params.maxtimestep = static_cast<rr_uint>(params.simulation_time / params.dt);
		if (params.maxtimestep % params.save_step != 0) { // fix last save step
			params.maxtimestep = params.maxtimestep + (params.save_step - params.maxtimestep % params.save_step);
		}

		if (params.step_treatment == STEPPING_TREATMENT_TIME) {
			model_params.save_step = params.save_step = params.save_time / params.dt;
			model_params.dump_step = params.dump_step = params.dump_time / params.dt;
		}

		if (params.step_time_estimate == 0) {
			model_params.step_time_estimate = params.step_time_estimate = params.save_step;
		}
	}

	switch (params.nwm) {
	case NWM_METHOD_DYNAMIC: 
		{
			params.nwm_wave_number = 2. * params.pi / params.nwm_wave_length;
			rr_float kd = params.nwm_wave_number * params.depth;
			params.nwm_freq = sqrt(params.nwm_wave_number * params.g * tanh(kd));
			params.nwm_piston_magnitude = params.nwm_wave_magnitude * 0.5f / sqr(sinh(kd)) * (sinh(kd) * cosh(kd) + kd);
		}	
		break;
	default:
		break;
	}
}

SPH2DParams make_SPH2DParams() {
	SPH2DParams sph2DParams;

#define set_param(param) do { sph2DParams.param = params.param; } while(false)
#define set_param_not_null(param) do {if (params.param != 0) sph2DParams.param = params.param;} while(false)

	params.format_line = "fmt: p vx vy ";
	set_param(format_line);
	set_param(starttimestep);
	set_param(pi);
	set_param(g);
	set_param(TYPE_BOUNDARY);
	set_param(TYPE_NON_EXISTENT);
	set_param(TYPE_WATER);
	set_param(cell_scale_k);
	set_param(max_cells);
	set_param(maxn);
	set_param(hsml);
	set_param(mass);
	set_param_not_null(maxtimestep);
	set_param_not_null(nwm_wave_number);
	set_param_not_null(nwm_freq);
	set_param_not_null(nwm_piston_magnitude);

#undef set_param
#undef set_param_not_null
	return sph2DParams;
}

void fileInput(
	heap_darray<rr_float2>& r,	// coordinates of all particles
	heap_darray<rr_float2>& v,	// velocities of all particles
	heap_darray<rr_float>& rho,	// particle densities
	heap_darray<rr_float>& p,	// particle pressure
	heap_darray<rr_int>& itype,// particle material type 
	rr_uint& ntotal, // total particle number
	rr_uint& nfluid, // total fluid particles
	rr_uint starttimestep,
	const std::filesystem::path& experiment_directory)
{
	setupOutput(experiment_directory);

	auto particle_params = load_particle_params(experiment_directory);
	auto model_params = load_model_params(experiment_directory);
	apply_particle_params(params, particle_params);
	apply_model_params(params, model_params);
	
	printlog("Experiment name: ")(experiment_directory.stem().string())();
	printlog()(__func__)();

	ntotal = params.ntotal;
	nfluid = params.nfluid;

	r = heap_darray<rr_float2>(params.maxn);
	v = heap_darray<rr_float2>(params.maxn);
	rho = heap_darray<rr_float>(params.maxn);
	p = heap_darray<rr_float>(params.maxn);
	itype = heap_darray<rr_int>(params.maxn);

	auto particles_data_path = experiment_directory / "dump" / fmt::format("{}.csv", starttimestep);
	params.starttimestep = std::stoi(particles_data_path.stem().string());
	fillInSPH2DParams();

	std::cout << "read data...";
	csv::CSVReader reader(particles_data_path.string());

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

	particle_params.depth = params.depth = find_depth(nfluid, r);
	postFillInModelParams(model_params);

	std::cout << "...success" << std::endl;

	params_make_SPH2D_json(experiment_directory, make_SPH2DParams());
	params_make_model_json(experiment_directory, model_params);
	params_make_particles_json(experiment_directory, particle_params);
	params_make_json(experiment_directory);
	printParams();
}
