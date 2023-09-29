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
	params.hsml = params.delta * params.intf_hsml_coef;

	params.maxn = 1 << (1 + intlog2(params.ntotal));

	params.max_cells = countCells(
		params.hsml,
		params.x_mingeom,
		params.y_mingeom,
		params.x_maxgeom,
		params.y_maxgeom);

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
void postFillInSPH2DParams()
{
	if (params.eos_sound_vel_method == EOS_SOUND_VEL_DAM_BREAK) {
		params.eos_sound_vel = sqrt(200 * params.g * params.depth * params.eos_sound_vel_coef);
	}

	if (params.dt_correction_method == DT_CORRECTION_CONST_CFL) {
		params.dt = params.CFL_coef * params.hsml / (params.eos_sound_vel * (1 + 1.2 * params.artificial_shear_visc));
	}
	else if (params.dt_correction_method == DT_CORRECTION_DYNAMIC) {
		throw std::runtime_error{ "not implemented: DT_CORRECTION_DYNAMIC" };
	}

	if (params.dt_correction_method == DT_CORRECTION_DYNAMIC) {
		params.maxtimestep = 0;
	}
	else {
		params.maxtimestep = static_cast<rr_uint>(params.simulation_time / params.dt);
		if (params.maxtimestep % params.save_step != 0) { // fix last save step
			params.maxtimestep = params.maxtimestep + (params.save_step - params.maxtimestep % params.save_step);
		}

		if (params.step_treatment == STEPPING_TREATMENT_TIME) {
			params.save_step = params.save_time / params.dt;
			params.dump_step = params.dump_time / params.dt;
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

void print_SPH2DParams(const std::filesystem::path& experiment_directory) {
	std::ofstream stream{ experiment_directory / "SPH2DParams.json" };

	nlohmann::json json;
#define print_param(param) json[#param] = params.param;
#define print_not_null(param) if (params.param != 0) json[#param] = params.param;

	params.format_line = "fmt: vx vy p ";
	print_param(format_line)
	print_param(starttimestep)
	print_param(pi)
	print_param(g)
	print_param(TYPE_BOUNDARY)
	print_param(TYPE_NON_EXISTENT)
	print_param(TYPE_WATER)
	print_param(cell_scale_k)
	print_param(max_cells)
	print_param(maxn)
	print_param(hsml)
	print_param(depth)
	print_not_null(maxtimestep)
	print_not_null(nwm_wave_number)
	print_not_null(nwm_freq)
	print_not_null(nwm_piston_magnitude)

#undef print_param
#undef print_not_null

	stream << json.dump(4) << std::endl;
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
	params = load_experiment_params(experiment_directory);
	
	printlog("Experiment name: ")(params.experiment_name)();
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

	params.depth = find_depth(nfluid, r);
	postFillInSPH2DParams();

	std::cout << "...success" << std::endl;

	print_SPH2DParams(experiment_directory);
	params_make_json(experiment_directory / "Params.json");
	printParams();
}
