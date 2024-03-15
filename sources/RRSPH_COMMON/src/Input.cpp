#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <csv-parser/csv.hpp>
#include <fmt/format.h>
#include <nlohmann/json.hpp>
#include <cfloat>

#include "CommonIncl.h"
#include "Input.h"
#include "Output.h"
#include "GridUtils.h"
#include "ParamsIO.h"
#include "RRSPHCOMMONVersion.h"

static rr_uint get_cell_scale_k(rr_uint skf) {
	switch (skf) {
	case SKF_CUBIC: return 2;
	case SKF_GAUSS: return 3;
	case SKF_WENDLAND: return 2;
	case SKF_DESBRUN: return 2;
	default: 
		throw std::runtime_error{fmt::format("get_cell_scale_k: unknown skf {}", skf)};
	}
}

static rr_float get_cell_scale_k(std::vector<rr_uint> skf) {
	assert(!skf.empty());

	std::vector<rr_uint> cell_scale_k;
	std::transform(skf.begin(), skf.end(), 
		std::back_inserter(cell_scale_k), 
		[](rr_uint skf_type){
			return get_cell_scale_k(skf_type);
		});

	auto[min, max] = std::minmax_element(cell_scale_k.begin(), cell_scale_k.end());
	if (*min != *max) { // cell_scale_k must be equal for all skf
		throw std::runtime_error{ "selected incompatible skf" };
	}
	return *max;
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

void fillInComputingParams() {
	params.hsml = params.delta * params.intf_hsml_coef;

	params.maxn = 1 << (1 + intlog2(params.ntotal));

	params.max_cells = countCells(
		params.hsml,
		params.x_mingeom,
		params.y_mingeom,
		params.x_maxgeom,
		params.y_maxgeom);

	params.mass = params.rho0 * params.delta * params.delta;

	params.cell_scale_k = get_cell_scale_k({
		params.artificial_pressure_skf,
		params.artificial_viscosity_skf,
		params.average_velocity_skf,
		params.intf_skf,
		params.density_skf
	});
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
		model_params.dt = 0;
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

ComputingParams make_ComputingParams() {
	ComputingParams computing_params;

#define set_param(param) do { computing_params.param = params.param; } while(false)
#define set_param_not_null(param) do {if (params.param != 0) computing_params.param = params.param;} while(false)

	set_param(start_simulation_time);
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
	set_param_not_null(nwm_wave_number);
	set_param_not_null(nwm_freq);
	set_param_not_null(nwm_piston_magnitude);

	computing_params.params_version_major = RRSPH_PARAMS_VERSION_MAJOR;
	computing_params.params_version_minor = RRSPH_PARAMS_VERSION_MINOR;
	computing_params.params_version_patch = RRSPH_PARAMS_VERSION_PATCH;

	computing_params.RRSPH_common_version_major = RRSPH_COMMON_VERSION_MAJOR;
	computing_params.RRSPH_common_version_minor = RRSPH_COMMON_VERSION_MINOR;
	computing_params.RRSPH_common_version_patch = RRSPH_COMMON_VERSION_PATCH;

	computing_params.RRSPH_version_major = RRSPH_VERSION_MAJOR;
	computing_params.RRSPH_version_minor = RRSPH_VERSION_MINOR;
	computing_params.RRSPH_version_patch = RRSPH_VERSION_PATCH;

	computing_params.RRSPH_specific_version_major = RRSPH_GetSpecificVersionMajor();
	computing_params.RRSPH_specific_version_minor = RRSPH_GetSpecificVersionMinor();
	computing_params.RRSPH_specific_version_patch = RRSPH_GetSpecificVersionPatch();
	computing_params.RRSPH_specific_version_name = RRSPH_GetSpecificVersionName();

#undef set_param
#undef set_param_not_null
	return computing_params;
}

void fileInput(
	heap_darray<rr_float2>& r,	// coordinates of all particles
	heap_darray<rr_float2>& v,	// velocities of all particles
	heap_darray<rr_float>& rho,	// particle densities
	heap_darray<rr_float>& p,	// particle pressure
	heap_darray<rr_int>& itype,// particle material type 
	rr_uint& ntotal, // total particle number
	rr_uint& nfluid, // total fluid particles
	const std::filesystem::path& initial_dump_path,
	const std::filesystem::path& experiment_directory)
{
	SPH2DOutput::instance().initialize(experiment_directory);

	auto particle_params = load_particle_params(experiment_directory);
	auto model_params = load_model_params(experiment_directory);
	apply_particle_params(params, particle_params);
	apply_model_params(params, model_params);
	
	printlog()("Experiment name: ")(experiment_directory.stem().string())();
	printlog(fmt::format("RRSPH v{}.{}.{}", 
		RRSPH_VERSION_MAJOR, 
		RRSPH_VERSION_MINOR,
		RRSPH_VERSION_PATCH))();
	printlog(fmt::format("Use {} v{}.{}.{}",
		RRSPH_GetSpecificVersionName(),
		RRSPH_GetSpecificVersionMajor(),
		RRSPH_GetSpecificVersionMinor(),
		RRSPH_GetSpecificVersionPatch()))();
	printlog(fmt::format("Params v{}.{}.{}", 
		RRSPH_PARAMS_VERSION_MAJOR,
		RRSPH_PARAMS_VERSION_MINOR,
		RRSPH_PARAMS_VERSION_PATCH))();
	printlog(fmt::format("Common v{}.{}.{}",
		RRSPH_COMMON_VERSION_MAJOR,
		RRSPH_COMMON_VERSION_MINOR,
		RRSPH_COMMON_VERSION_PATCH))();
	printlog()(__func__)();

	ntotal = params.ntotal;
	nfluid = params.nfluid;

	r = heap_darray<rr_float2>(params.maxn);
	v = heap_darray<rr_float2>(params.maxn);
	rho = heap_darray<rr_float>(params.maxn);
	p = heap_darray<rr_float>(params.maxn);
	itype = heap_darray<rr_int>(params.maxn);

	params.start_simulation_time = std::stod(initial_dump_path.stem().string());
	fillInComputingParams();

	std::cout << "read data...";
	csv::CSVReader reader(initial_dump_path.string());

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

	params_make_computing_json(experiment_directory, make_ComputingParams());
	params_make_model_json(experiment_directory, model_params);
	params_make_particles_json(experiment_directory, particle_params);
	params_make_json(experiment_directory);
	params_make_header(std::filesystem::current_path() / "cl" / "clparams.h");
}
