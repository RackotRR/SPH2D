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
#include "ConsistencyCheck.h"

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

static rr_uint countCells(rr_float hsml) {
	rr_uint x_id = get_cell_coord_from_particle_coord(params.x_maxgeom, params.x_mingeom);
	rr_uint y_id = get_cell_coord_from_particle_coord(params.y_maxgeom, params.y_mingeom);
	rr_uint z_id;
	rr_uint max_id;
	rr_uint id;

	if (params.dim == 3) {
		z_id = get_cell_coord_from_particle_coord(params.z_maxgeom, params.z_mingeom);
		max_id = 1 << 10;
		id = get_cell_idx_by_cell_coord3({ x_id, y_id, z_id });
	}
	else
	{
		z_id = 0;
		max_id = 1 << 16;
		id = get_cell_idx_by_cell_coord2({ x_id, y_id });
	}

	if (x_id > max_id || y_id > max_id || z_id >= max_id) {
		std::string x_constraint = fmt::format("{}/{}", x_id, max_id);
		std::string y_constraint = fmt::format("{}/{}", y_id, max_id);
		std::string z_constraint = fmt::format("{}/{}", z_id, max_id);
		std::string constraints = fmt::format("({}; {}; {})", x_constraint, y_constraint, z_constraint);
		throw std::runtime_error{ "can't make grid with so many cells: " + constraints };
	}

	rr_uint max_cells = 1 << (intlog2(id) + 1);

	return max_cells;
}

template<typename rr_floatn>
rr_float find_depth(const heap_darray<rr_floatn>& r)
{
	if (params.depth != 0) return params.depth;

	rr_float y_fluid_min = FLT_MAX;
	rr_float y_fluid_max = -FLT_MAX;

	for (rr_uint i = 0; i < params.nfluid; ++i) {
		y_fluid_max = std::max(r(i).y, y_fluid_max);
		y_fluid_min = std::min(r(i).y, y_fluid_min);
	}

	return y_fluid_max - y_fluid_min;
}

static void fillInComputingParams() {
	printlog()(__func__)();

	params.hsml = /*sqrt(params.dim) * */params.delta * params.intf_hsml_coef;

	params.maxn = 1 << (1 + intlog2(params.ntotal));

	params.max_cells = countCells(params.hsml);
	
	params.mass = params.rho0 * powun(params.delta, params.dim);

	params.cell_scale_k = get_cell_scale_k({
		params.artificial_pressure_skf,
		params.artificial_viscosity_skf,
		params.average_velocity_skf,
		params.intf_skf,
		params.density_skf
	});
}
static void postFillInModelParams(ModelParams& model_params)
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

static ComputingParams make_ComputingParams() {
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

static void printLogVersions() {
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
}

template<typename rr_floatn>
void loadArrays(
	const std::filesystem::path& initial_dump_path,
	heap_darray<rr_floatn>& r,
	heap_darray<rr_floatn>& v,
	heap_darray<rr_float>& rho,
	heap_darray<rr_float>& p,
	heap_darray<rr_int>& itype)
{
	std::cout << "read data...";
	csv::CSVReader reader(initial_dump_path.string());
	size_t j = 0;

	for (const auto& row : reader) {
		r(j).x = row["x"].get<rr_float>();
		r(j).y = row["y"].get<rr_float>();
		if constexpr (std::is_same<rr_floatn, rr_float3>::value == true) {
			r(j).z = row["z"].get<rr_float>();
		}
		itype(j) = row["itype"].get<rr_int>();
		v(j).x = row["vx"].get<rr_float>();
		v(j).y = row["vy"].get<rr_float>();
		if constexpr (std::is_same<rr_floatn, rr_float3>::value == true) {
			v(j).z = row["vz"].get<rr_float>();
		}
		rho(j) = row["rho"].get<rr_float>();
		p(j) = row["p"].get<rr_float>();
		++j;
	}

	if (j != params.ntotal) {
		throw std::runtime_error{ "dump corrupted: dump ntotal doesn't match params.ntotal!" };
	}

	std::cout << "...success" << std::endl;
}

void fileInput(
	vheap_darray_floatn& r_var,	// coordinates of all particles
	vheap_darray_floatn& v_var,	// velocities of all particles
	heap_darray<rr_float>& rho,	// particle densities
	heap_darray<rr_float>& p,	// particle pressure
	heap_darray<rr_int>& itype,// particle material type
	const std::filesystem::path& initial_dump_path,
	const std::filesystem::path& experiment_directory)
{
	params.experiment_name = experiment_directory.stem().string();
	RRSPHOutput::instance().initialize(experiment_directory);

	auto particle_params = load_particle_params(experiment_directory);
	auto model_params = load_model_params(experiment_directory);
	apply_particle_params(params, particle_params);
	apply_model_params(params, model_params);
	
	printlog()("Experiment name: ")(params.experiment_name)();
	printLogVersions();

	params.start_simulation_time = std::stod(initial_dump_path.stem().string());
	fillInComputingParams();

	// global init dimensions for variant-based arrays
	printlog("set experiment dimensions: ")(params.dim)();
	vheap_darray_floatn::set_dimenstions(params.dim);
	vheap_darray_floatn_md::set_dimenstions(params.dim);

	printlog("allocate memory for ")(params.maxn)(" particles")();
	r_var = vheap_darray_floatn(params.maxn);
	v_var = vheap_darray_floatn(params.maxn);
	rho = heap_darray<rr_float>(params.maxn);
	p = heap_darray<rr_float>(params.maxn);
	itype = heap_darray<rr_int>(params.maxn);

	if (params.dim == 2) {
		auto& r = r_var.get_flt2();
		auto& v = v_var.get_flt2();
		loadArrays(initial_dump_path, r, v, rho, p, itype);
		particle_params.depth = params.depth = find_depth(r);
	}
	else if (params.dim == 3) {
		auto& r = r_var.get_flt3();
		auto& v = v_var.get_flt3();
		loadArrays(initial_dump_path, r, v, rho, p, itype);
		particle_params.depth = params.depth = find_depth(r);
	}
	else {
		printlog("invalid params.dim value: ")(params.dim)();
		exit(-1);
	}

	if (!check_particles_are_within_boundaries(r_var, itype)) {
		throw std::runtime_error{ "failed consistency check on input" };
	}

	postFillInModelParams(model_params);

	params_make_computing_json(experiment_directory, make_ComputingParams());
	params_make_model_json(experiment_directory, model_params);
	params_make_particles_json(experiment_directory, particle_params);
	params_make_json(experiment_directory);

	if (RRSPH_GetSpecificVersionName() == "RRSPH_CL") {
		params_make_header(std::filesystem::current_path() / "cl" / "clparams.h");
	}
}
