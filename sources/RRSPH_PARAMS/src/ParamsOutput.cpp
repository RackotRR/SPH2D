#include <nlohmann/json.hpp>
#include <fmt/format.h>
#include <sstream>
#include <fstream>

#include "Params.h"
#include "ParamsGeneration.h"
#include "Version.h"


namespace ParamsGeneration {
	class ParamsSerializer {
	public:
		virtual void set_param(const char* param_name, const std::string& value) {}
		virtual void set_param(const char* param_name, rr_int value) {}
		virtual void set_param(const char* param_name, rr_uint value) {}
		virtual void set_param(const char* param_name, rr_float value) {}
		virtual void set_param(const char* param_name, bool value) {}
		virtual void print_before(const std::string& before) {}
		virtual void print_after(const std::string& after) {}

		void serialize(const ExperimentParams& params) {
			print_before("");
#define set(param) set_param(#param, params.param);
			set(params_version_major);
			set(params_version_minor);
			set(params_version_patch);
			set(RRSPH_version_major);
			set(RRSPH_version_minor);
			set(RRSPH_version_patch);
			set(dim);
			set(maxn);
			set(max_neighbours);
			set(max_cells);
			set(x_maxgeom);
			set(x_mingeom);
			set(y_maxgeom);
			set(y_mingeom);
			set(z_mingeom);
			set(z_maxgeom);
			set(x_fluid_particles);
			set(y_fluid_particles);
			set(x_fluid_min);
			set(y_fluid_min);
			set(x_fluid_max);
			set(y_fluid_max);
			set(x_boundary_min);
			set(y_boundary_min);
			set(x_boundary_max);
			set(y_boundary_max);
			set(nfluid);
			set(nvirt);
			set(ntotal);
			set(nwm_wave_length);
			set(depth);
			set(nwm_freq);
			set(nwm_piston_magnitude);
			set(nwm_wave_magnitude);
			set(nwm_wave_number);
			set(beach_x);
			set(nwm_particles_start);
			set(nwm_particles_end);
			set(nwm_time_start);
			set(CFL_coef);
			set(dt);
			set(dt_correction_method);
			set(simulation_time);
			set(local_threads);
			set(eos_sound_vel_coef);
			set(eos_sound_vel_method);
			set(eos_sound_vel);
			set(intf_sph_approximation);
			set(intf_hsml_coef);
			set(artificial_pressure);
			set(artificial_pressure_skf);
			set(artificial_pressure_coef);
			set(artificial_pressure_index);
			set(density_skf);
			set(density_delta_sph_coef);
			set(intf_skf);
			set(average_velocity_skf);
			set(artificial_viscosity_skf);
			set(cell_scale_k);
			set(nwm);
			set(boundary_layers_num);
			set(boundary_treatment);
			set(use_chess_order);
			set(hsml);
			set(delta);
			set(boundary_delta);
			set(density_treatment);
			set(density_normalization);
			set(average_velocity);
			set(average_velocity_coef);
			set(visc);
			set(visc_coef);
			set(artificial_viscosity);
			set(artificial_shear_visc);
			set(artificial_bulk_visc);
			set_param("TYPE_BOUNDARY", (rr_int)params.TYPE_BOUNDARY);
			set_param("TYPE_NON_EXISTENT", (rr_int)params.TYPE_NON_EXISTENT);
			set_param("TYPE_WATER", (rr_int)params.TYPE_WATER);
			set(mass);
			set(rho0);
			set(consistency_check);
			set(consistency_treatment);
			set(start_simulation_time);
			set(save_velocity);
			set(save_pressure);
			set(save_density);
			set(save_every_step);
			set(save_time);
			set(dump_time);
			set(step_time_estimate);
			set(pi);
			set(g);
			set(experiment_name);
#undef set
			print_after("");
		}
	};


	class ParamsHeader : public ParamsSerializer {
	public:
		void set_param(const char* param_name, const std::string& value) override {
			buffer << fmt::format("#define params_{} \"{}\"\n", param_name, value);
		}
		void set_param(const char* param_name, rr_int value) override {
			buffer << fmt::format("#define params_{} {}\n", param_name, value);
		}
		void set_param(const char* param_name, rr_uint value) override {
			buffer << fmt::format("#define params_{} {}\n", param_name, value);
		}
		void set_param(const char* param_name, rr_float value) override {
			std::string postfix = "";
			if (std::is_same_v<rr_float, float>) postfix = "f";
			buffer << fmt::format("#define params_{} {:.15f}{}\n", param_name, value, postfix);
		}
		void set_param(const char* param_name, bool value) override {
			if (value) {
				buffer << fmt::format("#define params_{}\n", param_name);
			}
			else {
				buffer << fmt::format("#undef params_{}\n", param_name);
			}
		}
		void print_before(const std::string& before) override {
			buffer << "#ifndef CL_PARAMS_H" << std::endl;
			buffer << "#define CL_PARAMS_H" << std::endl << std::endl;
		}
		void print_after(const std::string& after) override {
			buffer << std::endl << "#endif" << std::endl;
		}

		std::string string() {
			return buffer.str();
		}
	private:
		std::stringstream buffer;
	};


	class ParamsJson : public ParamsSerializer {
	public:
		void set_param(const char* param_name, const std::string& value) override {
			json[param_name] = value;
		}
		void set_param(const char* param_name, rr_int value) override {
			json[param_name] = value;
		}
		void set_param(const char* param_name, rr_uint value) override {
			json[param_name] = value;
		}
		void set_param(const char* param_name, rr_float value) override {
			json[param_name] = value;
		}
		void set_param(const char* param_name, bool value) override {
			json[param_name] = value;
		}
		std::string string() {
			return json.dump(4);
		}
	private:
		nlohmann::json json;
	};

	class ParamsGeneratorClass : public ParamsSerializer {
	public:
		void set_param(const char* param_name, const std::string& value) override {
			buffer << "public string " << param_name << " { get; set; }" << std::endl;
		}
		void set_param(const char* param_name, rr_int value) override {
			buffer << "public int " << param_name << " { get; set; }" << std::endl;
		}
		void set_param(const char* param_name, rr_uint value) override {
			buffer << "public uint " << param_name << " { get; set; }" << std::endl;
		}
		void set_param(const char* param_name, rr_float value) override {
			buffer << "public float " << param_name << " { get; set; }" << std::endl;
		}
		void set_param(const char* param_name, bool value) override {
			buffer << "public bool " << param_name << " { get; set; }" << std::endl;
		}
		void print_before(const std::string& before) override {
			buffer << "namespace SPH2DParamsGenerator {" << std::endl;
			buffer << "class ExperimentParams {" << std::endl;
		}
		void print_after(const std::string& after) override {
			buffer << "}" << std::endl << "}" << std::endl;
		}

		std::string string() {
			return buffer.str();
		}
	private:
		std::stringstream buffer;
	};


	void makeHeader(const std::filesystem::path& path, const ExperimentParams& params) {
		ParamsHeader header;
		header.serialize(params);
		std::string params_str = header.string();
		std::ofstream stream{ path };
		stream << params_str << std::endl;
	}
	void makeJson(const std::filesystem::path& path, const ExperimentParams& params) {
		ParamsJson json;
		json.serialize(params);
		std::string params_str = json.string();
		std::ofstream stream{ path };
		stream << params_str << std::endl;
	}
	void makeParamsGeneratorClass(const std::filesystem::path& path) {
		ParamsGeneratorClass generatorClass;
		generatorClass.serialize(ExperimentParams{});
		std::string params_str = generatorClass.string();
		std::ofstream stream{ path };
		stream << params_str << std::endl;
	}
}


void params_make_header(const std::filesystem::path& path) {
	ParamsGeneration::makeHeader(path, params);
}
void params_make_json(const std::filesystem::path& path) {
	ParamsGeneration::makeJson(path, params);
}