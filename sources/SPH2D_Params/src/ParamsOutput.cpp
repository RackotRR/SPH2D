#include <nlohmann/json.hpp>
#include <fmt/format.h>
#include <sstream>
#include <fstream>

#include "Params.h"
#include "Version.h"


namespace {
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
			set(version_major);
			set(version_minor);
			set(dim);
			set(maxn);
			set(max_neighbours);
			set(max_cells);
			set(x_maxgeom);
			set(x_mingeom);
			set(y_maxgeom);
			set(y_mingeom);
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
			set(fluid_particles_per_d);
			set(wave_length);
			set(depth);
			set(freq);
			set(piston_amp);
			set(wave_amp);
			set(wave_number);
			set(beach_x);
			set(left_wall_start);
			set(left_wall_end);
			set(generator_time_wait);
			set(dt);
			set(simulation_time);
			set(local_threads);
			set(eos_csqr_k);
			set(pa_sph);
			set(skf);
			set(int_force_kernel);
			set(nwm);
			set(boundary_layers_num);
			set(hsml);
			set(delta);
			set(boundary_delta);
			set(summation_density);
			set(nor_density);
			set(average_velocity);
			set(average_velocity_epsilon);
			set(visc);
			set(water_dynamic_visc);
			set_param("TYPE_BOUNDARY", params.TYPE_BOUNDARY);
			set_param("TYPE_NON_EXISTENT", params.TYPE_NON_EXISTENT);
			set_param("TYPE_WATER", params.TYPE_WATER);
			set(mass);
			set(enable_check_consistency);
			set(inf_stop);
			set(maxtimestep);
			set(normal_check_step);
			set(save_step);
			set(dump_step);
			set(print_time_est_step);
			set(pi);
			set(g);
			set(experiment_name);
			set(format_line);
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
			buffer << fmt::format("#define params_{} {:.10f}f\n", param_name, value);
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
}

void ExperimentParams::makeHeader(const std::string& path) {
	::ParamsHeader header;
	header.serialize(*this);
	std::string params_str = header.string();
	std::ofstream stream{ path };
	stream << params_str << std::endl;
}
void ExperimentParams::makeJson(const std::string& path) {
	::ParamsJson json;
	json.serialize(*this);
	std::string params_str = json.string();
	std::ofstream stream{ path };
	stream << params_str << std::endl;
}