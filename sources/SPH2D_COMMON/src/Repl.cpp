#include <iostream>
#include <string>
#include <filesystem>

#include "Input.h"
#include "Output.h"

static std::vector<std::string> findTimeLayersPath(const std::filesystem::path& directory, int save_step, int start = 0) {
	std::vector<std::string> meta;
	int current_step = start;
	while (true) {
		auto path = directory / std::to_string(current_step);
		if (std::filesystem::exists(path)) {
			meta.emplace_back(path.filename().string());
			current_step += save_step;
		}
		else {
			break;
		}
	}
	return meta;
}

bool overrideDirectory() {
	while (true) {
		std::cout << "Would you like to write here?: " << params.experiment_name << std::endl;
		std::cout << "[y/n]" << std::endl;
		std::string answer;
		std::cin >> answer;
		if (answer == "y" || answer == "yes") {
			std::filesystem::remove_all(std::filesystem::current_path().append(params.experiment_name));
			return true;
		}
		else if (answer == "n" || answer == "no") {
			return false;
		}
	}
}

// loading or generating initial particle information
void repl(
	heap_darray<rr_float2>& r,	// coordinates of all particles
	heap_darray<rr_float2>& v,	// velocities of all particles
	heap_darray<rr_float>& mass,	// particle masses
	heap_darray<rr_float>& rho,	// particle densities
	heap_darray<rr_float>& p,	// particle pressure
	heap_darray<rr_float>& u,	// particle internal energy
	heap_darray<rr_float>& c,	// sound velocity
	heap_darray<rr_int>& itype,	// particle material type 
	rr_uint& ntotal, // total particle number
	rr_uint& nfluid) // total fluid particles
{
	while (true) {
		try {
			std::cout << "[SPH2D] Experiment name: ";
			std::string line;
			std::getline(std::cin, params.experiment_name);

			auto experiment_directory = std::filesystem::current_path() / params.experiment_name;
			if (std::filesystem::exists(experiment_directory)) {

				auto params_path = experiment_directory / "Params.json";
				if (!std::filesystem::exists(params_path)) {
					std::cout << "Can't find Params.json" << std::endl;
					if (overrideDirectory()) {
						input(r, v, mass, rho, p, u, c, itype, ntotal, nfluid);
						return;
					}
					else {
						continue;
					}
				} // no Params.json

				std::cout << "Found experiment" << std::endl;
				params.load(params_path.string());
				auto data_path = experiment_directory / "data";
				auto dump_path = experiment_directory / "dump";
				auto time_layers_path = findTimeLayersPath(data_path, params.save_step);
				auto dumps_path = findTimeLayersPath(dump_path, params.dump_step, params.dump_step);
				std::cout << "found " << time_layers_path.size() << " time layers" << std::endl;
				std::cout << "found " << dumps_path.size() << " dumps:" << std::endl;
				for (int i = 0; i < dumps_path.size(); ++i) {
					std::cout << "[" << i + 1 << "]: " << dumps_path[i] << std::endl;
				}

				std::cout << "Write dump number to load or [0] to start from Params.json" << std::endl;
				std::cout << "Write [-1] to clear directory and start from generated Params" << std::endl;
				int num;
				std::cin >> num;
				std::cin.get();

				if (num == -1 && overrideDirectory()) {
					input(r, v, mass, rho, p, u, c, itype, ntotal, nfluid);
					return;
				} // clear and run default
				else if (num == 0 && overrideDirectory()) {
					input(r, v, mass, rho, p, u, c, itype, ntotal, nfluid, false);
					return;
				} // clear and run from Params.json
				else if (num <= dumps_path.size()) {
					auto& chosen_dump = dumps_path[num - 1];
					auto iter = std::find(time_layers_path.begin(), time_layers_path.end(), chosen_dump);
					for (iter++; iter != time_layers_path.end(); ++iter) {
						auto path_to_remove = experiment_directory / "data" / *iter;
						std::filesystem::remove(path_to_remove);
						std::cout << "layer " << *iter << " removed" << std::endl;
					}
					std::cout << "loading data... ";

					auto particles_data_path = std::filesystem::path(experiment_directory) / "dump" / chosen_dump;
					fileInput(r, v, mass, rho, p, u, c, itype, ntotal, nfluid, particles_data_path.string());
					std::cout << "completed" << std::endl;
					return;
				} // run from dump
				else {
					std::cout << "Wrong number provided!" << std::endl;
				} // wrong number
			}
			else {
				input(r, v, mass, rho, p, u, c, itype, ntotal, nfluid);
				return;
			}
		}
		catch (const std::exception& ex) {
			std::cout << "Something went wrong: " << std::endl << ex.what() << std::endl;
		}
	}
}