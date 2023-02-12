#include <string>
#include <fstream>
#include <iostream>
#include <thread>
#include <filesystem>
#include <nlohmann/json.hpp>

#include "Output.h"

namespace {
	std::string experimentRelativePath;
	std::string dataOutputRelativePath;

	// params and other things
	void printParams() {
		printlog()(__func__)();
		
		nlohmann::json json;
		json = {
			{"experiment_name", Params::experimentName},
			{"x_mingeom", Params::x_mingeom},
			{"y_mingeom", Params::y_mingeom},
			{"x_size", Params::x_maxgeom - Params::x_mingeom},
			{"y_size", Params::y_maxgeom - Params::y_mingeom},
			{"dx", Params::delta},
			{"dy", Params::delta},
			{"dt", Params::dt},
			{"hsml", Params::hsml},
			{"simulation_time", Params::simulationTime},
			{"save_step", Params::save_step},
			{"wave_length", Params::L },
			{"depth", Params::d},
			{"skf", Params::skf},
			{"particles_per_d", Params::fluid_particles_per_d},
			{"particles_fluid", Params::particles_fluid},
			{"particles_boundary", Params::particles_boundary},
			{"particles_total", Params::particles_total}
		};

		std::string path = experimentRelativePath + "Params.json";
		std::ofstream stream(path, std::ofstream::out);
		stream << json;
	}

	void printFast(
		const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
		const heap_array<rr_int, Params::maxn> itype,
		const rr_uint ntotal,
		const rr_uint itimestep)
	{
		std::ofstream stream(::dataOutputRelativePath + std::to_string(itimestep));
		stream << ntotal << std::endl;
		for (rr_uint i = 0; i < ntotal; i++) {
			stream << r(i).x << std::endl << r(i).y << std::endl;
			stream << itype(i) << std::endl;
		}
	}
}


void setupOutput() {
	printlog(__func__)();

	::experimentRelativePath = Params::experimentName + "\\";
	::dataOutputRelativePath = ::experimentRelativePath + "data\\";
	auto analysisResultsPath = ::experimentRelativePath + "analysis\\";

	std::filesystem::create_directory(std::filesystem::current_path().append(::experimentRelativePath));
	std::filesystem::create_directory(std::filesystem::current_path().append(::dataOutputRelativePath));
	std::filesystem::create_directory(std::filesystem::current_path().append(analysisResultsPath));

	printParams();
}

// save particle information to external disk file
void output(
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& mass,// particle masses
	const heap_array<rr_float, Params::maxn>& rho,// density
	const heap_array<rr_float, Params::maxn>& p,	// pressure
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	const heap_array<rr_float, Params::maxn>& c,	// sound velocity
	const heap_array<rr_int, Params::maxn>& itype,	// material type 
	const rr_uint ntotal,	// number of particles
	const rr_uint itimestep,// current time step
	const long long timePassedTotal,
	const long long timeEstimates)
{
	printlog()(__func__)();

	std::cout << itimestep << " / " << Params::maxtimestep << " \t (part: " << ntotal << ")";
	std::cout << "{ passed: " << timePassedTotal << "; w8 est." << timeEstimates << " }" << std::endl;

	std::thread(printFast, r.copy(), itype.copy(), ntotal, itimestep).detach();
}