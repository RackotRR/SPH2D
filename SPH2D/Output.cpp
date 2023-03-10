#include <string>
#include <fstream>
#include <iostream>
#include <thread>
#include <filesystem>
#include <nlohmann/json.hpp>

#include "Output.h"
#include "Input.h"
#include "CLCommon.h"

namespace {
	std::string experimentRelativePath;
	std::string dataOutputRelativePath;

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
	void printFull(
		const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
		const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
		const heap_array<rr_float, Params::maxn>& rho,// density
		const heap_array<rr_float, Params::maxn>& p,	// pressure
		const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
		const heap_array<rr_float, Params::maxn>& c,	// sound velocity
		const heap_array<rr_int, Params::maxn>& itype,	// material type 
		const rr_uint ntotal,	// number of particles
		const rr_uint itimestep) // current time step
	{
		std::ofstream stream(::dataOutputRelativePath + std::to_string(itimestep));
		stream << ntotal << std::endl;
		for (rr_uint i = 0; i < ntotal; i++) {
			stream << r(i).x << std::endl << r(i).y << std::endl;
			stream << v(i).x << std::endl << v(i).y << std::endl;
			stream << rho(i) << std::endl;
			stream << p(i) << std::endl;
			stream << c(i) << std::endl;
			stream << itype(i) << std::endl;
		}
	}

	void print_on_demand(
		std::unique_ptr<heap_array<rr_float2, Params::maxn>> r,	// coordinates of all particles
		std::unique_ptr<heap_array<rr_int, Params::maxn>> itype,	// material type 
		std::unique_ptr<heap_array<rr_float2, Params::maxn>> v,	// velocities of all particles
		std::unique_ptr<heap_array<rr_float, Params::maxn>> rho,// density
		std::unique_ptr<heap_array<rr_float, Params::maxn>> p,	// pressure
		std::unique_ptr<heap_array<rr_float, Params::maxn>> u,	// specific internal energy
		const rr_uint ntotal,	// number of particles
		const rr_uint itimestep) // current time step
	{
		if (!r || !itype) {
			throw std::runtime_error{ "print_on_demand error: r and itype were expected not to be null" };
		}

		std::ofstream stream(::dataOutputRelativePath + std::to_string(itimestep));
		stream << ntotal << std::endl;
		for (rr_uint i = 0; i < ntotal; i++) {
			stream << r->at(i).x << std::endl << r->at(i).y << std::endl;
			stream << itype->at(i) << std::endl;

			if (v) {
				stream << v->at(i).x << std::endl << v->at(i).y << std::endl;
			}
			if (rho) {
				stream << rho->at(i) << std::endl;
			}
			if (p) {
				stream << p->at(i) << std::endl;
			}
			if (u) {
				stream << u->at(i) << std::endl;
			}
		}
	}
}

// params and other things
void printParams() {
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
		{"simulation_time", Params::simulation_time},
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

	makeParamsHeader(Params::particles_total,
		Params::particles_fluid,
		Params::particles_boundary,
		experimentRelativePath + "clparams.h");
}


void setupOutput() {
	::experimentRelativePath = Params::experimentName + "\\";
	::dataOutputRelativePath = ::experimentRelativePath + "data\\";
	auto analysisResultsPath = ::experimentRelativePath + "analysis\\";

	std::filesystem::create_directory(std::filesystem::current_path().append(::experimentRelativePath));
	std::filesystem::create_directory(std::filesystem::current_path().append(::dataOutputRelativePath));
	std::filesystem::create_directory(std::filesystem::current_path().append(analysisResultsPath));

	//init_logger();
	init_logger(Params::experimentName);
	logCLInfo();
}

void output_on_demand(
	std::unique_ptr<heap_array<rr_float2, Params::maxn>> r,	// coordinates of all particles
	std::unique_ptr<heap_array<rr_int, Params::maxn>> itype,	// material type 
	std::unique_ptr<heap_array<rr_float2, Params::maxn>> v,	// velocities of all particles
	std::unique_ptr<heap_array<rr_float, Params::maxn>> rho,// density
	std::unique_ptr<heap_array<rr_float, Params::maxn>> p,	// pressure
	std::unique_ptr<heap_array<rr_float, Params::maxn>> u,	// specific internal energy
	const rr_uint ntotal,	// number of particles
	const rr_uint itimestep)// current time step
{
	printlog()(__func__)();

	std::thread(print_on_demand,
		std::move(r),
		std::move(itype),
		std::move(v),
		std::move(rho),
		std::move(p),
		std::move(u),
		ntotal,
		itimestep).detach();
}

// save particle information to external disk file
void output(
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_float2, Params::maxn>& v,	// velocities of all particles
	const heap_array<rr_float, Params::maxn>& rho,// density
	const heap_array<rr_float, Params::maxn>& p,	// pressure
	const heap_array<rr_float, Params::maxn>& u,	// specific internal energy
	const heap_array<rr_float, Params::maxn>& c,	// sound velocity
	const heap_array<rr_int, Params::maxn>& itype,	// material type 
	const rr_uint ntotal,	// number of particles
	const rr_uint itimestep)// current time step
{
	printlog()(__func__)();

	std::thread(printFull, 
		r.copy(), 
		v.copy(),
		rho.copy(),
		p.copy(),
		u.copy(),
		c.copy(),
		itype.copy(), 
		ntotal, 
		itimestep).detach();
}

void fast_output(
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_int, Params::maxn>& itype,	// material type 
	const rr_uint ntotal,	// number of particles
	const rr_uint itimestep)// current time step
{
	fast_output(r.copy(), itype, ntotal, itimestep);
}

void fast_output(
	heap_array<rr_float2, Params::maxn>&& r,	// coordinates of all particles
	const heap_array<rr_int, Params::maxn>& itype,	// material type 
	const rr_uint ntotal,	// number of particles
	const rr_uint itimestep) // current time step
{
	printlog()(__func__)();
	std::thread(printFast, std::move(r), itype.copy(), ntotal, itimestep).detach();
}

static std::string getTimeInAppropriateForm(long long timeSec) {
	if (timeSec > 120) {
		timeSec /= 60;
	}
	else {
		return std::to_string(timeSec) + "s";
	}

	if (timeSec > 120) {
		timeSec /= 60;
	}
	else {
		return std::to_string(timeSec) + "m";
	}

	return std::to_string(timeSec) + "h";
}

void printTimeEstimate(long long totalTime_ns, rr_uint timeStep) {
	if (timeStep && timeStep % Params::print_time_est_step == 0) {
		constexpr rr_float coefNanosecondsToSeconds = 1.E-9f;
		long long timeEstimates = static_cast<long long>((totalTime_ns / timeStep) * coefNanosecondsToSeconds * (Params::maxtimestep - timeStep));
		long long timePassed = static_cast<long long>(totalTime_ns * coefNanosecondsToSeconds);

		std::cout << timeStep << " / " << Params::maxtimestep << " \t (part: " << Params::particles_total << ")";
		std::cout << "{ passed: " << getTimeInAppropriateForm(timePassed);
		std::cout << "; w8 est." << getTimeInAppropriateForm(timeEstimates) << " }" << std::endl;
	}
}