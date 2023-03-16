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

	void print2(
		heap_array<rr_float2, Params::maxn>&& r,
		heap_array<rr_int, Params::maxn>&& itype,
		std::optional<heap_array<rr_float2, Params::maxn>> v,
		std::optional<heap_array<rr_float, Params::maxn>> rho,
		std::optional<heap_array<rr_float, Params::maxn>> u,
		std::optional<heap_array<rr_float, Params::maxn>> p,
		const rr_uint itimestep)
	{
		try {
			std::ofstream stream(::dataOutputRelativePath + std::to_string(itimestep));

			stream << "fmt: ";
			auto add_format_specifier = [&](const char* value, auto& opt) {
				if (opt.has_value()) {
					stream << value << " ";
				}
			};
			add_format_specifier("vx vy", v);
			add_format_specifier("rho", rho);
			add_format_specifier("u", u);
			add_format_specifier("p", p);
			stream << std::endl;

			stream << Params::particles_total << std::endl;
			for (rr_uint i = 0; i < Params::particles_total; i++) {
				stream << r(i).x << std::endl << r(i).y << std::endl;
				stream << itype(i) << std::endl;

				if (v) {
					stream << v.value().at(i).x << std::endl << v.value().at(i).y << std::endl;
				}
				if (rho) {
					stream << rho.value().at(i) << std::endl;
				}
				if (u) {
					stream << u.value().at(i) << std::endl;
				}
				if (p) {
					stream << p.value().at(i) << std::endl;
				}
			}
		}
		catch (std::exception& ex) {
			printlog("print2 error :")(ex.what())();
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

void output2(
	heap_array<rr_float2, Params::maxn>&& r,
	heap_array<rr_int, Params::maxn>&& itype,
	std::optional<heap_array<rr_float2, Params::maxn>> v,
	std::optional<heap_array<rr_float, Params::maxn>> rho,
	std::optional<heap_array<rr_float, Params::maxn>> u,
	std::optional<heap_array<rr_float, Params::maxn>> p,
	const rr_uint itimestep)
{
	printlog_debug()(__func__)();

	std::thread(print2,
		std::move(r),
		std::move(itype),
		std::move(v),
		std::move(rho),
		std::move(u),
		std::move(p),
		itimestep
	).detach();

	std::cout << "output2: " << itimestep << std::endl;
}

void fast_output(
	const heap_array<rr_float2, Params::maxn>& r,	// coordinates of all particles
	const heap_array<rr_int, Params::maxn>& itype,	// material type 
	const rr_uint ntotal,	// number of particles
	const rr_uint itimestep)// current time step
{
	fast_output(r.copy(), itype, ntotal, itimestep);

	std::cout << "fast output: " << itimestep << std::endl;
}

void fast_output(
	heap_array<rr_float2, Params::maxn>&& r,	// coordinates of all particles
	const heap_array<rr_int, Params::maxn>& itype,	// material type 
	const rr_uint ntotal,	// number of particles
	const rr_uint itimestep) // current time step
{
	printlog_debug()(__func__)();
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

		std::cout << timeStep << " / " << Params::maxtimestep << "   (part: " << Params::particles_total << ")";
		std::cout << "{ passed: " << getTimeInAppropriateForm(timePassed);
		std::cout << "; w8 est." << getTimeInAppropriateForm(timeEstimates) << " }" << std::endl;
	}
}