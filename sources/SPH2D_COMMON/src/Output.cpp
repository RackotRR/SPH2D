#include <string>
#include <fstream>
#include <iostream>
#include <thread>
#include <filesystem>

#include "Output.h"
#include "Input.h"


namespace {
	std::string experimentRelativePath;
	std::string dataOutputRelativePath;
	std::string dumpRelativePath;

	void printFast(
		const heap_darray<rr_float2>& r,	// coordinates of all particles
		const heap_darray<rr_int> itype,
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

	void printDump(
		heap_darray<rr_float2>&& r,
		heap_darray<rr_int>&& itype,
		heap_darray<rr_float2>&& v,
		heap_darray<rr_float>&& rho,
		heap_darray<rr_float>&& u,
		heap_darray<rr_float>&& p,
		const rr_uint itimestep)
	{
		std::ofstream stream(::dumpRelativePath + std::to_string(itimestep));
		stream << params.ntotal << std::endl;
		for (rr_uint i = 0; i < params.ntotal; i++) {
			stream << r(i).x << std::endl << r(i).y << std::endl;
			stream << itype(i) << std::endl;
			stream << v(i).x << std::endl << v(i).y << std::endl;
			stream << rho(i) << std::endl;
			stream << u(i) << std::endl;
			stream << p(i) << std::endl;
		}
	}

	void print(
		heap_darray<rr_float2>&& r,
		heap_darray<rr_int>&& itype,
		std::optional<heap_darray<rr_float2>> v,
		std::optional<heap_darray<rr_float>> rho,
		std::optional<heap_darray<rr_float>> u,
		std::optional<heap_darray<rr_float>> p,
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

			stream << params.ntotal << std::endl;
			for (rr_uint i = 0; i < params.ntotal; i++) {
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
			printlog("print error :")(ex.what())();
		}
	}
}


void setupOutput() {
	::experimentRelativePath = params.experiment_name + "/";
	::dataOutputRelativePath = ::experimentRelativePath + "data/";
	::dumpRelativePath = ::experimentRelativePath + "dump/";
	auto analysisResultsPath = ::experimentRelativePath + "analysis/";

	std::filesystem::create_directory(std::filesystem::current_path().append(::experimentRelativePath));
	std::filesystem::create_directory(std::filesystem::current_path().append(::dataOutputRelativePath));
	std::filesystem::create_directory(std::filesystem::current_path().append(::dumpRelativePath));
	std::filesystem::create_directory(std::filesystem::current_path().append(analysisResultsPath));

	//init_logger();
	init_logger(params.experiment_name + "/log.txt");
}

void dump(
	heap_darray<rr_float2>&& r,
	heap_darray<rr_int>&& itype,
	heap_darray<rr_float2>&& v,
	heap_darray<rr_float>&& rho,
	heap_darray<rr_float>&& u,
	heap_darray<rr_float>&& p,
	const rr_uint itimestep) 
{
	printlog_debug()(__func__)();

	std::thread(printDump,
		std::move(r),
		std::move(itype),
		std::move(v),
		std::move(rho),
		std::move(u),
		std::move(p),
		itimestep
	).detach();

	std::cout << "dump: " << itimestep << std::endl;
}

void output(
	heap_darray<rr_float2>&& r,
	heap_darray<rr_int>&& itype,
	std::optional<heap_darray<rr_float2>> v,
	std::optional<heap_darray<rr_float>> rho,
	std::optional<heap_darray<rr_float>> u,
	std::optional<heap_darray<rr_float>> p,
	const rr_uint itimestep)
{
	printlog_debug()(__func__)();

	std::thread(print,
		std::move(r),
		std::move(itype),
		std::move(v),
		std::move(rho),
		std::move(u),
		std::move(p),
		itimestep
	).detach();

	std::cout << "output: " << itimestep << std::endl;
}

void fast_output(
	const heap_darray<rr_float2>& r,	// coordinates of all particles
	const heap_darray<rr_int>& itype,	// material type 
	const rr_uint ntotal,	// number of particles
	const rr_uint itimestep)// current time step
{
	fast_output(r.copy(), itype, ntotal, itimestep);

	std::cout << "fast output: " << itimestep << std::endl;
}

void fast_output(
	heap_darray<rr_float2>&& r,	// coordinates of all particles
	const heap_darray<rr_int>& itype,	// material type 
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
	rr_uint time_steps_passed = timeStep - params.starttimestep;
	if (time_steps_passed && timeStep % params.print_time_est_step == 0) {
		constexpr rr_float coefNanosecondsToSeconds = 1.E-9f;
		long long timeEstimates = static_cast<long long>(
			(totalTime_ns / time_steps_passed) * coefNanosecondsToSeconds * (params.maxtimestep - timeStep)
			);
		long long timePassed = static_cast<long long>(totalTime_ns * coefNanosecondsToSeconds);

		std::cout << timeStep << " / " << params.maxtimestep << "   (part: " << params.ntotal << ")";
		std::cout << "{ passed: " << getTimeInAppropriateForm(timePassed);
		std::cout << "; w8 est." << getTimeInAppropriateForm(timeEstimates) << " }" << std::endl;
	}
}

// params and other things
void printParams() {
	params.makeJson(experimentRelativePath + "Params.json");
	params.makeHeader("cl/clparams.h");
}