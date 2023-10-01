#include <string>
#include <fstream>
#include <iostream>
#include <thread>
#include <filesystem>
#include <csv-parser/csv.hpp>
#include <array>
#include <fmt/format.h>

#include <RR/Threads/ThreadPool.h>

#include "Output.h"
#include "Input.h"
#include "ParamsIO.h"

static std::filesystem::path experimentPath;
static std::filesystem::path dataPath;
static std::filesystem::path dumpPath;

namespace {
	RR::Threads::ThreadPool threadPool;

	void printFast(
		const heap_darray<rr_float2>& r,
		const heap_darray<rr_int>& itype,
		const rr_uint ntotal,
		const rr_uint itimestep)
	{
		std::ofstream stream{ dataPath / fmt::format("{}.csv", itimestep) };
		auto writer = csv::make_csv_writer(stream);
		writer << std::vector{ "x", "y", "itype" };
		for (rr_uint i = 0; i < ntotal; i++) {
			writer << std::array{
				std::to_string(r(i).x),
				std::to_string(r(i).y),
				std::to_string(itype(i))
			};
		}
	}

	void printDump(
		heap_darray<rr_float2>&& r,
		heap_darray<rr_int>&& itype,
		heap_darray<rr_float2>&& v,
		heap_darray<rr_float>&& rho,
		heap_darray<rr_float>&& p,
		const rr_uint itimestep)
	{
		std::ofstream stream{ dumpPath / fmt::format("{}.csv", itimestep) };
		auto writer = csv::make_csv_writer(stream);
		auto header = std::array{
			"x", "y", "itype", "vx", "vy", "rho", "p"
		};

		writer << header;
		for (rr_uint i = 0; i < params.ntotal; i++) {
			writer << std::array{
				std::to_string(r(i).x),
				std::to_string(r(i).y),
				std::to_string(itype(i)),
				std::to_string(v(i).x),
				std::to_string(v(i).y),
				std::to_string(rho(i)),
				std::to_string(p(i))
			};
		}
	}

	void print(
		heap_darray<rr_float2>&& r,
		heap_darray<rr_int>&& itype,
		std::optional<heap_darray<rr_float2>> v,
		std::optional<heap_darray<rr_float>> rho,
		std::optional<heap_darray<rr_float>> p,
		const std::filesystem::path& path)
	{
		try {
			std::ofstream stream{ path };
			auto writer = csv::make_csv_writer(stream);

			std::vector<std::string> header{
				"x", "y", "itype"
			};
			if (v) {
				header.emplace_back("vx");
				header.emplace_back("vy");
			}
			if (rho) {
				header.emplace_back("rho");
			}
			if (p) {
				header.emplace_back("p");
			}
			writer << header;

			std::vector<std::string> row;
			row.reserve(header.size());
			for (rr_uint i = 0; i < params.ntotal; i++) {
				row.push_back(std::to_string(r(i).x));
				row.push_back(std::to_string(r(i).y));
				row.push_back(std::to_string(itype(i)));
				if (v.has_value()) {
					row.push_back(std::to_string(v.value().at(i).x));
					row.push_back(std::to_string(v.value().at(i).y));
				}
				if (rho.has_value()) {
					row.push_back(std::to_string(rho.value().at(i)));
				}
				if (p.has_value()) {
					row.push_back(std::to_string(p.value().at(i)));
				}

				writer << row;
				row.clear();
			}
		}
		catch (std::exception& ex) {
			printlog("print error: ")(ex.what())();
		}
	}
}


void setupOutput(const std::filesystem::path& experiment_path) {
	experimentPath = experiment_path;
	dataPath = experimentPath / "data";
	dumpPath = experimentPath / "dump";
	auto analysisResultsPath = ::experimentPath / "analysis";

	std::filesystem::create_directory(experimentPath);
	std::filesystem::create_directory(dataPath);
	std::filesystem::create_directory(dumpPath);
	std::filesystem::create_directory(analysisResultsPath);

	//init_logger();
	init_logger((experimentPath / "log.txt").string());
}

void dump(
	heap_darray<rr_float2>&& r,
	heap_darray<rr_int>&& itype,
	heap_darray<rr_float2>&& v,
	heap_darray<rr_float>&& rho,
	heap_darray<rr_float>&& p,
	const rr_uint itimestep) 
{
	printlog_debug()(__func__)();

	threadPool.add_thread(
		std::thread(printDump,
			std::move(r),
			std::move(itype),
			std::move(v),
			std::move(rho),
			std::move(p),
			itimestep));

	std::cout << "dump: " << itimestep << std::endl;
}

void crash_dump(
	heap_darray<rr_float2>&& r,
	heap_darray<rr_int>&& itype,
	std::optional<heap_darray<rr_float2>> v,
	std::optional<heap_darray<rr_float>> rho,
	std::optional<heap_darray<rr_float>> p,
	const rr_uint itimestep)
{
	printlog_debug()(__func__)();

	threadPool.add_thread(
		std::thread(print,
			std::move(r),
			std::move(itype),
			std::move(v),
			std::move(rho),
			std::move(p),
			experimentPath / fmt::format("crash_dump_{}.csv", itimestep)));

	std::cout << "crash dump: " << itimestep << std::endl;
}

void output(
	heap_darray<rr_float2>&& r,
	heap_darray<rr_int>&& itype,
	std::optional<heap_darray<rr_float2>> v,
	std::optional<heap_darray<rr_float>> rho,
	std::optional<heap_darray<rr_float>> p,
	const rr_uint itimestep)
{
	printlog_debug()(__func__)();

	threadPool.add_thread(
		std::thread(print,
			std::move(r),
			std::move(itype),
			std::move(v),
			std::move(rho),
			std::move(p),
			dataPath / fmt::format("{}.csv", itimestep)));

	std::cout << "output: " << itimestep << std::endl;
}

void fast_output(
	heap_darray<rr_float2>&& r,	// coordinates of all particles
	const heap_darray<rr_int>& itype,	// material type 
	const rr_uint ntotal,	// number of particles
	const rr_uint itimestep) // current time step
{
	printlog_debug()(__func__)();
	threadPool.add_thread(
		std::thread(
			printFast,
			std::move(r),
			itype.copy(),
			ntotal,
			itimestep));
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
	rr_uint step_time_estimate = params.use_custom_time_estimate_step ?
		params.step_time_estimate : params.save_step;

	if (time_steps_passed && timeStep % step_time_estimate == 0) {
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
	params_make_header(std::filesystem::current_path() / "cl" / "clparams.h");
}