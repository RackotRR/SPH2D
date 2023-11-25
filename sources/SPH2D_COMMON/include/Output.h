#pragma once
#include <filesystem>
#include <chrono>
#include <functional>

#include "CommonIncl.h"
#include <RR/Threads/ThreadPool.h>
#include <RR/Time/Timer.h>


class SPH2DOutput {
	using func_load_arr_float2 = std::function<shared_darray<rr_float2>()>;
	using func_load_arr_float = std::function<shared_darray<rr_float>()>;
	using func_load_arr_int = std::function<shared_darray<rr_int>()>;
public:
	static SPH2DOutput& instance() {
		static SPH2DOutput the_instance;
		return the_instance;
	}

	void initialize(const std::filesystem::path& experiment_path);
	void setup_output(
		func_load_arr_float2 load_r,
		func_load_arr_int load_itype,
		func_load_arr_float2 load_v,
		func_load_arr_float load_p,
		func_load_arr_float load_rho);

	void start_step();
	void finish_step();
	void update_step(rr_float time, rr_uint itimestep);

	void output(
		shared_darray<rr_float2> r,
		shared_darray<rr_int> itype,
		shared_darray<rr_float2> v,
		shared_darray<rr_float> rho,
		shared_darray<rr_float> p,
		rr_uint itimestep,
		rr_float time);
	void dump(
		shared_darray<rr_float2> r,
		shared_darray<rr_int> itype,
		shared_darray<rr_float2> v,
		shared_darray<rr_float> rho,
		shared_darray<rr_float> p,
		rr_uint itimestep,
		rr_float time);
	void crash_dump(
		shared_darray<rr_float2> r,
		shared_darray<rr_int> itype,
		shared_darray<rr_float2> v,
		shared_darray<rr_float> rho,
		shared_darray<rr_float> p,
		rr_uint itimestep,
		rr_float time);
private:
	std::vector<std::string> make_csv_header(
		bool r,
		bool itype,
		bool v,
		bool rho,
		bool p);
	void print_dump(
		shared_darray<rr_float2> r,
		shared_darray<rr_int> itype,
		shared_darray<rr_float2> v,
		shared_darray<rr_float> rho,
		shared_darray<rr_float> p,
		const std::filesystem::path& path);
	void print(
		shared_darray<rr_float2> r,
		shared_darray<rr_int> itype,
		shared_darray<rr_float2> v,
		shared_darray<rr_float> rho,
		shared_darray<rr_float> p,
		const std::filesystem::path& path);

	SPH2DOutput() = default;
	SPH2DOutput(const SPH2DOutput&) = delete;
	SPH2DOutput& operator=(const SPH2DOutput&) = delete;
private:
	RR::Threads::ThreadPool thread_pool;
	RR::Timer<std::chrono::nanoseconds, std::chrono::nanoseconds> timer;

	std::filesystem::path experiment_path;
	std::filesystem::path data_path;
	std::filesystem::path dump_path;
	std::filesystem::path analysis_path;

	rr_float last_save_time;
	rr_float last_dump_time;

	func_load_arr_float2 load_r;
	func_load_arr_int load_itype;
	func_load_arr_float2 load_v;
	func_load_arr_float load_p;
	func_load_arr_float load_rho;
};