#pragma once
#include <filesystem>
#include <chrono>
#include <functional>

#include "CommonIncl.h"
#include <RR/Threads/ThreadPool.h>
#include <RR/Time/Timer.h>


class RRSPHOutput {
	using func_load_arr_floatn = std::function<shared_vheap_darray_floatn()>;
	using func_load_arr_float = std::function<shared_darray<rr_float>()>;
	using func_load_arr_int = std::function<shared_darray<rr_int>()>;
public:
	static RRSPHOutput& instance() {
		static RRSPHOutput the_instance;
		return the_instance;
	}

	void initialize(const std::filesystem::path& experiment_path);
	void setup_output(
		func_load_arr_floatn load_r,
		func_load_arr_int load_itype,
		func_load_arr_floatn load_v,
		func_load_arr_float load_p,
		func_load_arr_float load_rho);

	void start_step(rr_float time);
	void finish_step();
	void update_step(rr_float time, rr_uint itimestep);

	void output(
		shared_vheap_darray_floatn r,
		shared_darray<rr_int> itype,
		shared_vheap_darray_floatn v,
		shared_darray<rr_float> rho,
		shared_darray<rr_float> p,
		rr_uint itimestep,
		rr_float time);
	void dump(
		shared_vheap_darray_floatn r,
		shared_darray<rr_int> itype,
		shared_vheap_darray_floatn v,
		shared_darray<rr_float> rho,
		shared_darray<rr_float> p,
		rr_uint itimestep,
		rr_float time);
	void crash_dump(
		shared_vheap_darray_floatn r,
		shared_darray<rr_int> itype,
		shared_vheap_darray_floatn v,
		shared_darray<rr_float> rho,
		shared_darray<rr_float> p,
		rr_uint itimestep,
		rr_float time);
private:
	RRSPHOutput() = default;
	RRSPHOutput(const RRSPHOutput&) = delete;
	RRSPHOutput& operator=(const RRSPHOutput&) = delete;
private:
	RR::Threads::ThreadPool thread_pool;
	RR::Timer<std::chrono::nanoseconds, std::chrono::nanoseconds> timer;

	std::filesystem::path experiment_path;
	std::filesystem::path data_path;
	std::filesystem::path dump_path;
	std::filesystem::path analysis_path;

	rr_float last_save_time{};
	rr_float last_dump_time{};

	func_load_arr_floatn load_r;
	func_load_arr_int load_itype;
	func_load_arr_floatn load_v;
	func_load_arr_float load_p;
	func_load_arr_float load_rho;

	int save_time_decimal_digits{};
};