#include <fstream>
#include <iostream>
#include <thread>
#include <csv-parser/csv.hpp>
#include <array>
#include <fmt/format.h>

#include "Output.h"
#include "ConsistencyCheck.h"
#include "TimeEstimate.h"
#include "TimeFormat.h"

template<typename T>
static void check_passed_consistency(const shared_darray<T>& arr, bool need) {
	assert(need && arr.get() && arr->size() == params.maxn);
}

static void check_passed_consistency(const shared_vheap_darray_floatn& arr, bool need) {
	assert(need && arr.get());
	if (params.dim == 3) {
		assert(arr->get_flt3().size() == params.maxn);
	}
	else {
		assert(arr->get_flt2().size() == params.maxn);
	}
}

void RRSPHOutput::initialize(const std::filesystem::path& experiment_path) {
	this->experiment_path = experiment_path;
	if (!std::filesystem::exists(experiment_path)) {
		throw std::runtime_error{ "No experiment directory exists: " + experiment_path.string() };
	}

	this->data_path = experiment_path / "data";
	this->dump_path = experiment_path / "dump";
	this->analysis_path = experiment_path / "analysis";
	std::filesystem::create_directory(data_path);
	std::filesystem::create_directory(dump_path);
	std::filesystem::create_directory(analysis_path);

	//init_logger();
	init_logger((experiment_path / "log.txt").string());
}

void RRSPHOutput::setup_output(
	func_load_arr_floatn load_r,
	func_load_arr_int load_itype,
	func_load_arr_floatn load_v,
	func_load_arr_float load_p,
	func_load_arr_float load_rho)
{
	printlog(__func__)();

	this->load_r = load_r;
	this->load_itype = load_itype;
	this->load_v = load_v;
	this->load_p = load_p;
	this->load_rho = load_rho;
	this->save_time_decimal_digits = format_count_decimal_digits(params.save_time);

	this->last_save_time = params.start_simulation_time;
	this->last_dump_time = params.start_simulation_time;

	// first layer output
	auto r = load_r();
	auto itype = load_itype();
	shared_vheap_darray_floatn v = params.save_velocity ? load_v() : nullptr;
	shared_darray<rr_float> rho = params.save_density ? load_rho() : nullptr;
	shared_darray<rr_float> p = params.save_pressure ? load_p() : nullptr;

	output(
		r,
		itype,
		v,
		rho,
		p,
		0,
		params.start_simulation_time);
}

static std::vector<std::string> make_csv_header(
	bool r,
	bool itype,
	bool v,
	bool rho,
	bool p)
{
	std::vector<std::string> header;

	if (r) {
		header.emplace_back(NAME_VARIABLE_X);
		header.emplace_back(NAME_VARIABLE_Y);
		if (params.dim == 3) {
			header.emplace_back(NAME_VARIABLE_Z);
		}
	}
	if (itype) {
		header.emplace_back(NAME_VARIABLE_ITYPE);
	}
	if (v) {
		header.emplace_back(NAME_VARIABLE_VX);
		header.emplace_back(NAME_VARIABLE_VY);
		if (params.dim == 3) {
			header.emplace_back(NAME_VARIABLE_VZ);
		}
	}
	if (rho) {
		header.emplace_back(NAME_VARIABLE_RHO);
	}
	if (p) {
		header.emplace_back(NAME_VARIABLE_P);
	}

	return header;
}

template<typename rr_floatn>
void print_dump(
	shared_vheap_darray_floatn r_var,
	shared_darray<rr_int> itype,
	shared_vheap_darray_floatn v_var,
	shared_darray<rr_float> rho,
	shared_darray<rr_float> p,
	const std::filesystem::path& path)
{
	std::ofstream stream{ path };
	auto writer = csv::make_csv_writer(stream);
	writer << make_csv_header(r_var.get(), itype.get(), v_var.get(), rho.get(), p.get());

	const heap_darray<rr_floatn>* r{};
	const heap_darray<rr_floatn>* v{};
	if constexpr (is_using_float3<rr_floatn>()) {
		r = &r_var->get_flt3();
		v = &v_var->get_flt3();
	}
	else {
		r = &r_var->get_flt2();
		v = &v_var->get_flt2();
	}

	for (rr_uint i = 0; i < params.ntotal; i++) {
		if constexpr (is_using_float3<rr_floatn>()) {
			writer << std::array{
				std::to_string(r->at(i).x),
				std::to_string(r->at(i).y),
				std::to_string(r->at(i).z),
				std::to_string(itype->at(i)),
				std::to_string(v->at(i).x),
				std::to_string(v->at(i).y),
				std::to_string(v->at(i).z),
				std::to_string(rho->at(i)),
				std::to_string(p->at(i))
			};
		}
		else {
			writer << std::array{
				std::to_string(r->at(i).x),
				std::to_string(r->at(i).y),
				std::to_string(itype->at(i)),
				std::to_string(v->at(i).x),
				std::to_string(v->at(i).y),
				std::to_string(rho->at(i)),
				std::to_string(p->at(i))
			};
		}
	}
}

template<typename rr_floatn>
void print(
	shared_vheap_darray_floatn r_var,
	shared_darray<rr_int> itype,
	shared_vheap_darray_floatn v_var,
	shared_darray<rr_float> rho,
	shared_darray<rr_float> p,
	const std::filesystem::path& path)
{
	std::ofstream stream{ path };
	auto writer = csv::make_csv_writer(stream);
	auto header = make_csv_header(r_var.get(), itype.get(), v_var.get(), rho.get(), p.get());
	writer << header;

	const heap_darray<rr_floatn>* r{};
	const heap_darray<rr_floatn>* v{};
	if constexpr (is_using_float3<rr_floatn>()) {
		r = &r_var->get_flt3();
		v = &v_var->get_flt3();
	}
	else {
		r = &r_var->get_flt2();
		v = &v_var->get_flt2();
	}

	std::vector<std::string> row;
	row.reserve(header.size());
	for (rr_uint i = 0; i < params.ntotal; i++) {
		row.push_back(std::to_string(r->at(i).x));
		row.push_back(std::to_string(r->at(i).y));
		if constexpr (is_using_float3<rr_floatn>()) {
			row.push_back(std::to_string(r->at(i).z));
		}

		row.push_back(std::to_string(itype->at(i)));
		if (v) {
			row.push_back(std::to_string(v->at(i).x));
			row.push_back(std::to_string(v->at(i).y));
			if constexpr (is_using_float3<rr_floatn>()) {
				row.push_back(std::to_string(v->at(i).z));
			}
		}
		if (rho) {
			row.push_back(std::to_string(rho->at(i)));
		}
		if (p) {
			row.push_back(std::to_string(p->at(i)));
		}

		writer << row;
		row.clear();
	}
}

void RRSPHOutput::dump(
	shared_vheap_darray_floatn r_var,
	shared_darray<rr_int> itype,
	shared_vheap_darray_floatn v_var,
	shared_darray<rr_float> rho,
	shared_darray<rr_float> p,
	rr_uint itimestep,
	rr_float time)
{
	printlog_debug()(__func__)();

	check_passed_consistency(r_var, true);
	check_passed_consistency(itype, true);
	check_passed_consistency(v_var, true);
	check_passed_consistency(rho, true);
	check_passed_consistency(p, true);

	std::string t = format_save_time(time, params.dump_time);
	auto target_path = dump_path / fmt::format("{}.csv", t);

	if (params.dim == 3) {
		thread_pool.add_thread(
			std::thread(print_dump<rr_float3>,
				r_var,
				itype,
				v_var,
				rho,
				p,
				target_path));
	}
	else {
		thread_pool.add_thread(
			std::thread(print_dump<rr_float2>,
				r_var,
				itype,
				v_var,
				rho,
				p,
				target_path));
	}

	std::cout << fmt::format("dump: {} ({}s)", itimestep, t) << std::endl;
}

void RRSPHOutput::crash_dump(
	shared_vheap_darray_floatn r_var,
	shared_darray<rr_int> itype,
	shared_vheap_darray_floatn v_var,
	shared_darray<rr_float> rho,
	shared_darray<rr_float> p,
	rr_uint itimestep,
	rr_float time)
{
	printlog_debug()(__func__)();

	check_passed_consistency(r_var, true);
	check_passed_consistency(itype, true);
	check_passed_consistency(v_var, true);
	check_passed_consistency(rho, true);
	check_passed_consistency(p, true);

	std::string t = format_time_digits(time, save_time_decimal_digits);
	auto target_path = experiment_path / fmt::format("crash_dump_{}.csv", t);

	if (params.dim == 3) {
		thread_pool.add_thread(
			std::thread(print_dump<rr_float3>,
				r_var,
				itype,
				v_var,
				rho,
				p,
				target_path));
	}
	else {
		thread_pool.add_thread(
			std::thread(print_dump<rr_float2>,
				r_var,
				itype,
				v_var,
				rho,
				p,
				target_path));
	}

	std::cout << fmt::format("crash dump: {} ({}s)", itimestep, t) << std::endl;
}

void RRSPHOutput::output(
	shared_vheap_darray_floatn r_var,
	shared_darray<rr_int> itype,
	shared_vheap_darray_floatn v_var,
	shared_darray<rr_float> rho,
	shared_darray<rr_float> p,
	rr_uint itimestep,
	rr_float time)
{
	printlog_debug()(__func__)();

	check_passed_consistency(r_var, true);
	check_passed_consistency(itype, true);
	check_passed_consistency(v_var, params.save_velocity);
	check_passed_consistency(rho, params.save_density);
	check_passed_consistency(p, params.save_pressure);

	std::string t = format_time_digits(time, save_time_decimal_digits);
	auto target_path = data_path / fmt::format("{}.csv", t);

	if (params.dim == 3) {
		thread_pool.add_thread(
			std::thread(print<rr_float3>,
				r_var,
				itype,
				v_var,
				rho,
				p,
				target_path));
	}
	else {
		thread_pool.add_thread(
			std::thread(print<rr_float2>,
				r_var,
				itype,
				v_var,
				rho,
				p,
				target_path));
	}

	std::cout << fmt::format("output: {} ({}s)", itimestep, t) << std::endl;
}

void RRSPHOutput::start_step(rr_float time) {
	printlog()(fmt::format("time: {}/{} s", 
		format_save_time(time, params.save_time), 
		params.simulation_time))();

	timer.start();
}
void RRSPHOutput::finish_step() {
	timer.finish();
}

void RRSPHOutput::update_step(rr_float time, rr_uint itimestep) {
	printlog_debug()(fmt::format("{} at {} ({}s)", __func__, itimestep, time))();
	
	assert(load_r);
	assert(load_itype);
	assert(load_v);
	assert(load_p);
	assert(load_rho);

	bool should_save = (time - last_save_time) > params.save_time;
	if (params.save_every_step) should_save = true;

	bool should_dump = params.use_dump && (time - last_dump_time) > params.dump_time;
	bool should_check = params.consistency_check;

	if (should_save || should_dump || should_check) {

		bool need_r = should_save || should_dump || should_check;
		bool need_itype = should_save || should_dump || should_check;
		bool need_v = (should_save && params.save_velocity) || should_dump;
		bool need_p = (should_save && params.save_pressure) || should_dump;
		bool need_rho = (should_save && params.save_density) || should_dump;

		shared_vheap_darray_floatn r_temp = need_r ? load_r() : nullptr;
		shared_darray<rr_int> itype_temp = need_itype ? load_itype() : nullptr;
		shared_vheap_darray_floatn v_temp = need_v ? load_v() : nullptr;
		shared_darray<rr_float> p_temp = need_p ? load_p() : nullptr;
		shared_darray<rr_float> rho_temp = need_rho ? load_rho() : nullptr;

		if (should_check) {
			try {
				check_particles_are_within_boundaries(r_temp, itype_temp, params.consistency_treatment);
			}
			catch (...) {
				if (params.use_crash_dump) {
					if (!v_temp) v_temp = load_v();
					if (!p_temp) p_temp = load_p();
					if (!rho_temp) rho_temp = load_rho();

					crash_dump(
						r_temp,
						itype_temp,
						v_temp,
						rho_temp,
						p_temp,
						itimestep,
						time);
				}
				throw;
			}
		}

		if (should_save) {
			output(
				r_temp,
				itype_temp,
				params.save_velocity ? v_temp : nullptr,
				params.save_density ? rho_temp : nullptr,
				params.save_pressure ? p_temp : nullptr,
				itimestep,
				time);
			last_save_time = time;
		}

		if (should_dump) {
			dump(
				r_temp,
				itype_temp,
				v_temp,
				rho_temp,
				p_temp,
				itimestep,
				time);
			last_dump_time = time;
		}
	}

	bool should_estimate_by_default = !params.use_custom_time_estimate_step && should_save;
	if (check_custom_time_estimate_step(itimestep) || should_estimate_by_default) {
		print_time_estimate(itimestep, timer.total(), time);
	}
}