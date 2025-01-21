#pragma once
#include <filesystem>
#include "CommonIncl.h"

void fileInput(
	vheap_darray_floatn& r_var,	// coordinates of all particles
	vheap_darray_floatn& v_var,	// velocities of all particles
	heap_darray<rr_float>& rho,	// particle densities
	heap_darray<rr_float>& p,	// particle pressure
	heap_darray<rr_int>& itype,// particle material type 
	const std::filesystem::path& initial_dump_path,
	const std::filesystem::path& experiment_directory);

// loading or generating initial particle information
void cli(
	vheap_darray_floatn& r_var,	// coordinates of all particles
	vheap_darray_floatn& v_var,	// velocities of all particles
	heap_darray<rr_float>& rho,	// particle densities
	heap_darray<rr_float>& p,	// particle pressure
	heap_darray<rr_int>& itype);// particle material type 
