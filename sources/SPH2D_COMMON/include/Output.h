#pragma once
#include "CommonIncl.h"
#include <memory>
#include <optional>

void dump(
	heap_darray<rr_float2>&& r,
	heap_darray<rr_int>&& itype,
	heap_darray<rr_float2>&& v,
	heap_darray<rr_float>&& rho,
	heap_darray<rr_float>&& u,
	heap_darray<rr_float>&& p,
	const rr_uint itimestep);

void output(
	heap_darray<rr_float2>&& r,
	heap_darray<rr_int>&& itype,
	std::optional<heap_darray<rr_float2>> v,
	std::optional<heap_darray<rr_float>> rho,
	std::optional<heap_darray<rr_float>> u,
	std::optional<heap_darray<rr_float>> p,
	const rr_uint itimestep);

void fast_output(
	const heap_darray<rr_float2>& r, // coordinates of all particles
	const heap_darray<rr_int>& itype, // material type 
	const rr_uint ntotal,	// number of particles
	const rr_uint itimestep);// current time step
void fast_output(
	heap_darray<rr_float2>&& r,	// coordinates of all particles
	const heap_darray<rr_int>& itype,	// material type 
	const rr_uint ntotal,	// number of particles
	const rr_uint itimestep); // current time step

// call once at start
void setupOutput();

// call once at start after initialization
void printParams();

void printTimeEstimate(long long totalTime_ns, rr_uint timeStep);

void makeParamsHeader(std::string path = "cl/clparams.h");