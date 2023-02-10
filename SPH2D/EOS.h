#pragma once
  
// artificial equation of state for the artificial compressibility
void p_art_water(
	const rr_float rho, // density
	const rr_float u,	// internal energy
	rr_float& p,  // out, pressure
	rr_float& c); // out, sound velocity