#pragma once
  
// artificial equation of state for the artificial compressibility
void p_art_water(
	const double rho, // density
	const double u,	// internal energy
	double& p,  // out, pressure
	double& c); // out, sound velocity