#pragma once
  
// artificial equation of state for the artificial compressibility
rr_float p_art_water(rr_float rho);
rr_float rho_from_p_art_water(rr_float p);

// artificial sound velocity
rr_float c_art_water();
rr_float c_sqr_art_water();
