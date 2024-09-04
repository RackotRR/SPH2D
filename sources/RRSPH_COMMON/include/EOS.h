#pragma once
  
// artificial equation of state for the artificial compressibility
rr_float eos_art_p(rr_float rho);
rr_float eos_art_rho(rr_float p);

// artificial sound velocity
rr_float eos_art_c();
rr_float eos_art_c_sqr();
