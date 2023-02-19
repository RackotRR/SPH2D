#include <iostream>
#include "Test.h"


#define test_module(testing_function)   \
    result = testing_function();        \
    std::cout << #testing_function << (result ? " passed" : " failed") << std::endl;

Test::Test() {
    bool result;
    //test_module(test_smoothing_kernel);
    //test_module(test_eos);
    //test_module(test_predict_step);
    //test_module(test_grid_find);
    //test_module(test_sum_density);
    //test_module(test_find_stress_tensor);
    //test_module(test_update_internal_state);
    //test_module(test_find_internal_changes_pidrho2i_pjdrho2j);
    //test_module(test_find_internal_changes_pij_d_rhoij);
    //test_module(test_external_force);
    //test_module(test_artificial_viscosity);
    //test_module(test_average_velocity);
    test_module(test_single_step);
}
