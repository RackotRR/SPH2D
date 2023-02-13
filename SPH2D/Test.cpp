#include <iostream>
#include "Test.h"


#define test_module(testing_function)   \
    result = testing_function();        \
    std::cout << #testing_function << (result ? " passed" : " failed") << std::endl;

Test::Test() {
    bool result;
    test_module(test_smoothing_kernel);
    test_module(test_eos);
}
