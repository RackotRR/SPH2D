#pragma once
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <nlohmann/json.hpp>

struct HeightTestingParams {
    std::string mode;

    double t;
    double x0 = 0;
    double x_k = 1;
    double y0 = 0;
    double y_k = 1;
    double search_n = 3;

	static HeightTestingParams load(const std::filesystem::path& filePath) {

#define TRY_LOAD(param) do \
		if (json.contains(#param)) { json.at(#param).get_to(testing_params.param); } \
		while(false)


		nlohmann::json json;
		std::ifstream stream{ filePath };
		stream >> json;

		HeightTestingParams testing_params;
		json.at("mode").get_to(testing_params.mode);
		json.at("t").get_to(testing_params.t);

		TRY_LOAD(x0);
		TRY_LOAD(x_k);
		TRY_LOAD(y0);
		TRY_LOAD(y_k);
		TRY_LOAD(search_n);
		return testing_params;

#undef TRY_LOAD
	}

	void print() {
        std::cout << "Testing params:" << std::endl;
        std::cout << "mode: " << mode << std::endl;

        std::cout << "t: " << t << std::endl;
        std::cout << "x0: " << x0 << std::endl;
        std::cout << "x_k: " << x_k << std::endl;
        std::cout << "y0: " << y0 << std::endl;
        std::cout << "y_k: " << y_k << std::endl;        

        std::cout << "search_n: " << search_n << std::endl;
	}
};

