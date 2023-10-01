#include "HeightTestingParams.h"

void SpaceTestingParams::print() {
    std::cout << "SpaceTesting params:" << std::endl;

    std::cout << "t: " << t << std::endl;
    std::cout << "x0: " << x0 << std::endl;
    std::cout << "x_k: " << x_k << std::endl;
    std::cout << "y0: " << y0 << std::endl;
    std::cout << "y_k: " << y_k << std::endl;

    std::cout << "search_n: " << search_n << std::endl;
}

void TimeTestingParams::print() {
	std::cout << "TimeTesting params:" << std::endl;

	std::cout << "x: " << x << std::endl;
	std::cout << "t0: " << t0 << std::endl;
	std::cout << "t_k: " << t_k << std::endl;
	std::cout << "y0: " << y0 << std::endl;
	std::cout << "y_k: " << y_k << std::endl;

	std::cout << "search_n: " << search_n << std::endl;
}


std::shared_ptr<HeightTestingParams> HeightTestingParams::load(const std::filesystem::path& experiment_dir) {
	nlohmann::json json;
	std::ifstream stream{ experiment_dir / HeightTestingParams::filename };
	stream >> json;

	std::string mode;
	json.at("mode").get_to(mode);

	std::shared_ptr<HeightTestingParams> testing_params;
	if (mode == "space") {
		auto space_testing_params = std::make_shared<SpaceTestingParams>();
		json.at("t").get_to(space_testing_params->t);

		if (json.contains("x0")) json.at("x0").get_to(space_testing_params->x0);
		if (json.contains("x_k")) json.at("x_k").get_to(space_testing_params->x_k);

		testing_params = space_testing_params;
	}
	else if (mode == "time") {
		auto time_testing_params = std::make_shared<TimeTestingParams>();
		json.at("x").get_to(time_testing_params->x);

		if (json.contains("t0")) json.at("t0").get_to(time_testing_params->t0);
		if (json.contains("t_k")) json.at("t_k").get_to(time_testing_params->t_k);

		testing_params = time_testing_params;
	}
	else {
		throw std::runtime_error{ "Wrong mode passed: " + mode };
	}

	testing_params->mode = mode;
	if (json.contains("y0")) json.at("y0").get_to(testing_params->y0);
	if (json.contains("y_k")) json.at("y_k").get_to(testing_params->y_k);
	if (json.contains("search_n")) json.at("search_n").get_to(testing_params->search_n);
	return testing_params;
}

void HeightTestingParams::generate_default(const std::filesystem::path& experiment_dir) {
	std::filesystem::path path = experiment_dir / HeightTestingParams::filename;	

	if (!std::filesystem::exists(path)) {
		std::ofstream stream{ path };
		nlohmann::json json;
		json["mode"] = "space";
		json["t"] = 0;
		json["x"] = nullptr;
		json["x0"] = 0;
		json["x_k"] = 1;
		json["y0"] = 1;
		json["y_k"] = 1;
		json["t0"] = nullptr;
		json["t_k"] = nullptr;
		stream << json.dump(4) << std::endl;
	}
}