#include "HeightTestingParams.h"

static void print_vec(const std::vector<double>& vec, std::ostream& stream) {
	int i = 0;
	for (double val : vec) {
		stream << val;
		++i;
		if (i < vec.size()) {
			stream << ", ";
		}
	}
}

void HeightTestingParams::print_common(std::ostream& stream) {
	stream << "y0: " << y0 << std::endl;
	stream << "y_k: " << y_k << std::endl;

	stream << "search_n: " << search_n << std::endl;

	if (particles_type.has_value()) {
		stream << "particles_type: " << particles_type.value() << std::endl;
	}

	if (!postfix.empty()) {
		stream << "postfix: " << postfix << std::endl;
	}
}

void SpaceTestingParams::print(std::ostream& stream) {
    stream << "SpaceTesting params:" << std::endl;

	stream << "t: ";
	print_vec(t, stream);
	stream << std::endl;

    stream << "x0: " << x0 << std::endl;
    stream << "x_k: " << x_k << std::endl;

	print_common(stream);
}

void TimeTestingParams::print(std::ostream& stream) {
	stream << "TimeTesting params:" << std::endl;

	stream << "x: ";
	print_vec(x, stream);
	stream << std::endl;

	stream << "t0: " << t0 << std::endl;
	stream << "t_k: " << t_k << std::endl;

	print_common(stream);
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
	if (json.contains("postfix")) json.at("postfix").get_to(testing_params->postfix);
	if (json.contains("particles_type") && json.at("particles_type").is_number()) testing_params->particles_type = json.at("particles_type").get<int>();
	return testing_params;
}

void HeightTestingParams::generate_default(const std::filesystem::path& experiment_dir) {
	std::filesystem::path path = experiment_dir / HeightTestingParams::filename;	

	if (!std::filesystem::exists(path)) {
		std::ofstream stream{ path };
		nlohmann::json json;

		SpaceTestingParams default_space_testing;

		json["mode"] = "space";
		json["t"] = std::vector{ 0.0 };
		json["x"] = nullptr;
		json["x0"] = default_space_testing.x0;
		json["x_k"] = default_space_testing.x_k;
		json["y0"] = default_space_testing.y0;
		json["y_k"] = default_space_testing.y_k;
		json["t0"] = nullptr;
		json["t_k"] = nullptr;
		json["search_n"] = default_space_testing.search_n;
		json["particles_type"] = nullptr;
		json["postfix"] = "";
		stream << json.dump(4) << std::endl;
	}
}