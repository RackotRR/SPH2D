#pragma once
#include <string>
#include "Params.h"

namespace ParamsGeneration {
	void makeHeader(const std::string& path, const ExperimentParams& params);
	void makeJson(const std::string& path, const ExperimentParams& params);

	void makeParamsGeneratorClass(const std::string& path);
}