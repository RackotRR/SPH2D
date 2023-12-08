#pragma once
#include <string>
#include <filesystem>
#include <unordered_set>
#include "Params.h"

#include "Grid.h"
#include "Directories.h"

namespace sphfio {

	struct Square {
		Square() = default;
		Square(const ParamsPtr& params);

		bool contains(rr_float x, rr_float y) const;
		bool contains(rr_float2 r) const;

		rr_float origin_x;
		rr_float origin_y;
		rr_float size_x;
		rr_float size_y;
	};

	class SPHFIO {
	public:
		SPHFIO();
		SPHFIO(const std::filesystem::path& experiment_dir);

		ParamsPtr getParams() const;

		Grid makeGrid() const;
		LazyGrid makeLazyGrid() const;

		bool isAdditionalValuePresented(const std::string& value) const;
	public:
		const Directories directories;
	private:
		ParamsPtr loadExperimentParams();
		LayersPathPtr findTimeLayersPath();
		std::unordered_set<std::string> findAvailableVariables(ParamsPtr params);
	private:
		ParamsPtr params;
		LayersPathPtr available_layers_path;
		std::unordered_set<std::string> available_variables;
	};

}
