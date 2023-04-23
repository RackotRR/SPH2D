#pragma once
#include <string>
#include <Params.h>

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
		SPHFIO(const std::string& experiment_name);

		ParamsPtr getParams() const;

		Grid makeGrid() const;
		LazyGrid makeLazyGrid() const;

		bool isAdditionalValuePresented(const std::string& value) const;
	public:
		const Directories directories;
	private:
		ParamsPtr loadExperimentParams();
		static LayersPathPtr findTimeLayersPath(ParamsPtr params);
	private:
		ParamsPtr params;
		LayersPathPtr available_layers_path;
	};

}
