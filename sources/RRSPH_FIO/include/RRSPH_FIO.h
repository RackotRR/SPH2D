#pragma once
#include <string>
#include <filesystem>
#include <unordered_set>
#include "Params.h"

#include "ExperimentDirectory.h"
#include "ExperimentDirectories.h"

#include "Grid.h"
#include "Directories.h"

namespace sphfio {

	struct Square {
		Square() = default;
		Square(const ParamsPtr& params) : 
			origin{ 
				params->x_mingeom,
				params->y_mingeom,
				params->z_mingeom, // equals 0 in 2D
			},
			size{
				params->x_maxgeom - params->x_mingeom,
				params->y_maxgeom - params->y_mingeom,
				params->z_maxgeom - params->z_mingeom, // equals 0 in 2D
			}
		{
		}

		inline bool contains2D(const rr_float2& r) const {
			return origin.x <= r.x && origin.x + size.x >= r.x &&
				origin.y <= r.y && origin.y + size.y >= r.y;
		}
		inline bool contains3D(const rr_float3& r) const {
			return contains2D({r.x, r.y}) && 
				origin.z <= r.z && origin.z + size.z >= r.z;
		}

		template<typename rr_floatn>
		bool contains(rr_floatn r) const {
			if constexpr (is_using_float3<rr_floatn>()) {
				return contains3D(r);
			}
			else {
				return contains2D(r);
			}
		}

		rr_float3 origin{};
		rr_float3 size{};
	};

	/**
	 * @brief user interface for ExperimentDirectory selection with specific properties
	*/
	ExperimentDirectory::Ptr CLI(
		const ExperimentDirectory::properties_t& properties,
		std::filesystem::path experiments_directory = std::filesystem::current_path());

	class SPHFIO {
	public:
		/**
		 * @brief user interface for ExperimentDirectory selection with default settings
		*/
		SPHFIO();

		/**
		 * @brief loading ExperimentDirectory from path
		*/
		SPHFIO(const std::filesystem::path& experiment_dir);

		/**
		 * @brief common ctor
		*/
		SPHFIO(ExperimentDirectory::Ptr experiment_dir);

		ParamsPtr getParams() const;

		Grid makeGrid() const;
		LazyGrid makeLazyGrid() const;

		bool isAdditionalValuePresented(const std::string& value) const;
	public:
		const Directories directories;
	private:
		ParamsPtr loadExperimentParams();
		ExperimentLayers::Ptr find_time_layers_path() const;
		std::unordered_set<std::string> findAvailableVariables(ParamsPtr params);
	private:
		ExperimentDirectory::Ptr experiment;
		ParamsPtr params;
		std::unordered_set<std::string> available_variables;
	};

}
