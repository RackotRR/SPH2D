#pragma once
#include <vector>
#include <memory>
#include <filesystem>

#include <Params.h>
#include <RR/Memory/HeapDArray.h>

#include "ExperimentLayer.h"
#include "ExperimentLayers.h"

namespace sphfio {
	using ParamsPtr = std::shared_ptr<ExperimentParams>;

	struct TimeLayer {
		TimeLayer() = default;
		TimeLayer(const ExperimentLayer& experiment_layer, ParamsPtr params);

		vheap_darray_floatn r_var;
		RR::Memory::heap_darray<rr_int> itype;
		vheap_darray_floatn v_var;
		RR::Memory::heap_darray<rr_float> p;
		RR::Memory::heap_darray<rr_float> rho;

		rr_uint ntotal = 0;
		rr_float time = 0;

		rr_float getByTag(const std::string& value, rr_uint i) const;
		ParamsPtr params;
	};

	/**
	 * @brief grid with layers loading on demand
	*/
	class LazyGrid {
	public:
		using time_points_t = std::vector<rr_float>;

		/**
		 * @brief wrapper aroung ExperimentLayers::iterator for TimeLayer construction on access
		*/
		class Iterator {
		public:
			Iterator(const LazyGrid& lazy_grid, ExperimentLayers::iterator iter);

			[[nodiscard]] TimeLayer operator*() const;
			Iterator& operator++();
			Iterator operator++(int);
			bool operator==(const Iterator& other) const;
			bool operator!=(const Iterator& other) const;
		private:
			const LazyGrid& lazy_grid;
			ExperimentLayers::iterator iter;
		};
	public:
		LazyGrid(const ExperimentLayers::Ptr experiment_layers, ParamsPtr params);

		Iterator begin() const;
		Iterator end() const;
		Iterator find(rr_float time) const;

		size_t size() const;
		bool empty() const;

		[[nodiscard]] time_points_t time_points() const;
	private:
		const ExperimentLayers::Ptr experiment_layers;
		ParamsPtr params;
	};

	/**
	 * @brief grid with layers loading on construction
	*/
	class Grid {
	public:
		using iterator = std::vector<TimeLayer>::const_iterator;
		using time_layers_t = std::vector<TimeLayer>;
		using time_points_t = std::vector<rr_float>;
	public:
		Grid(const ExperimentLayers::Ptr experiment_layers, ParamsPtr params);

		iterator begin() const;
		iterator end() const;
		iterator find(rr_float time) const;

		const TimeLayer& at(size_t layer_num) const;
		size_t size() const;
		bool empty() const;

		[[nodiscard]] time_points_t time_points() const;
	private:
		time_layers_t grid;
		const ExperimentLayers::Ptr experiment_layers;
	};
}