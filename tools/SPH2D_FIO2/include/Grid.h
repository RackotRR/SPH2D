#pragma once
#include <vector>
#include <memory>

#include <Params.h>
#include <RR/Memory/HeapDArray.h>

namespace sphfio {
	using LayersPath = std::vector<std::string>;
	using LayersPathPtr = std::shared_ptr<LayersPath>;
	using ParamsPtr = std::shared_ptr<ExperimentParams>;

	struct TimeLayer {
		TimeLayer() = default;
		TimeLayer(const std::string& filename, rr_uint maxn);

		RR::Memory::heap_darray<rr_float2> r;
		RR::Memory::heap_darray<rr_int> itype;
		RR::Memory::heap_darray<rr_float2> v;
		RR::Memory::heap_darray<rr_float> p;

		rr_uint ntotal = 0;

		rr_float getByTag(const std::string& value, rr_uint i) const;
	};

	struct LazyGrid {
		LazyGrid(LayersPathPtr available_layers_path, ParamsPtr params);

		class Iterator {
		public:
			Iterator(const LazyGrid& lazy_grid, size_t current = 0);

			[[nodiscard]] TimeLayer operator*() const;
			Iterator& operator++();
			Iterator operator++(int);
			bool operator==(const Iterator& other) const;
			bool operator!=(const Iterator& other) const;
		private:
			const LazyGrid& lazy_grid;
			size_t current;
		};

		Iterator begin() const;
		Iterator end() const;

		Iterator find(size_t layer_num) const;

		size_t size() const;
		bool empty() const;
	private:
		LayersPathPtr available_layers_path;
		ParamsPtr params;
	};

	struct Grid {
		Grid(LayersPathPtr available_layers_path, ParamsPtr params);

		std::vector<TimeLayer>::const_iterator begin() const;
		std::vector<TimeLayer>::const_iterator end() const;
		std::vector<TimeLayer>::iterator begin();
		std::vector<TimeLayer>::iterator end();

		std::vector<TimeLayer>::const_iterator find(size_t layer_num) const;

		const TimeLayer& at(size_t layer_num) const;
		size_t size() const;
		bool empty() const;
	private:
		ParamsPtr params;
		std::vector<TimeLayer> grid;
	};
}