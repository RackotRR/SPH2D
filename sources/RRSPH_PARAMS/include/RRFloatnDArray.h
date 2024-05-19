#pragma once
#include "RRFloat2.h"
#include "RRFloat3.h"
#include "RR/Memory/HeapDArray.h"
#include <variant>

class vheap_darray_floatn {
	using flt2 = RR::Memory::heap_darray<rr_float2>;
	using flt3 = RR::Memory::heap_darray<rr_float3>;

public:
	std::variant<flt2, flt3> data;
	static inline rr_uint dimensions = 0;
	static void set_dimenstions(rr_uint dimensions) {
		vheap_darray_floatn::dimensions = dimensions;
	}

	vheap_darray_floatn() = default;
	vheap_darray_floatn(rr_uint n) {
		if (dimensions == 2) {
			data.emplace<flt2>(n);
		}
		else if (dimensions == 3) {
			data.emplace<flt3>(n);
		}
		else {
			throw std::exception{ "vheap_darray_floatn ctor error: invalid dimensions value" };
		}
	}

	flt2& get_flt2() {
		return std::get<flt2>(data);
	}
	const flt2& get_flt2() const {
		return std::get<flt2>(data);
	}

	flt3& get_flt3() {
		return std::get<flt3>(data);
	}
	const flt3& get_flt3() const {
		return std::get<flt3>(data);
	}
};


class vheap_darray_floatn_md {
	using flt2 = RR::Memory::heap_darray_md<rr_float2>;
	using flt3 = RR::Memory::heap_darray_md<rr_float3>;

	static inline rr_uint dimensions = 0;
	std::variant<flt2, flt3> data;
public:
	static void set_dimenstions(rr_uint dimensions) {
		vheap_darray_floatn_md::dimensions = dimensions;
	}

	vheap_darray_floatn_md(rr_uint n, rr_uint m) {
		if (dimensions == 2) {
			data = flt2(n, m);
		}
		else if (dimensions == 3) {
			data = flt3(n, m);
		}
		else {
			assert(0);
		}
	}

	flt2& get_flt2() {
		return std::get<flt2>(data);
	}
	const flt2& get_flt2() const {
		return std::get<flt2>(data);
	}

	flt3& get_flt3() {
		return std::get<flt3>(data);
	}
	const flt3& get_flt3() const {
		return std::get<flt3>(data);
	}
};

using shared_vheap_darray_floatn = std::shared_ptr<vheap_darray_floatn>;
