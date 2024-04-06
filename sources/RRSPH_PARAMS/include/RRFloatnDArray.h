#pragma once
#include "RRFloat2.h"
#include "RRFloat3.h"
#include "RR/Memory/HeapDArray.h"
#include<variant>

class vheap_darray_floatn {
	using flt2 = RR::Memory::heap_darray<rr_float2>;
	using flt3 = RR::Memory::heap_darray<rr_float3>;

	static inline rr_uint dimensions = 0;
	std::variant<flt2, flt3> data;
public:
	static void set_dimenstions(rr_uint dimensions) {
		vheap_darray_floatn::dimensions = dimensions;
	}

	vheap_darray_floatn(rr_uint n) {
		if (params.dim == 2) {
			data = flt2(n);
		}
		else if (params.dim == 3) {
			data = flt3(n);
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
		if (params.dim == 2) {
			data = flt2(n, m);
		}
		else if (params.dim == 3) {
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