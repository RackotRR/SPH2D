#pragma once
#include <utility>
#include <array>

namespace RR::Memory {
	template<typename T>
	class heap_darray_md {
	private:
		size_t dimensions;
		size_t size_in_dim;
		T* ptr;

		T* create(size_t dimensions, size_t size_in_dim, T initValue) {
			if (dimensions == 0 || size_in_dim == 0) return nullptr;
			size_t elements = dimensions * size_in_dim;
			T* ptr = new T[elements];
			for (size_t i{}; i < elements; i++) {
				ptr[i] = initValue;
			}
			return ptr;
		}

	public:
		~heap_darray_md() {
			delete[] ptr;
		}


		heap_darray_md(heap_darray_md&& arr) noexcept {
			dimensions = std::exchange(arr.dimensions, 0);
			size_in_dim = std::exchange(arr.size_in_dim, 0);
			ptr = std::exchange(arr.ptr, nullptr);
		}
		heap_darray_md& operator=(heap_darray_md&& arr) noexcept {
			std::swap(dimensions, arr.dimensions);
			std::swap(size_in_dim, arr.size_in_dim);
			std::swap(ptr, arr.ptr);
			return *this;
		}

		heap_darray_md(heap_darray_md&) = delete;
		auto operator=(heap_darray_md&) = delete;

		heap_darray_md() : heap_darray_md(0, 0, T{})
		{
		}

		heap_darray_md(size_t dimensions, size_t size_in_dim) : heap_darray_md(dimensions, size_in_dim, T{})
		{
		}

		heap_darray_md(size_t dimensions, size_t size_in_dim, T init_value) :
			dimensions{ dimensions },
			size_in_dim{ size_in_dim },
			ptr{ create(dimensions, size_in_dim, init_value) }
		{
		}

		auto copy() const {
			heap_darray_md arr{ dimensions, size_in_dim };
			for (size_t i{}; i < size(); i++) {
				arr.ptr[i] = ptr[i];
			}
			return arr;
		}


		const T& at(size_t i, size_t j) const {
			assert(i < dimensions && j < size_in_dim);
			return (*this)(i, j);
		}
		T& at(size_t i, size_t j) {
			assert(i < dimensions&& j < size_in_dim);
			return (*this)(i, j);
		}
		const T& operator() (size_t i, size_t j) const {
			assert(i < dimensions&& j < size_in_dim);
			return ptr[i + dimensions * j];
		}
		T& operator() (size_t i, size_t j) {
			assert(i < dimensions&& j < size_in_dim);
			return ptr[i + dimensions * j];
		}
		const T& operator[] (size_t k) const {
			assert(k < size());
			return ptr[k];
		}
		T& operator[] (size_t k) {
			assert(k < size());
			return ptr[k];
		}

		T* data() {
			return ptr;
		}
		const T* data() const {
			return ptr;
		}
		T* begin() {
			return ptr;
		}
		T* end() {
			return ptr + size();
		}
		size_t size() const {
			return dimensions * size_in_dim;
		}

		void fill(const T& value) {
			for (size_t i = 0; i < size(); ++i) {
				ptr[i] = value;
			}
		}
	};


	template<typename T>
	class heap_darray {
	private:
		size_t elements;
		T* ptr;

		static T* create(size_t elements, T initValue) {
			if (elements == 0) return nullptr;
			T* ptr = new T[elements];
			for (size_t i{}; i < elements; i++) {
				ptr[i] = initValue;
			}
			return ptr;
		}
	public:
		~heap_darray() {
			delete[] ptr;
		}

		heap_darray(heap_darray&& arr) noexcept {
			elements = std::exchange(arr.elements, 0);
			ptr = std::exchange(arr.ptr, nullptr);
		}
		heap_darray& operator=(heap_darray&& arr) noexcept {
			std::swap(elements, arr.elements);
			std::swap(ptr, arr.ptr);
			return *this;
		}

		heap_darray(heap_darray&) = delete;
		auto operator=(heap_darray&) = delete;

		auto copy() const {
			heap_darray arr{ elements };
			for (size_t i{}; i < elements; i++) {
				arr.ptr[i] = ptr[i];
			}
			return arr;
		}

		heap_darray() : heap_darray{ 0, T{} } 
		{
		}

		explicit heap_darray(size_t elements) : heap_darray(elements, T{}) 
		{
		}

		heap_darray(size_t elements, T initValue) : 
			elements{ elements },
			ptr{ create(elements, initValue) }
		{
		}

		const T& at(size_t i) const {
			assert(i < size());
			return (*this)(i);
		}
		T& at(size_t i) {
			assert(i < size());
			return (*this)(i);
		}
		const T& operator() (size_t i) const {
			assert(i < size());
			return ptr[i];
		}
		T& operator() (size_t i) {
			assert(i < size());
			return ptr[i];
		}
		const T& operator[] (size_t k) const {
			assert(k < size());
			return ptr[k];
		}
		T& operator[] (size_t k) {
			assert(k < size());
			return ptr[k];
		}

		T* data() {
			return ptr;
		}
		const T* data() const {
			return ptr;
		}
		T* begin() {
			return ptr;
		}
		T* end() {
			return ptr + elements;
		}
		size_t size() const {
			return elements;
		}

		void fill(const T& value) {
			for (size_t i = 0; i < elements; ++i) {
				ptr[i] = value;
			}
		}
	};
}