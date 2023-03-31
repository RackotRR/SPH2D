#pragma once
#include <utility>
#include <array>

namespace RR::Memory {
	template <typename T, size_t dimensions, size_t size_in_dim>
	class heap_array_md {
	private:
		T* ptr;
		void create(T initValue) {
			ptr = new T[dimensions * size_in_dim];
			for (size_t i{}; i < dimensions * size_in_dim; i++) {
				ptr[i] = initValue;
			}
		}
	public:
		~heap_array_md() { 
			delete[] ptr;
		}


		heap_array_md(heap_array_md&& arr) noexcept {
			ptr = std::exchange(arr.ptr, nullptr);
		}
		heap_array_md& operator=(heap_array_md&& arr) noexcept {
			std::swap(ptr, arr.ptr);
			return *this;
		}

		heap_array_md(heap_array_md&) = delete;
		auto operator=(heap_array_md&) = delete;

		heap_array_md() : heap_array_md(T{}) {

		}

		explicit heap_array_md(T initValue) {
			create(initValue);
		}

		auto copy() const {
			heap_array_md arr;
			for (size_t i{}; i < size(); i++) {
				arr.ptr[i] = ptr[i];
			}
			return arr;
		}


		const T& at(size_t i, size_t j) const {
			return (*this)(i, j);
		}
		T& at(size_t i, size_t j) {
			return (*this)(i, j);
		}
		const T& operator() (size_t i, size_t j) const {
			return ptr[i + dimensions * j];
		}
		T& operator() (size_t i, size_t j) {
			return ptr[i + dimensions * j];
		}
		const T& operator[] (size_t k) const {
			return ptr[k];
		}
		T& operator[] (size_t k) {
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
		constexpr size_t size() const {
			return dimensions * size_in_dim;
		}

		void fill(const T& value) {
			for (size_t i = 0; i < size(); ++i) {
				ptr[i] = value;
			}
		}
	};

	template<typename T, size_t elements>
	class heap_array {
	private:
		T* ptr;
		void create(T initValue) {
			ptr = new T[elements];
			for (size_t i{}; i < elements; i++) {
				ptr[i] = initValue;
			}
		}
	public:
		~heap_array() {
			delete[] ptr;
		}

		heap_array(heap_array&& arr) noexcept {
			ptr = std::exchange(arr.ptr, nullptr);
		}
		heap_array& operator=(heap_array&& arr) noexcept {
			std::swap(ptr, arr.ptr);
			return *this;
		}

		heap_array(heap_array&) = delete;
		auto operator=(heap_array&) = delete;

		auto copy() const {
			heap_array arr;
			for (size_t i{}; i < elements; i++) {
				arr.ptr[i] = ptr[i];
			}
			return arr;
		}

		heap_array() : heap_array(T{}) {

		}

		explicit heap_array(T initValue) {
			create(initValue);
		}

		const T& at(size_t i) const {
			return (*this)(i);
		}
		T& at(size_t i) {
			return (*this)(i);
		}
		const T& operator() (size_t i) const {
			return ptr[i];
		} 
		T& operator() (size_t i) {
			return ptr[i];
		}
		const T& operator[] (size_t k) const {
			return ptr[k];
		}
		T& operator[] (size_t k) {
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
		constexpr size_t size() const {
			return elements;
		}

		void fill(const T& value) {
			for (size_t i = 0; i < elements; ++i) {
				ptr[i] = value;
			}
		}
	};

	template<typename T, size_t size>
	class stack_array : public std::array<T, size> {
	public:
		T& operator() (size_t i) {
			return (*this)[i];
		}
		const T& operator() (size_t i) const {
			return (*this)[i];
		}
	};

}