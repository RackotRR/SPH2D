#pragma once
#include <utility>
#include <array>
 
template <typename T, size_t dimensions, size_t size>
class heap_array_md {
private:
	T* ptr;
	void create(T initValue) {
		ptr = new T[dimensions * size]{ initValue }; 
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
		size_t count{ dimensions * size };
		for (size_t i{}; i < count; i++) {
			arr.ptr[i] = ptr[i];
		}
		return arr;
	}

	const T& operator() (size_t i, size_t j) const {
		return ptr[i + dimensions * j];
	}
	T& operator() (size_t i, size_t j) {
		return ptr[i + dimensions * j];
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
		size_t count{ dimensions * size };
		return ptr + count;
	}

	void fill(const T& value) {
		size_t count{ dimensions * size };
		for (size_t i = 0; i < count; ++i) {
			ptr[i] = value;
		}
	}
};

template<typename T, size_t size>
class heap_array {
private:
	T* ptr;
	void create(T initValue) {
		ptr = new T[size];
		for (size_t i{}; i < size; i++) {
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
		for (size_t i{}; i < size; i++) {
			arr.ptr[i] = ptr[i];
		}
		return arr;
	}

	heap_array() : heap_array(T{}) {

	}

	explicit heap_array(T initValue) {
		create(initValue);
	}

	const T& operator() (size_t i) const {
		return ptr[i];
	} 
	T& operator() (size_t i) {
		return ptr[i];
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
		return ptr + size;
	}

	void fill(const T& value) {
		for (size_t i = 0; i < size; ++i) {
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