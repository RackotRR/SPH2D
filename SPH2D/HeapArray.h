#pragma once
 
template <typename T, size_t dimensions, size_t size>
class heap_array_md {
private:
	T* ptr;
	void Create(T initValue) {
		ptr = new T[dimensions * size]{ initValue }; 
	}
public:
	~heap_array_md() { 
		delete[] ptr;
	}


	heap_array_md<T, dimensions, size>(heap_array_md<T, dimensions, size>&& arr) noexcept {
		ptr = std::exchange(arr.ptr, nullptr);
	}
	heap_array_md<T, dimensions, size>& operator=(heap_array_md<T, dimensions, size>&& arr) noexcept {
		std::swap(ptr, arr.ptr);
		return *this;
	}

	heap_array_md<T, dimensions, size>(heap_array_md<T, dimensions, size>&) = delete;
	auto operator=(heap_array_md&) = delete;

	heap_array_md() : heap_array_md(T{}) {

	}

	explicit heap_array_md(T initValue) {
		Create(initValue);
	}

	auto MakeCopy() const {
		heap_array_md<T, dimensions, size> arr;
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
};

template<typename T, size_t size>
class heap_array {
private:
	T* ptr;
	void Create(T initValue) {
		ptr = new T[size];
		for (size_t i{}; i < size; i++) {
			ptr[i] = initValue;
		}
	}
public:
	~heap_array() {
		delete[] ptr;
	}

	heap_array<T, size>(heap_array<T, size>&& arr) noexcept {
		ptr = std::exchange(arr.ptr, nullptr);
	}
	heap_array<T, size>& operator=(heap_array<T, size>&& arr) noexcept {
		std::swap(ptr, arr.ptr);
		return *this;
	}

	heap_array<T, size>(heap_array<T, size>&) = delete;
	auto operator=(heap_array<T, size>&) = delete;

	auto MakeCopy() const {
		heap_array<T, size> arr;
		for (size_t i{}; i < size; i++) {
			arr.ptr[i] = ptr[i];
		}
		return arr;
	}

	heap_array() : heap_array(T{}) {

	}

	explicit heap_array(T initValue) {
		Create(initValue);
	}

	const T& operator() (size_t i) const {
		return ptr[i];
	} 
	T& operator() (size_t i) {
		return ptr[i];
	} 
};