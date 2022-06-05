#pragma once
#include <thread>
#include <functional>

void DoMTh(int begin, int end, int nThreads, std::function<void(int begin, int end)> functor) {
	std::vector<std::jthread> threads;

	int byThread = (end - begin) / nThreads;
	for (int i = 0; i < nThreads; i++) {
		threads.emplace_back(functor, begin, begin + byThread);
		begin += byThread;
	}

	if ((end - begin) % nThreads != 0) {
		threads.emplace_back(functor, begin, end);
	}
}