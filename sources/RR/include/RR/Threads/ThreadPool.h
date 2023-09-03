#pragma once
#include <thread>
#include <queue>

namespace RR::Threads {

    class ThreadPool {
    public:
        ~ThreadPool() {
            while (!threads.empty()) {
                threads.front().join();
                threads.pop();
            }
        }

        void add_thread(std::thread&& thread) {
            refresh();
            threads.push(std::move(thread));
        }
        void refresh() {
            while (!threads.empty() && 
                threads.front().get_id() == std::thread::id{})
            {
                threads.front().join();
                threads.pop();
            }
        }
    private:
        std::queue<std::thread> threads;
    };


}