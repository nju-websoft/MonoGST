#ifndef Timer_H
#define Timer_H

#include <chrono>
#include <stdexcept>
class Timer 
{
    private:
        std::chrono::high_resolution_clock::time_point start_time;
        bool running = false;

    public:
        // 开始计时
        void start() {
            start_time = std::chrono::high_resolution_clock::now();
            running = true;
        }

        // 结束计时并返回时间 (ms)
        double stop() {
            if (!running) {
                throw std::runtime_error("Timer has not been started.");
            }
            auto end_time = std::chrono::high_resolution_clock::now();
            running = false;
            std::chrono::duration<double, std::milli> elapsed = end_time - start_time;
            return elapsed.count();
        }
};

#endif