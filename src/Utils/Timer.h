#ifndef TIMER_H
#define TIMER_H

#include <chrono>

namespace AnasenSim {

    class Timer
    {
    public:
		Timer();
		~Timer();
		void Start();
		void Stop();
		double GetElapsedSeconds();
		double GetElapsedMilliseconds();

	private:
		using Time = std::chrono::high_resolution_clock::time_point;
		using Clock = std::chrono::high_resolution_clock;

		Time m_startTime;
        Time m_stopTime;
    };

}

#endif