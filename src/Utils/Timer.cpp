#include "Timer.h"

namespace AnasenSim {

	Timer::Timer()
	{
		m_startTime = Clock::now();
		m_stopTime = m_startTime;
	}

	Timer::~Timer() {}

	void Timer::Start()
	{
		m_startTime = Clock::now();
	}

	void Timer::Stop()
	{
		m_stopTime = Clock::now();
	}

	double Timer::GetElapsedSeconds()
	{
		return std::chrono::duration_cast<std::chrono::duration<double>>(m_stopTime-m_startTime).count();
	}

	double Timer::GetElapsedMilliseconds()
	{
		return std::chrono::duration_cast<std::chrono::duration<double>>(m_stopTime-m_startTime).count()*1000.0;
	}

}