#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H

#include <random>

namespace AnasenSim {

	class RandomGenerator
	{
	public:
		~RandomGenerator();
		std::mt19937& GetGenerator() { return rng; }
		static RandomGenerator& GetInstance()
		{
			static thread_local RandomGenerator s_instance;
			return s_instance;
		}

	private:
		RandomGenerator();
		
		std::mt19937 rng;
	};

}

#endif