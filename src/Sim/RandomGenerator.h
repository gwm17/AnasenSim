#ifndef RANDOMGENERATOR_H
#define RANDOMGENERATOR_H

#include <random>
#include <iostream>

namespace AnasenSim {

	class RandomGenerator
	{
	public:
		//Valid for float, double, long double
		template<typename T>
		static T GetUniformReal(T min, T max)
		{
			std::uniform_real_distribution<T> distribution(min, max);
			return distribution(GetGenerator());
		}

		//Valid for any integer type (signed or unsigned)
		template<typename T>
		static T GetUniformInt(T min, T max)
		{
			std::uniform_int_distribution<T> distribution(min, max);
			return distribution(GetGenerator());
		}

		//Valid for float, double, long double
		template<typename T>
		static T GetNormal(T mean, T sigma)
		{
			std::normal_distribution<T> distribution(mean, sigma);
			return distribution(GetGenerator());
		}

		//This is the most common use case, so we eliminate recreation of distribution, templating. 
		//For randomization of decimals in conversion from integer -> floating point for histograming
		static double GetUniformFraction()
		{
			static std::uniform_real_distribution<double> distribution(0.0, 1.0);
			return distribution(GetGenerator());
		}

	private:
		static std::mt19937_64& GetGenerator()
		{
			static thread_local auto seed = std::random_device()();
			static bool print = false; //For debugging, pickout a failing seed
			if(print)
			{
				std::cout << "Using seed: " << seed << std::endl;
				print = false;
			}
			static thread_local std::mt19937_64 generator(seed);
			return generator;
		}
	};

}

#endif