#include "RandomGenerator.h"

namespace AnasenSim {
	
	RandomGenerator::RandomGenerator()
	{
		std::random_device rd;
		rng.seed(rd());
	}

	RandomGenerator::~RandomGenerator() {}
}