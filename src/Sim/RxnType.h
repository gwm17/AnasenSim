#ifndef RXNTYPE_H
#define RXNTYPE_H

#include <string>

namespace AnasenSim {

	enum class RxnType
	{
		Decay,
		Reaction,
		None
	};

	enum RxnSize
	{
		DecaySize = 2,
		OneStepSize = 3,
		TwoStepSize = 4
	};

	static RxnType StringToRxnType(const std::string& type)
	{
		if (type == "Decay")
			return RxnType::Decay;
		else if (type == "Reaction")
			return RxnType::Reaction;
		else
			return RxnType::None;
	}
}

#endif