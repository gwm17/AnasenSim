#include "ReactionSystem.h"
#include "DecaySystem.h"
#include "OneStepSystem.h"
#include "TwoStepSystem.h"

namespace AnasenSim {

	ReactionSystem* CreateSystem(const SystemParameters& params)
	{
		switch(params.stepParams.size())
		{
			case 1:
			{
				if(params.stepParams[0].rxnType == RxnType::Decay)
					return new DecaySystem(params);
				else if (params.stepParams[0].rxnType == RxnType::Reaction)
					return new OneStepSystem(params);
			}
			case 2: return new TwoStepSystem(params);
		}

		return nullptr;
	}

	ReactionSystem::ReactionSystem(const SystemParameters& params) :
		m_params(params), m_isValid(true), m_sysEquation("")
	{
	}
	
	ReactionSystem::~ReactionSystem()
	{
	}
}