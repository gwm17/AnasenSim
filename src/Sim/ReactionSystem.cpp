#include "ReactionSystem.h"
#include "DecaySystem.h"
#include "OneStepSystem.h"
#include "TwoStepSystem.h"

namespace AnasenSim {

	ReactionSystem::ReactionSystem(const Target& target) :
		m_isValid(true), m_sysEquation(""), m_phiDist(0.0, M_PI*2.0), m_cosThetaDist(-1.0, 1.0),
		m_target(target)
	{
	}
	
	ReactionSystem::~ReactionSystem()
	{
	}

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

	//Each reaction step can generate an excited nucleus (for a reaction can make an excited residual, decay can make an excited
	//breakup2 or "heavy")
	void ReactionSystem::AddExcitationDistribution(double mean, double sigma)
	{
		if(mean == -1.0 || sigma == -1.0)
		{
			m_isValid = false;
			std::cerr << "Invalid parameters at ReactionSystem::AddExcitationDistribution(), distribution is invalid -> mean: " << mean
					  << " sigma: " << sigma << std::endl;
			return;
		}
		m_exDistributions.emplace_back(mean, sigma);
	}
}