#ifndef DECAYSYSTEM_H
#define DECAYSYSTEM_H

#include "ReactionSystem.h"

namespace AnasenSim {

	class DecaySystem: public ReactionSystem
	{
	public:
		DecaySystem(const SystemParameters& params);
		~DecaySystem();
	
		virtual void RunSystem() override;
	
	private:
		void Init();
		void SetSystemEquation() override;
		void SampleParameters();
	
		Reaction m_step1;
		double m_rxnTheta;
		double m_rxnPhi;
		double m_ex;
	};

}

#endif