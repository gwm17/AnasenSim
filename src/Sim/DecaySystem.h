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
		void Init(const std::vector<StepParameters>& params);
		void SetSystemEquation() override;
	
		Reaction m_step1;
	};

}

#endif