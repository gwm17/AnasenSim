#ifndef ONESTEPSYSTEM_H
#define ONESTEPSYSTEM_H

#include "ReactionSystem.h"

namespace AnasenSim {

	class OneStepSystem: public ReactionSystem
	{
	public:
		OneStepSystem(const SystemParameters& params);
		~OneStepSystem();
	
		void RunSystem() override;
	
	private:
		void Init(const std::vector<StepParameters>& params);
		virtual void SetSystemEquation() override;

		double m_rxnPathLength;
		double m_beamStraggling;
			
		Reaction m_step1;
	};

}

#endif