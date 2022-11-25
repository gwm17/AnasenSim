#ifndef TWOSTEPSYSTEM_H
#define TWOSTEPSYSTEM_H

#include "ReactionSystem.h"

namespace AnasenSim {

	class TwoStepSystem : public ReactionSystem
	{
	public:
		TwoStepSystem(const SystemParameters& params);
		~TwoStepSystem();

		virtual void RunSystem() override;
	
	private:
		void Init(const std::vector<StepParameters>& params);
		void SetSystemEquation() override;

		double m_rxnPathLength;
		double m_beamStraggling;

		Reaction m_step1, m_step2;
	};

}

#endif