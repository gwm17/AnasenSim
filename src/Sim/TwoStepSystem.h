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
		void Init();
		void SetSystemEquation() override;
		void SampleParameters();

		//reaction parameters
		double m_rxnPathLength;
		double m_beamStraggling;
		double m_rxnBeamEnergy;
		double m_rxnTheta;
		double m_rxnPhi;
		double m_decay1Theta;
		double m_decay1Phi;
		double m_residEx;
		double m_decay2Ex;
		double m_beamTheta;
		double m_beamPhi;

		Reaction m_step1, m_step2;
	};

}

#endif