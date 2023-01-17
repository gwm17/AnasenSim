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
		void Init();
		virtual void SetSystemEquation() override;
		void SampleParameters();

		double m_rxnPathLength;
		double m_beamStraggling;
		double m_rxnBeamEnergy;
		double m_rxnTheta;
		double m_rxnPhi;
		double m_residEx;
		double m_beamTheta;
		double m_beamPhi;
			
		Reaction m_step1;
	};

}

#endif