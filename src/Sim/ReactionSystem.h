#ifndef REACTIONSYSTEM_H
#define REACTIONSYSTEM_H

#include "SimBase.h"
#include "RxnType.h"
#include "Reaction.h"
#include "Target.h"
#include <vector>
#include <random>

namespace AnasenSim {

	struct StepParameters
	{
		RxnType rxnType = RxnType::None;
		std::vector<int> Z;
		std::vector<int> A;
		double meanResidualEx = -1.0;
		double sigmaResidualEx = -1.0;
	};

	struct SystemParameters
	{
		Target target;
		double initialBeamEnergy = 0.0;
		double rxnBeamEnergy = 0.0;
		std::vector<StepParameters> stepParams;
	};

	class ReactionSystem
	{
	public:
		ReactionSystem(const SystemParameters& params);
		virtual ~ReactionSystem();

		virtual void RunSystem() = 0;

		std::vector<Nucleus>* GetNuclei() { return &m_nuclei; }
		const std::string& GetSystemEquation() const { return m_sysEquation; }
		bool IsValid() const { return m_isValid; }
		//Need to reset the detected status of the nulcei after they're written to disk
		void ResetNucleiDetected();

	protected:
		virtual void SetSystemEquation() = 0;

		SystemParameters m_params;

		bool m_isValid;

		std::string m_sysEquation;
		std::vector<Nucleus> m_nuclei;

		static constexpr double s_deg2rad = M_PI/180.0;
		static constexpr double s_cosThetaMin = -1.0;
		static constexpr double s_cosThetaMax = 1.0;
		static constexpr double s_phiMin = 0.0;
		static constexpr double s_phiMax = 2.0*M_PI;
	};

	ReactionSystem* CreateSystem(const SystemParameters& params);
}

#endif