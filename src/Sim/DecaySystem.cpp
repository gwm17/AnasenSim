#include "DecaySystem.h"
#include "RandomGenerator.h"

#include <sstream>

namespace AnasenSim {

	DecaySystem::DecaySystem(const SystemParameters& params) :
		ReactionSystem(params)
	{
		m_nuclei.resize(3);
		Init(m_params.stepParams);
	}
	
	DecaySystem::~DecaySystem() {}
	
	void DecaySystem::Init(const std::vector<StepParameters>& params)
	{
		if(params.size() != 1 || params[0].rxnType != RxnType::Decay ||
		   params[0].Z.size() != 2 || params[0].A.size() != 2)
		{
			m_isValid = false;
			std::cerr << "Invalid parameters at DecaySystem::Init(), does not match Decay signature!" << std::endl;
			return;
		}

		const StepParameters& step1Params = params[0];

		int zr = step1Params.Z[0] - step1Params.Z[1];
		int ar = step1Params.A[0] - step1Params.A[1];

		m_nuclei[0] = CreateNucleus(step1Params.Z[0], step1Params.A[0]); //target
		m_nuclei[1] = CreateNucleus(step1Params.Z[1], step1Params.A[1]); //breakup1
		m_nuclei[2] = CreateNucleus(zr, ar); //breakup2

		m_step1.BindNuclei(&(m_nuclei[0]), nullptr, &(m_nuclei[1]), &(m_nuclei[2]));
		SetSystemEquation();
		return;
	}
	
	void DecaySystem::SetSystemEquation()
	{
		std::stringstream stream;
		stream << m_nuclei[0].isotopicSymbol << "->"
			   << m_nuclei[1].isotopicSymbol << "+"
			   << m_nuclei[2].isotopicSymbol;
		m_sysEquation = stream.str();
	}
	
	void DecaySystem::RunSystem()
	{
		double rxnTheta = std::acos(RandomGenerator::GetUniformReal(s_cosThetaMin, s_cosThetaMax));
		double rxnPhi = RandomGenerator::GetUniformReal(s_phiMin, s_phiMax);
		double ex = RandomGenerator::GetNormal(m_params.stepParams[0].meanResidualEx,
											   m_params.stepParams[0].sigmaResidualEx);

		m_step1.SetPolarRxnAngle(rxnTheta);
		m_step1.SetAzimRxnAngle(rxnPhi);
		m_step1.SetExcitation(ex);
		m_step1.Calculate();
	}

}