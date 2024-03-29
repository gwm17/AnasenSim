#include "TwoStepSystem.h"
#include "RandomGenerator.h"

#include <sstream>

namespace AnasenSim {
	
	TwoStepSystem::TwoStepSystem(const SystemParameters& params) :
		ReactionSystem(params)
	{
		int zp = m_params.stepParams[0].Z[1];
		int ap = m_params.stepParams[0].A[1];
		
		Init();
	}
	
	TwoStepSystem::~TwoStepSystem() {}
	
	void TwoStepSystem::Init()
	{
		if(m_params.stepParams.size() != 2 || m_params.stepParams[0].rxnType != RxnType::Reaction || m_params.stepParams[1].rxnType != RxnType::Decay ||
		   m_params.stepParams[0].Z.size() != 3 || m_params.stepParams[0].A.size() != 3 || m_params.stepParams[1].Z.size() != 2 || 
		   m_params.stepParams[1].A.size() != 2)
		{
			m_isValid = false;
			std::cerr << "Invalid parameters at TwoStepSystem::Init(), does not match TwoStep signature!" << std::endl;
			return;
		}

		const StepParameters& step1Params = m_params.stepParams[0];
		const StepParameters& step2Params = m_params.stepParams[1];

		//Setup nuclei
		int zr = step1Params.Z[0] + step1Params.Z[1] - step1Params.Z[2];
		int ar = step1Params.A[0] + step1Params.A[1] - step1Params.A[2];
		if(zr != step2Params.Z[0] || ar != step2Params.A[0])
		{
			m_isValid = false;
			std::cerr << "Invalid parameters at TwoStepSystem::Init(), step one and step two are not sequential! Step one recoil (Z,A): ("
					  << zr << "," << ar << ") Step two target (Z,A): (" << step2Params.Z[0] << "," << step2Params.A[0] << ")" <<std::endl;
			return;  
		}
		int zb = step2Params.Z[0] - step2Params.Z[1];
		int ab = step2Params.A[0] - step2Params.A[1];

		m_nuclei.resize(6);
		m_nuclei[0] = CreateNucleus(step1Params.Z[0], step1Params.A[0], Nucleus::ReactionRole::Target); //target
		m_nuclei[1] = CreateNucleus(step1Params.Z[1], step1Params.A[1], Nucleus::ReactionRole::Projectile); //projectile
		m_nuclei[2] = CreateNucleus(step1Params.Z[2], step1Params.A[2], Nucleus::ReactionRole::Ejectile); //ejectile
		m_nuclei[3] = CreateNucleus(zr, ar, Nucleus::ReactionRole::Residual); //residual
		m_nuclei[4] = CreateNucleus(step2Params.Z[1], step2Params.A[1], Nucleus::ReactionRole::Breakup1); //breakup1
		m_nuclei[5] = CreateNucleus(zb, ab, Nucleus::ReactionRole::Breakup2); //breakup2

		m_step1.BindNuclei(&(m_nuclei[0]), &(m_nuclei[1]), &(m_nuclei[2]), &(m_nuclei[3]));
		m_step2.BindNuclei(&(m_nuclei[3]), nullptr, &(m_nuclei[4]), &(m_nuclei[5]));
		SetSystemEquation();

		if(!m_params.sampleBeam)
		{
			m_rxnBeamEnergy = m_params.rxnBeamEnergy;
			m_rxnPathLength = m_params.target.GetPathLength(m_nuclei[1].Z, m_nuclei[1].A, m_params.initialBeamEnergy, m_rxnBeamEnergy);
			m_beamStraggling = m_params.target.GetAngularStraggling(m_nuclei[1].Z, m_nuclei[1].A, m_params.initialBeamEnergy, m_rxnPathLength);
		}
	}
	
	void TwoStepSystem::SetSystemEquation()
	{
		std::stringstream stream;
		stream << m_nuclei[0].isotopicSymbol << "("
			   << m_nuclei[1].isotopicSymbol << ", "
			   << m_nuclei[2].isotopicSymbol << ")"
			   << m_nuclei[3].isotopicSymbol << "->"
			   << m_nuclei[4].isotopicSymbol << "+"
			   << m_nuclei[5].isotopicSymbol;
		m_sysEquation = stream.str();
	}

	void TwoStepSystem::SampleParameters()
	{
		m_rxnTheta = std::acos(RandomGenerator::GetUniformReal(s_cosThetaMin, s_cosThetaMax));
		m_rxnPhi = RandomGenerator::GetUniformReal(s_phiMin, s_phiMax);
		m_decay1Theta = std::acos(RandomGenerator::GetUniformReal(s_cosThetaMin, s_cosThetaMax));
		m_decay1Phi = RandomGenerator::GetUniformReal(s_phiMin, s_phiMax);
		m_residEx = RandomGenerator::GetNormal(m_params.stepParams[0].meanResidualEx, m_params.stepParams[0].sigmaResidualEx);
		m_decay2Ex = RandomGenerator::GetNormal(m_params.stepParams[1].meanResidualEx, m_params.stepParams[1].sigmaResidualEx);
		//m_beamTheta = RandomGenerator::GetUniformReal(0.0, m_beamStraggling);
		m_beamPhi = RandomGenerator::GetUniformReal(s_phiMin, s_phiMax);
		if(m_params.sampleBeam)
		{
			m_rxnBeamEnergy = RandomGenerator::GetUniformReal(0.0, m_params.initialBeamEnergy);
			m_rxnPathLength = m_params.target.GetPathLength(m_nuclei[1].Z, m_nuclei[1].A, m_params.initialBeamEnergy, m_rxnBeamEnergy);
			m_beamStraggling = m_params.target.GetAngularStraggling(m_nuclei[1].Z, m_nuclei[1].A, m_params.initialBeamEnergy, m_rxnPathLength);
			m_beamTheta = RandomGenerator::GetUniformReal(0.0, m_beamStraggling);
		}
		//Testing against Nabin
		//m_residEx = m_step1.SampleExcitationPhaseSpace(m_rxnBeamEnergy, m_beamTheta, m_beamPhi, m_rxnTheta, m_rxnPhi);
	}

	void TwoStepSystem::RunSystem()
	{
		
		static ROOT::Math::XYZPoint rxnPoint;

		SampleParameters();
		//Check to make sure that the sampled configuration is valid (energy is conserved)
		while(!(m_step1.CheckReactionThreshold(m_rxnBeamEnergy, m_residEx) && m_step2.CheckDecayThreshold(m_residEx, m_decay2Ex)))
		{
			SampleParameters();
		}

		rxnPoint.SetXYZ(std::sin(m_beamTheta)*std::cos(m_beamPhi)*m_rxnPathLength,
						std::sin(m_beamTheta)*std::sin(m_beamPhi)*m_rxnPathLength,
						std::cos(m_beamTheta)*m_rxnPathLength);
	
		m_step1.SetPolarRxnAngle(m_rxnTheta);
		m_step1.SetAzimRxnAngle(m_rxnPhi);
		m_step1.SetExcitation(m_residEx);
		m_step1.SetBeamKE(m_rxnBeamEnergy);
		m_step1.SetBeamTheta(m_beamTheta);
		m_step1.SetBeamPhi(m_beamPhi);
	
		m_step2.SetPolarRxnAngle(m_decay1Theta);
		m_step2.SetAzimRxnAngle(m_decay1Phi);
		m_step2.SetExcitation(m_decay2Ex);
		
		m_step1.Calculate();
		m_step2.Calculate();

		for(auto& nucleus : m_nuclei)
			nucleus.rxnPoint = rxnPoint;
	
	}

}
