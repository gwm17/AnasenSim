#include "OneStepSystem.h"
#include "RandomGenerator.h"

#include <sstream>

namespace AnasenSim {

	OneStepSystem::OneStepSystem(const SystemParameters& params) :
		ReactionSystem(params)
	{
		Init();
	}
	
	OneStepSystem::~OneStepSystem() {}
	
	void OneStepSystem::Init()
	{
		if(m_params.stepParams.size() != 1 || m_params.stepParams[0].rxnType != RxnType::Reaction ||
		   m_params.stepParams[0].Z.size() != 3 || m_params.stepParams[0].A.size() != 3)
		{
			m_isValid = false;
			std::cerr << "Invalid parameters at OneStepSystem::Init(), does not match OneStep signature!" << std::endl;
			return;
		}

		const StepParameters& step1Params = m_params.stepParams[0];

		//Set nuclei

		int zr = step1Params.Z[0] + step1Params.Z[1] - step1Params.Z[2];
		int ar = step1Params.A[0] + step1Params.A[1] - step1Params.A[2];

		m_nuclei.resize(4);
		m_nuclei[0] = CreateNucleus(step1Params.Z[0], step1Params.A[0], Nucleus::ReactionRole::Target); //target
		m_nuclei[1] = CreateNucleus(step1Params.Z[1], step1Params.A[1], Nucleus::ReactionRole::Projectile); //projectile
		m_nuclei[2] = CreateNucleus(step1Params.Z[2], step1Params.A[2], Nucleus::ReactionRole::Ejectile); //ejectile
		m_nuclei[3] = CreateNucleus(zr, ar, Nucleus::ReactionRole::Residual); //residual

		m_step1.BindNuclei(&(m_nuclei[0]), &(m_nuclei[1]), &(m_nuclei[2]), &(m_nuclei[3]));
		SetSystemEquation();

		if(!m_params.sampleBeam)
		{
			m_rxnBeamEnergy = m_params.rxnBeamEnergy;
			m_rxnPathLength = m_params.target.GetPathLength(m_nuclei[1].Z, m_nuclei[1].A, m_params.initialBeamEnergy, m_rxnBeamEnergy);
			m_beamStraggling = m_params.target.GetAngularStraggling(m_nuclei[1].Z, m_nuclei[1].A, m_params.initialBeamEnergy, m_rxnPathLength);
		}

	}
	
	void OneStepSystem::SetSystemEquation()
	{
		std::stringstream stream;
		stream << m_nuclei[0].isotopicSymbol << "("
			   << m_nuclei[1].isotopicSymbol << ", "
			   << m_nuclei[2].isotopicSymbol << ")"
			   << m_nuclei[3].isotopicSymbol;
		m_sysEquation = stream.str();
	}
	
	void OneStepSystem::RunSystem()
	{
		static double rxnTheta;
		static double rxnPhi;
		static double residEx;
		static double beamTheta;
		static double beamPhi;
		static ROOT::Math::XYZPoint rxnPoint;

		//Sample parameters
		rxnTheta = std::acos(RandomGenerator::GetUniformReal(s_cosThetaMin, s_cosThetaMax));
		rxnPhi = RandomGenerator::GetUniformReal(s_phiMin, s_phiMax);
		residEx = RandomGenerator::GetNormal(m_params.stepParams[0].meanResidualEx, m_params.stepParams[0].sigmaResidualEx);
		if(m_params.sampleBeam)
		{
			m_rxnBeamEnergy = RandomGenerator::GetUniformReal(0.0, m_params.initialBeamEnergy);
			m_rxnPathLength = m_params.target.GetPathLength(m_nuclei[1].Z, m_nuclei[1].A, m_params.initialBeamEnergy, m_params.rxnBeamEnergy);
			m_beamStraggling = m_params.target.GetAngularStraggling(m_nuclei[1].Z, m_nuclei[1].A, m_params.initialBeamEnergy, m_rxnPathLength);
		}
		beamTheta = RandomGenerator::GetUniformReal(0.0, m_beamStraggling);
		beamPhi = RandomGenerator::GetUniformReal(s_phiMin, s_phiMax);
		rxnPoint.SetXYZ(std::sin(beamTheta)*std::cos(beamPhi)*m_rxnPathLength,
						std::sin(beamTheta)*std::sin(beamPhi)*m_rxnPathLength,
						std::cos(beamTheta)*m_rxnPathLength);


		m_step1.SetPolarRxnAngle(rxnTheta);
		m_step1.SetAzimRxnAngle(rxnPhi);
		m_step1.SetExcitation(residEx);
		m_step1.SetBeamKE(m_rxnBeamEnergy);
		m_step1.SetBeamTheta(beamTheta);
		m_step1.SetBeamPhi(beamPhi);
		
		m_step1.Calculate();

		for(auto& nucleus : m_nuclei)
			nucleus.rxnPoint = rxnPoint;
	}

}