#include "OneStepSystem.h"
#include "RandomGenerator.h"

#include <sstream>

namespace AnasenSim {

	OneStepSystem::OneStepSystem(const SystemParameters& params) :
		ReactionSystem(params.target), m_rxnBeamEnergy(params.rxnBeamEnergy)
	{
		int zp = params.stepParams[0].Z[1];
		int ap = params.stepParams[0].A[1];
		m_rxnPathLength = m_target.GetPathLength(zp, ap, params.initialBeamEnergy, m_rxnBeamEnergy);
		m_beamStraggling = m_target.GetAngularStraggling(zp, ap, params.initialBeamEnergy, m_rxnPathLength);
		m_nuclei.resize(4);
		Init(params.stepParams);
	}
	
	OneStepSystem::~OneStepSystem() {}
	
	void OneStepSystem::Init(const std::vector<StepParameters>& params)
	{
		if(params.size() != 1 || params[0].rxnType != RxnType::Reaction ||
		   params[0].Z.size() != 3 || params[0].A.size() != 3)
		{
			m_isValid = false;
			std::cerr << "Invalid parameters at OneStepSystem::Init(), does not match OneStep signature!" << std::endl;
			return;
		}

		const StepParameters& step1Params = params[0];

		//Set nuclei

		int zr = step1Params.Z[0] + step1Params.Z[1] - step1Params.Z[2];
		int ar = step1Params.A[0] + step1Params.A[1] - step1Params.A[2];

		m_nuclei[0] = CreateNucleus(step1Params.Z[0], step1Params.A[0]); //target
		m_nuclei[1] = CreateNucleus(step1Params.Z[1], step1Params.A[1]); //projectile
		m_nuclei[2] = CreateNucleus(step1Params.Z[2], step1Params.A[2]); //ejectile
		m_nuclei[3] = CreateNucleus(zr, ar); //residual

		m_step1.BindNuclei(&(m_nuclei[0]), &(m_nuclei[1]), &(m_nuclei[2]), &(m_nuclei[3]));
		SetSystemEquation();

		//Set sampling parameters

		AddExcitationDistribution(step1Params.meanResidualEx, step1Params.sigmaResidualEx);
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
		std::mt19937& gen = RandomGenerator::GetInstance().GetGenerator();
		rxnTheta = std::acos(m_cosThetaDist(gen));
		rxnPhi = m_phiDist(gen);
		residEx = (m_exDistributions[0])(gen);
		std::normal_distribution<double> beamDist(0.0, m_beamStraggling);
		beamTheta = beamDist(gen);
		beamPhi = m_phiDist(gen);
		
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
	}

}