/*ayer
	Reaction.cpp
	Reaction is a class which implements either a decay or scattering reaction. As such it requires either
	3 (decay) or 4 (scattering) nuclei to perform any calcualtions. I also links together the target, which provides
	energy loss calculations, with the kinematics. Note that Reaction does not own the LayeredTarget.

	--GWM Jan. 2021
*/
#include "Reaction.h"
#include "SimBase.h"

#include "Math/Boost.h"

namespace AnasenSim {

	Reaction::Reaction() :
		m_target(nullptr), m_projectile(nullptr), m_ejectile(nullptr), m_residual(nullptr),
		m_bke(0), m_theta(0), m_phi(0), m_ex(0), m_beamTheta(0), m_beamPhi(0), m_rxnLayer(0), m_isInit(false)
	{
	}
	
	Reaction::Reaction(Nucleus* target, Nucleus* projectile, Nucleus* ejectile, Nucleus* residual) :
		m_target(nullptr), m_projectile(nullptr), m_ejectile(nullptr), m_residual(nullptr),
		m_bke(0), m_theta(0), m_phi(0), m_ex(0), m_rxnLayer(0)
	{
		BindNuclei(target, projectile, ejectile, residual);
	}
	
	Reaction::~Reaction()
	{
	}
	
	bool Reaction::Calculate()
	{
		if(!m_isInit) 
			return false;
	
		if(m_isDecay) {
			CalculateDecay();
			return true;
		} else {
			CalculateReaction();
			return true;
		}
	}
	
	void Reaction::BindNuclei(Nucleus* target, Nucleus* projectile, Nucleus* ejectile, Nucleus* residual)
	{
		m_target = target;
		m_projectile = projectile;
		m_ejectile = ejectile;
		m_residual = residual;

		if(m_projectile == nullptr)
			m_isDecay = true;
		else
			m_isDecay = false;

		if(m_target == nullptr || m_ejectile ==  nullptr || m_residual == nullptr)
			m_isInit = false;
		else
			m_isInit = true;
	}
	
	
	//Methods given by Iliadis in Nuclear Physics of Stars, Appendix C
	//For use with lab frame restricted angles. May not give appropriate disribution for ejectile
	void Reaction::CalculateReactionThetaLab()
	{
		m_target->vec4.SetPxPyPzE(0.,0.,0.,m_target->groundStateMass);
		double beamP = std::sqrt(m_bke*(m_bke + 2.0 * m_projectile->groundStateMass));
		double beamE = m_bke + m_projectile->groundStateMass;
		m_projectile->SetVec4Spherical(m_beamTheta, m_beamPhi, beamP, beamE);
	
		double Q = m_target->groundStateMass + m_projectile->groundStateMass - (m_ejectile->groundStateMass + m_residual->groundStateMass + m_ex);
	
		double Ethresh = -Q*(m_ejectile->groundStateMass+m_residual->groundStateMass) / 
						 (m_ejectile->groundStateMass + m_residual->groundStateMass - m_projectile->groundStateMass);

		ASIM_ASSERT(m_bke > Ethresh, "Reaction energy not above threshold");
	
		double term1 = std::sqrt(m_projectile->groundStateMass * m_ejectile->groundStateMass * m_bke)/
					   (m_ejectile->groundStateMass + m_residual->groundStateMass) * std::cos(m_theta);
		double term2 = (m_bke * (m_residual->groundStateMass - m_projectile->groundStateMass) + m_residual->groundStateMass*Q) / 
					   (m_residual->groundStateMass + m_ejectile->groundStateMass);
		double sqrt_pos_ejectKE = term1 + std::sqrt(term1*term1 + term2);
		double sqrt_neg_ejectKE = term1 - std::sqrt(term1*term1 + term2);
		double ejectKE;
		if(sqrt_pos_ejectKE > 0)
			ejectKE = sqrt_pos_ejectKE*sqrt_pos_ejectKE;
		else
			ejectKE = sqrt_neg_ejectKE*sqrt_neg_ejectKE;

		double ejectP = std::sqrt(ejectKE * (ejectKE + 2.0 * m_ejectile->groundStateMass));
		double ejectE = ejectKE + m_ejectile->groundStateMass;
	
		m_ejectile->SetVec4Spherical(m_theta, m_phi, ejectP, ejectE);
	
		m_residual->vec4 = m_target->vec4 + m_projectile->vec4 - m_ejectile->vec4;
	
		ejectP = std::sqrt(ejectKE*(ejectKE + 2.0 * m_ejectile->groundStateMass));
		ejectE = ejectKE + m_ejectile->groundStateMass;
		m_ejectile->SetVec4Spherical(m_theta, m_phi, ejectP, ejectE);
	
	}
	
	//Methods from original ANASEN. Gives proper distribution for inverse kinematics.
	void Reaction::CalculateReactionThetaCM()
	{
		//Target assumed at rest, with 0 excitation energy
		m_target->vec4.SetPxPyPzE(0.,0.,0.,m_target->groundStateMass);
		double beamP = std::sqrt(m_bke*(m_bke + 2.0 * m_projectile->groundStateMass));
		double beamE = m_bke + m_projectile->groundStateMass;
		m_projectile->SetVec4Spherical(m_beamTheta, m_beamPhi, beamP, beamE);
	
		double Q = m_target->groundStateMass + m_projectile->groundStateMass - (m_ejectile->groundStateMass + m_residual->groundStateMass + m_ex);
	
		double Ethresh = -Q*(m_ejectile->groundStateMass + m_residual->groundStateMass) / 
						 (m_ejectile->groundStateMass + m_residual->groundStateMass - m_projectile->groundStateMass);
		ASIM_ASSERT(m_bke > Ethresh, "Reaction energy not above threshold");
		
		
		ROOT::Math::PxPyPzEVector parent = m_target->vec4 + m_projectile->vec4;
		ROOT::Math::Boost boost(parent.BoostToCM());
		parent = boost*parent;
		double ejectE_cm = (std::pow(m_ejectile->groundStateMass, 2.0) - 
							std::pow(m_residual->groundStateMass + m_ex, 2.0) + std::pow(parent.E(),2.0))/
							(2.0*parent.E());
		double ejectP_cm = std::sqrt(ejectE_cm*ejectE_cm - std::pow(m_ejectile->groundStateMass, 2.0));
		m_ejectile->SetVec4Spherical(m_theta, m_phi, ejectP_cm, ejectE_cm);
		m_ejectile->vec4 = boost.Inverse() * m_ejectile->vec4;
		m_residual->vec4 = m_target->vec4 + m_projectile->vec4 - m_ejectile->vec4;
	
		double ejectKE = m_ejectile->GetKE();
		double ejectP = m_ejectile->vec4.P();
		double ejectE = m_ejectile->vec4.E();
		//energy loss for ejectile (after reaction!)
		ejectP = std::sqrt(ejectKE*(ejectKE + 2.0 * m_ejectile->groundStateMass));
		ejectE = ejectKE + m_ejectile->groundStateMass;
		m_ejectile->SetVec4Spherical(m_ejectile->vec4.Theta(), m_ejectile->vec4.Phi(), ejectP, ejectE);
	}
	
	void Reaction::CalculateReaction()
	{
		return CalculateReactionThetaCM();
	}
	
	//Calculate in CM, where decay is isotropic
	void Reaction::CalculateDecay()
	{
		double residualMass = m_residual->groundStateMass + m_ex;
		double Q = m_target->vec4.M() - m_ejectile->groundStateMass - residualMass;
		ASIM_ASSERT(Q > 0, "Decay not above threshold");
	
		ROOT::Math::Boost boost(m_target->vec4.BoostToCM());
		m_target->vec4 = boost*m_target->vec4;
		double ejectE_cm = (m_ejectile->groundStateMass*m_ejectile->groundStateMass - 
						   residualMass*residualMass + m_target->vec4.E()*m_target->vec4.E()) /
					       (2.0*m_target->vec4.E());
		double ejectP_cm = std::sqrt(ejectE_cm*ejectE_cm - m_ejectile->groundStateMass*m_ejectile->groundStateMass);
	
		m_ejectile->SetVec4Spherical(m_theta, m_phi, ejectP_cm, ejectE_cm);
		m_ejectile->thetaCM = m_theta;
	
		m_target->vec4 = boost.Inverse() * m_target->vec4;
		m_ejectile->vec4 = boost.Inverse() * m_ejectile->vec4;
	
		m_residual->vec4 = m_target->vec4 - m_ejectile->vec4;
	}

}



