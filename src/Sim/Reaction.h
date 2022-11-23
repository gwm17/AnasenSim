/*
	Reaction.h
	Reaction is a class which implements either a decay or scattering reaction. As such it requires either
	3 (decay) or 4 (scattering) nuclei to perform any calcualtions. I also links together the target, which provides
	energy loss calculations, with the kinematics. Note that Reaction does not own the LayeredTarget.

	--GWM Jan. 2021
*/
#ifndef REACTION_H
#define REACTION_H

#include "Dict/Nucleus.h"

namespace AnasenSim {

	class Reaction
	{
	public:
		Reaction();
		Reaction(Nucleus* target, Nucleus* projectile, Nucleus* ejectile, Nucleus* residual);
		~Reaction();
		bool Calculate(); //do sim

		//Bind system nuclei to the specific reaction. Reaction does NOT own nuclei
		void BindNuclei(Nucleus* target, Nucleus* projectile, Nucleus* ejectile, Nucleus* residual);
		void SetBeamKE(double bke) { m_bke = bke; }
		void SetBeamTheta(double theta) { m_beamTheta = theta; }
		void SetBeamPhi(double phi) { m_beamPhi = phi; }

		void SetPolarRxnAngle(double theta) { m_theta = theta; };
		void SetAzimRxnAngle(double phi) { m_phi = phi; };
		void SetExcitation(double ex) { m_ex = ex; };

		//Can rebind individuals if needed
		void BindTarget(Nucleus* nuc) { m_target = nuc; };
		void BindProjectile(Nucleus* nuc) { m_projectile = nuc; };
		void BindEjectile(Nucleus* nuc) { m_ejectile = nuc; };
		void BindResidual(Nucleus* nuc) { m_residual = nuc; };

		bool IsDecay() const { return m_isDecay; };

	private:
		void CalculateDecay(); //target -> light_decay (eject) + heavy_decay(resid)
		void CalculateReaction(); //target + project -> eject + resid
		void CalculateReactionThetaLab();
		void CalculateReactionThetaCM();
	
		//Reactants -> NOT OWNED BY RXN
		Nucleus* m_target;
		Nucleus* m_projectile;
		Nucleus* m_ejectile;
		Nucleus* m_residual;

		double m_bke, m_theta, m_phi, m_ex;
		double m_beamTheta, m_beamPhi;
	
		int m_rxnLayer;
	
		bool m_isDecay, m_isInit;
	};

}

#endif