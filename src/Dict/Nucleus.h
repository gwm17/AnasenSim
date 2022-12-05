/*
	Nucleus.h
	Nucleus is a derived class of Vec4. A nucleus is the kinematics is essentially a 4 vector with the
	additional properties of the number of total nucleons (A), the number of protons (Z), a ground state mass,
	an exctitation energy, and an isotopic symbol.

	--GWM Jan 2021
*/
#ifndef NUCLEUS_H
#define NUCLEUS_H

#include <string>
#include <vector>
#include "Math/Vector4D.h"
#include "Math/Point3D.h"

namespace AnasenSim {

	struct Nucleus
	{
		enum class ReactionRole
		{
			Target = 0,
			Projectile = 1,
			Ejectile = 2,
			Residual = 3,
			Breakup1 = 4,
			Breakup2 = 5,
			None = -1
		};

		void SetVec4Spherical(double theta, double phi, double p, double E)
		{
			vec4.SetPxPyPzE(std::sin(theta)*std::cos(phi)*p,
							std::sin(theta)*std::sin(phi)*p,
							std::cos(theta)*p,
							E
						   );
		}

		double GetKE() const //MeV
		{
			return vec4.E() - vec4.M();
		}

		double GetExcitationEnergy() const //MeV
		{
			return vec4.M() - groundStateMass;
		}

		uint32_t Z = 0; 
		uint32_t A = 0;
		double groundStateMass = 0.0; //MeV
		std::string isotopicSymbol = "";
		double thetaCM = 0.0; // rad
		ROOT::Math::PxPyPzEVector vec4;
		ROOT::Math::XYZPoint rxnPoint;
		ReactionRole role = ReactionRole::None;

		bool isDetected = false;
		double siliconDetKE = 0.0; //MeV
		ROOT::Math::XYZPoint siVector;
        double pcDetE = 0.0; //MeV
		ROOT::Math::XYZPoint pcVector;
		std::string siDetectorName = "";
	};

	bool EnforceDictionaryLinked();

};

#endif
