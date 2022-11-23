/*

Target.h
A basic target unit for use in the Mask environment. A target
is defined as a single compound with elements Z,A of a given stoichiometry
Holds an energy loss class

Based on code by D.W. Visser written at Yale for the original SPANC

Written by G.W. McCann Aug. 2020

*/
#ifndef TARGET_H
#define TARGET_H

#include <string>
#include <vector>
#include <cmath>
#include "catima/gwm_integrators.h"
#include "MassLookup.h"

namespace AnasenSim {

	class Target {
	
	public:
	 	Target() = default;
	 	Target(const std::vector<int>& z, const std::vector<int>& a, const std::vector<int>& stoich, double density);
	 	~Target();

	 	double GetEnergyLoss(int zp, int ap, double startEnergy, double pathLength);
	 	double GetReverseEnergyLoss(int zp, int ap, double finalEnergy, double pathLength);
		double GetPathLength(int zp, int ap, double startEnergy, double finalEnergy); //Returns pathlength for a particle w/ startE to reach finalE (cm)
		double GetAngularStraggling(int zp, int ap, double energy, double pathLength); //Returns planar angular straggling in radians for a particle with energy and pathLength
	 	inline double GetDensity() { return m_material.density(); } //g/cm^3
	
	private:
		catima::Material m_material;

		static constexpr double s_rangeErrorTol = 1.0e-6;
	};

}

#endif 
