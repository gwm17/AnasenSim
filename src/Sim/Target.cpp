/*

Target.cpp
A basic target unit for use in the SPANCRedux environment. A target
is defined as a single compound with elements Z,A of a given stoichiometry 
Holds an energy loss class

Based on code by D.W. Visser written at Yale for the original SPANC

Written by G.W. McCann Aug. 2020

*/
#include "Target.h"
#include "catima/nucdata.h"

namespace AnasenSim {

	/*Targets must be of known thickness*/
	Target::Target(const std::vector<uint32_t>& z, const std::vector<uint32_t>& a, const std::vector<int>& stoich, double density)
	{
		MassLookup& masses = MassLookup::GetInstance();
		for(size_t i=0; i<z.size(); i++)
		{
			m_material.add_element(masses.FindMassU(z[i], a[i]), z[i], stoich[i]);
		}
		m_material.density(density);
	}
	
	Target::~Target() {}
	
	/*Calculates energy loss for travelling all the way through the target*/
	double Target::GetEnergyLoss(int zp, int ap, double startEnergy, double pathLength)
	{
		catima::Projectile proj(MassLookup::GetInstance().FindMassU(zp, ap), zp, 0.0, 0.0);
		proj.T = startEnergy/proj.A;
		m_material.thickness_cm(pathLength); //Takes in a path length and calculates density corrected thickness to thickness param
		return catima::integrate_energyloss(proj, m_material);
	}

	/*Calculates reverse energy loss for travelling all the way through the target*/
	double Target::GetReverseEnergyLoss(int zp, int ap, double finalEnergy, double pathLength)
	{
		catima::Projectile proj(MassLookup::GetInstance().FindMassU(zp, ap), zp, 0.0, 0.0);
		proj.T = finalEnergy/proj.A;
		m_material.thickness_cm(pathLength);
		return catima::reverse_integrate_energyloss(proj, m_material);
	}

	//Use golden section search to find path length to reach finalEnergy
	double Target::GetPathLength(int zp, int ap, double startEnergy, double finalEnergy)
	{
		static double goldenRatioInv = 1.0/((std::sqrt(5.0) + 1.0) * 0.5);

		

		catima::Projectile proj(MassLookup::GetInstance().FindMassU(zp, ap), zp, 0.0, 0.0);
		proj.T = startEnergy/proj.A;

		//Our function to minimize, want enough energy loss that startEnergy = finalEnergy + energyLoss
		auto function = [this, startEnergy, finalEnergy](catima::Projectile& proj, double pathLength) -> double
		{
			proj.T = startEnergy/proj.A; //Reset projectile KE for energy loss algorithm
			m_material.thickness_cm(pathLength);
			return startEnergy - (finalEnergy + catima::integrate_energyloss(proj, m_material));
		};

		double stopRange = catima::range(proj, m_material); //get the total range for startEnergy -> 0
		//early out for complete stopping
		if(finalEnergy == 0.0)
			return stopRange;

		//Our search points
		double testLow, testHigh, boundLow, boundHigh;
		double valueLow, valueHigh, testValueLow, testValueHigh;
		boundLow = 0.0;
		boundHigh = stopRange;

		testHigh = boundHigh - (boundHigh - boundLow) * goldenRatioInv;
		testLow = boundLow + (boundHigh - boundLow) * goldenRatioInv;
	
		while(std::fabs(boundHigh - boundLow) > s_rangeErrorTol)
		{
			valueLow = function(proj, testLow);
			valueHigh = function(proj, testHigh);
			if(valueLow < valueHigh)
				boundHigh = testLow;
			else
				boundLow = testHigh;

			//Recalc both test points to avoid floating point errors
			testHigh = boundHigh - (boundHigh - boundLow) * goldenRatioInv;
			testLow = boundLow + (boundHigh - boundLow) * goldenRatioInv;
		}

		return (boundHigh - boundLow) * 0.5;
	}

	double Target::GetAngularStraggling(int zp, int ap, double energy, double pathLength)
	{
		catima::Projectile proj(MassLookup::GetInstance().FindMassU(zp, ap), zp, 0.0, 0.0);
		proj.T = energy/proj.A;
		m_material.thickness_cm(pathLength);

		return catima::angular_straggling(proj, m_material);
	}

}
