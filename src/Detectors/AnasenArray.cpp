#include "AnasenArray.h"
#include "PCDetector.h"
#include <fstream>
#include <iomanip>
#include <iostream>

namespace AnasenSim {

	AnasenArray::AnasenArray(const Target& gas) :
		m_detectorEloss({14}, {28}, {1}, s_detectorDensity), m_gasEloss(gas), m_nullPoint(0., 0., 0.)
	{
		for(int i=0; i<s_nSX3PerBarrel; i++)
		{
			m_Ring1.emplace_back(s_barrelPhiList[i], s_barrel1Z, s_barrelRhoList[i]);
			m_Ring1[i].SetPixelSmearing(true);
			m_Ring2.emplace_back(s_barrelPhiList[i], s_barrel2Z, s_barrelRhoList[i]);
			m_Ring2[i].SetPixelSmearing(true);
		}
		for(int i=0; i<s_nQQQ; i++)
		{
			m_forwardQQQs.emplace_back(s_qqqPhiList[i], s_qqqZList[i]);
			m_forwardQQQs[i].SetSmearing(true);
		}
	}

	AnasenArray::~AnasenArray() {}


	void AnasenArray::DrawDetectorSystem(const std::string& filename)
	{
		std::ofstream output(filename);

		std::vector<double> x, y, z;
		std::vector<double> cx, cy, cz;
		ROOT::Math::XYZPoint coords;
		for(int i=0; i<s_nSX3PerBarrel; i++)
		{
			for(int j=0; j<4; j++)
			{
				for(int k=0; k<4; k++)
				{
					coords = m_Ring1[i].GetRotatedFrontStripCoordinates(j, k);
					x.push_back(coords.X());
					y.push_back(coords.Y());
					z.push_back(coords.Z());
					coords = m_Ring1[i].GetRotatedBackStripCoordinates(j, k);
					x.push_back(coords.X());
					y.push_back(coords.Y());
					z.push_back(coords.Z());
				}
				coords = m_Ring1[i].GetHitCoordinates(j, 0);
				cx.push_back(coords.X());
				cy.push_back(coords.Y());
				cz.push_back(coords.Z());
			}
		}
		for(int i=0; i<s_nSX3PerBarrel; i++)
		{
			for(int j=0; j<4; j++)
			{
				for(int k=0; k<4; k++)
				{
					coords = m_Ring2[i].GetRotatedFrontStripCoordinates(j, k);
					x.push_back(coords.X());
					y.push_back(coords.Y());
					z.push_back(coords.Z());
					coords = m_Ring2[i].GetRotatedBackStripCoordinates(j, k);
					x.push_back(coords.X());
					y.push_back(coords.Y());
					z.push_back(coords.Z());
				}
				coords = m_Ring2[i].GetHitCoordinates(j, 0);
				cx.push_back(coords.X());
				cy.push_back(coords.Y());
				cz.push_back(coords.Z());
			}
		}
		for(int i=0; i<s_nQQQ; i++)
		{
			for(int j=0; j<16; j++)
			{
				for(int k=0; k<4; k++)
				{
					coords = m_forwardQQQs[i].GetRingCoordinates(j, k);
					x.push_back(coords.X());
					y.push_back(coords.Y());
					z.push_back(coords.Z());
					coords = m_forwardQQQs[i].GetWedgeCoordinates(j, k);
					x.push_back(coords.X());
					y.push_back(coords.Y());
					z.push_back(coords.Z());
				}
				for(int k=0; k<16; k++)
				{
					coords = m_forwardQQQs[i].GetHitCoordinates(j, k);
					cx.push_back(coords.X());
					cy.push_back(coords.Y());
					cz.push_back(coords.Z());
				}
			}
		}

		output<<"ANASEN Geometry File -- Coordinates for Detectors"<<std::endl;
		for(std::size_t i=0; i<x.size(); i++)
			output<<x[i]<<" "<<y[i]<<" "<<z[i]<<std::endl;
		for(std::size_t i=0; i<cx.size(); i++)
			output<<cx[i]<<" "<<cy[i]<<" "<<cz[i]<<std::endl;

		output.close();
	}

	double AnasenArray::RunConsistencyCheck()
	{
		std::vector<ROOT::Math::XYZPoint> r1_points;
		std::vector<ROOT::Math::XYZPoint> r2_points;
		std::vector<ROOT::Math::XYZPoint> fqqq_points;
		std::vector<ROOT::Math::XYZPoint> bqqq_points;
		for(int i=0; i<s_nSX3PerBarrel; i++)
		{
			for(int j=0; j<4; j++)
				r1_points.push_back(m_Ring1[i].GetHitCoordinates(j, 0));
		}
		for(int i=0; i<s_nSX3PerBarrel; i++)
		{
			for(int j=0; j<4; j++)
				r2_points.push_back(m_Ring2[i].GetHitCoordinates(j, 0));
		}
		for(int i=0; i<s_nQQQ; i++)
		{
			for(int j=0; j<16; j++)
			{
				for(int k=0; k<16; k++)
					fqqq_points.push_back(m_forwardQQQs[i].GetHitCoordinates(j, k));
			}
		}

		std::size_t npoints = r1_points.size() + r2_points.size() + fqqq_points.size() + bqqq_points.size();
		std::size_t count = 0;
		ROOT::Math::XYZPoint coords;
		for(auto& point : r1_points)
		{
			for(auto& sx3 : m_Ring1)
			{
				auto result = sx3.GetChannelRatio({0., 0., 0.}, point.Theta(), point.Phi());
				coords = sx3.GetHitCoordinates(result.front_strip_index, result.front_ratio);
				if(IsDoubleEqual(point.X(), coords.X()) && IsDoubleEqual(point.Y(), coords.Y()) && IsDoubleEqual(point.Z(), coords.Z()))
				{
					count++;
					break;
				}
			}
		}
		for(auto& point : r2_points)
		{
			for(auto& sx3 : m_Ring2)
			{
				auto result = sx3.GetChannelRatio({0., 0., 0.}, point.Theta(), point.Phi());
				coords = sx3.GetHitCoordinates(result.front_strip_index, result.front_ratio);
				if(IsDoubleEqual(point.X(), coords.X()) && IsDoubleEqual(point.Y(), coords.Y()) && IsDoubleEqual(point.Z(), coords.Z()))
				{
					count++;
					break;
				}
			}
		}
		for(auto& point : fqqq_points)
		{
			for(auto& qqq : m_forwardQQQs)
			{
				auto result = qqq.GetTrajectoryRingWedge({0., 0., 0.}, point.Theta(), point.Phi());
				coords = qqq.GetHitCoordinates(result.first, result.second);
				if(IsDoubleEqual(point.X(), coords.X()) && IsDoubleEqual(point.Y(), coords.Y()) && IsDoubleEqual(point.Z(), coords.Z()))
				{
					count++;
					break;
				}
			}
		}

		double ratio = ((double)count)/((double)npoints);

		return ratio;

	}

	void AnasenArray::IsRing1(Nucleus& nucleus)
	{
		static double thetaIncident;
		static double effectiveThickness;
		static double energyAtSi;
		for(int i=0; i<s_nSX3PerBarrel; i++)
		{
			auto result = m_Ring1[i].GetChannelRatio(nucleus.rxnPoint, nucleus.vec4.Theta(), nucleus.vec4.Phi());
			if(result.front_strip_index != -1) 
			{
				nucleus.pcVector = PCDetector::AssignPC(nucleus.rxnPoint, nucleus.vec4.Theta(), nucleus.vec4.Phi(), nucleus.Z);
				nucleus.isDetected = true;
				nucleus.siVector = m_Ring1[i].GetHitCoordinates(result.front_strip_index, result.front_ratio);

				thetaIncident = std::acos(nucleus.siVector.Dot(m_Ring1[i].GetNormRotated())/nucleus.siVector.R());
				effectiveThickness = s_detectorThickness_cm/std::fabs(std::cos(thetaIncident));
				nucleus.pcDetE = m_gasEloss.GetEnergyLoss(nucleus.Z, nucleus.A, nucleus.GetKE(), (nucleus.pcVector - nucleus.rxnPoint).R());
				energyAtSi = nucleus.GetKE() - m_gasEloss.GetEnergyLoss(nucleus.Z, nucleus.A, nucleus.GetKE(), (nucleus.siVector - nucleus.rxnPoint).R());
				nucleus.siliconDetKE = m_detectorEloss.GetEnergyLoss(nucleus.Z, nucleus.A, energyAtSi, effectiveThickness);

				nucleus.siDetectorName = "R1";
				return;
			}
		}
	}

	void AnasenArray::IsRing2(Nucleus& nucleus)
	{
		static double thetaIncident;
		static double effectiveThickness;
		static double energyAtSi;
		for(int i=0; i<s_nSX3PerBarrel; i++)
		{
			auto result = m_Ring2[i].GetChannelRatio(nucleus.rxnPoint, nucleus.vec4.Theta(), nucleus.vec4.Phi());
			if(result.front_strip_index != -1) 
			{
				nucleus.pcVector = PCDetector::AssignPC(nucleus.rxnPoint, nucleus.vec4.Theta(), nucleus.vec4.Phi(), nucleus.Z);
				nucleus.isDetected = true;
				nucleus.siVector = m_Ring2[i].GetHitCoordinates(result.front_strip_index, result.front_ratio);

				thetaIncident = std::acos(nucleus.siVector.Dot(m_Ring2[i].GetNormRotated())/nucleus.siVector.R());
				effectiveThickness = s_detectorThickness_cm/std::fabs(std::cos(thetaIncident));
				nucleus.pcDetE = m_gasEloss.GetEnergyLoss(nucleus.Z, nucleus.A, nucleus.GetKE(), (nucleus.pcVector - nucleus.rxnPoint).R());
				energyAtSi = nucleus.GetKE() - m_gasEloss.GetEnergyLoss(nucleus.Z, nucleus.A, nucleus.GetKE(), (nucleus.siVector - nucleus.rxnPoint).R());
				nucleus.siliconDetKE = m_detectorEloss.GetEnergyLoss(nucleus.Z, nucleus.A, energyAtSi, effectiveThickness);

				nucleus.siDetectorName = "R2";
				return;
			}
		}
	}

	void AnasenArray::IsQQQ(Nucleus& nucleus)
	{
		double thetaIncident;
		double effectiveThickness;
		static double energyAtSi;
		for(int i=0; i<s_nQQQ; i++)
		{
			auto result = m_forwardQQQs[i].GetTrajectoryRingWedge(nucleus.rxnPoint, nucleus.vec4.Theta(), nucleus.vec4.Phi());
			if(result.first != -1) 
			{
				nucleus.pcVector = PCDetector::AssignPC(nucleus.rxnPoint, nucleus.vec4.Theta(), nucleus.vec4.Phi(), nucleus.Z);
				nucleus.isDetected = true;
				nucleus.siVector = m_forwardQQQs[i].GetHitCoordinates(result.first, result.second);

				thetaIncident = std::acos(nucleus.siVector.Dot(m_forwardQQQs[i].GetNorm())/nucleus.siVector.R());
				effectiveThickness = s_detectorThickness_cm / std::fabs(std::cos(thetaIncident));
				nucleus.pcDetE = m_gasEloss.GetEnergyLoss(nucleus.Z, nucleus.A, nucleus.GetKE(), (nucleus.pcVector - nucleus.rxnPoint).R());
				energyAtSi = nucleus.GetKE() - m_gasEloss.GetEnergyLoss(nucleus.Z, nucleus.A, nucleus.GetKE(), (nucleus.siVector - nucleus.rxnPoint).R());
				nucleus.siliconDetKE = m_detectorEloss.GetEnergyLoss(nucleus.Z, nucleus.A, nucleus.GetKE(), effectiveThickness);
				nucleus.siDetectorName = "FQQQ";
				return;
			}
		}
	}

	void AnasenArray::IsDetected(Nucleus& nucleus)
	{
		if(nucleus.GetKE() <= s_energyThreshold)
			return;

		if(!nucleus.isDetected)
			IsRing1(nucleus);
		if(!nucleus.isDetected)
			IsRing2(nucleus);
		if(!nucleus.isDetected)
			IsQQQ(nucleus);
	}

}