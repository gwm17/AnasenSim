#include "QQQDetector.h"

namespace AnasenSim {

	QQQDetector::QQQDetector(double phiCentral, double zOffset, double xOffset, double yOffset) :
		m_centralPhi(phiCentral), m_translation(xOffset,yOffset,zOffset), m_norm(0.0,0.0,1.0), m_isSmearing(false)
	{
		m_zRotation.SetAngle(m_centralPhi);
		m_ringCoords.resize(s_nRings);
		m_wedgeCoords.resize(s_nWedges);
		for(auto& ring : m_ringCoords)
			ring.resize(4);
		for(auto& wedge : m_wedgeCoords)
			wedge.resize(4);

		CalculateCorners();
	}

	QQQDetector::~QQQDetector() {}

	void QQQDetector::CalculateCorners()
	{
		double x0, x1, x2, x3;
		double y0, y1, y2, y3;
		double z0, z1, z2, z3;

		//Generate flat ring corner coordinates
		for(int i=0; i<s_nRings; i++)
		{
			x0 = (s_innerR + s_deltaR*(i+1))*std::cos(-s_deltaPhiTotal/2.0);
			y0 = (s_innerR + s_deltaR*(i+1))*std::sin(-s_deltaPhiTotal/2.0);
			z0 = 0.0;
			m_ringCoords[i][0].SetXYZ(x0, y0, z0);

			x1 = (s_innerR + s_deltaR*(i))*std::cos(-s_deltaPhiTotal/2.0);
			y1 = (s_innerR + s_deltaR*(i))*std::sin(-s_deltaPhiTotal/2.0);
			z1 = 0.0;
			m_ringCoords[i][1].SetXYZ(x1, y1, z1);

			x2 = (s_innerR + s_deltaR*(i))*std::cos(s_deltaPhiTotal/2.0);
			y2 = (s_innerR + s_deltaR*(i))*std::sin(s_deltaPhiTotal/2.0);
			z2 = 0.0;
			m_ringCoords[i][2].SetXYZ(x2, y2, z2);

			x3 = (s_innerR + s_deltaR*(i+1))*std::cos(s_deltaPhiTotal/2.0);
			y3 = (s_innerR + s_deltaR*(i+1))*std::sin(s_deltaPhiTotal/2.0);
			z3 = 0.0;
			m_ringCoords[i][3].SetXYZ(x3, y3, z3);
		}

		//Generate flat wedge corner coordinates
		for(int i=0; i<s_nWedges; i++)
		{
			x0 = s_outerR * std::cos(-s_deltaPhiTotal/2.0 + i*s_deltaPhi);
			y0 = s_outerR * std::sin(-s_deltaPhiTotal/2.0 + i*s_deltaPhi);
			z0 = 0.0;
			m_wedgeCoords[i][0].SetXYZ(x0, y0, z0);

			x1 = s_innerR * std::cos(-s_deltaPhiTotal/2.0 + i*s_deltaPhi);
			y1 = s_innerR * std::sin(-s_deltaPhiTotal/2.0 + i*s_deltaPhi);
			z1 = 0.0;
			m_wedgeCoords[i][1].SetXYZ(x1, y1, z1);

			x2 = s_innerR * std::cos(-s_deltaPhiTotal/2.0 + (i+1)*s_deltaPhi);
			y2 = s_innerR * std::sin(-s_deltaPhiTotal/2.0 + (i+1)*s_deltaPhi);
			z2 = 0.0;
			m_wedgeCoords[i][2].SetXYZ(x2, y2, z2);

			x3 = s_outerR * std::cos(-s_deltaPhiTotal/2.0 + (i+1)*s_deltaPhi);
			y3 = s_outerR * std::sin(-s_deltaPhiTotal/2.0 + (i+1)*s_deltaPhi);
			z3 = 0.0;
			m_wedgeCoords[i][3].SetXYZ(x3, y3, z3);
		}

		for(int i=0; i<s_nRings; i++)
		{
			for(int j=0; j<4; j++)
				m_ringCoords[i][j] = TransformCoordinates(m_ringCoords[i][j]);
		}

		for(int i=0; i<s_nWedges; i++)
		{
			for(int j=0; j<4; j++)
				m_wedgeCoords[i][j] = TransformCoordinates(m_wedgeCoords[i][j]);
		}

	}

	ROOT::Math::XYZPoint QQQDetector::GetTrajectoryCoordinates(const ROOT::Math::XYZPoint& rxnPoint, double theta, double phi)
	{
		double z_to_detector = m_translation.Vect().Z() - rxnPoint.Z();
		double rho_traj = z_to_detector*std::tan(theta);
		if(rho_traj == 0.0)
			rho_traj = rxnPoint.Rho();
		double r_traj = std::sqrt(rho_traj*rho_traj + z_to_detector*z_to_detector);
		double min_rho, max_rho, min_phi, max_phi;

		ROOT::Math::XYZPoint result;

		for(auto& ring : m_ringCoords)
		{
			min_rho = ring[1].Rho();
			max_rho = ring[0].Rho();
			if(Precision::IsFloatLessOrAlmostEqual(rho_traj, max_rho, s_epsilon) &&
			   Precision::IsFloatGreaterOrAlmostEqual(rho_traj, min_rho, s_epsilon))
			{
				for(auto& wedge : m_wedgeCoords)
				{
					min_phi = wedge[0].Phi();
					max_phi = wedge[3].Phi();
					if(Precision::IsFloatGreaterOrAlmostEqual(phi, min_phi, s_epsilon) && 
					   Precision::IsFloatLessOrAlmostEqual(phi, max_phi, s_epsilon))
					{
						result.SetXYZ(std::sin(theta)*std::cos(phi)*r_traj, 
									  std::sin(theta)*std::sin(phi)*r_traj, 
									  std::cos(theta)*r_traj);
						break;
					}
				}
			}
		}

		return result;
	}

	std::pair<int,int> QQQDetector::GetTrajectoryRingWedge(const ROOT::Math::XYZPoint& rxnPoint, double theta, double phi)
	{
		double z_to_detector = m_translation.Vect().Z() - rxnPoint.Z();
		double rho_traj = z_to_detector*std::tan(theta);
		if(rho_traj == 0.0)
			rho_traj = rxnPoint.Rho();
		double min_rho, max_rho, min_phi, max_phi;


		for(int r=0; r<s_nRings; r++)
		{
			auto& ring = m_ringCoords[r];
			min_rho = ring[1].Rho();
			max_rho = ring[0].Rho();
			if(Precision::IsFloatLessOrAlmostEqual(rho_traj, max_rho, s_epsilon) &&
			   Precision::IsFloatGreaterOrAlmostEqual(rho_traj, min_rho, s_epsilon))
			{
				for(int w=0; w<s_nWedges; w++)
				{
					auto& wedge = m_wedgeCoords[w];
					min_phi = wedge[0].Phi();
					max_phi = wedge[3].Phi();
					if(Precision::IsFloatGreaterOrAlmostEqual(phi, min_phi, s_epsilon) && 
					   Precision::IsFloatLessOrAlmostEqual(phi, max_phi, s_epsilon))
						return std::make_pair(r, w);
				}
			}
		}

		return std::make_pair(-1, -1);
	}

	ROOT::Math::XYZPoint QQQDetector::GetHitCoordinates(int ringch, int wedgech)
	{
		if(!CheckChannel(ringch) || !CheckChannel(wedgech))
			return ROOT::Math::XYZPoint();

		double r_center, phi_center;
		if(m_isSmearing)
		{
			r_center  = s_innerR + (RandomGenerator::GetUniformFraction()+ringch)*s_deltaR;
			phi_center = -s_deltaPhiTotal/2.0 + (RandomGenerator::GetUniformFraction()+wedgech)*s_deltaPhi;
		}
		else
		{
			r_center  = s_innerR + (0.5+ringch)*s_deltaR;
			phi_center = -s_deltaPhiTotal/2.0 + (0.5+wedgech)*s_deltaPhi;
		}
		double x = r_center*std::cos(phi_center);
		double y = r_center*std::sin(phi_center);
		double z = 0;

		ROOT::Math::XYZPoint hit(x, y, z);

		return TransformCoordinates(hit);
	}

}