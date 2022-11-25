#include "SX3Detector.h"

/*
  Corner layout for each strip in the un-rotated frame
  0--------------------------1
  |                          |              
  |                          |              
  |                          |              
  |                          |              
  |                          |              x
  2--------------------------3    z<--------X
											|
											|
											|
											y
*/

namespace AnasenSim {

	SX3Detector::SX3Detector(double centerPhi, double centerZ, double centerRho) :
		m_centerPhi(centerPhi), m_centerZ(centerZ), m_centerRho(centerRho), m_norm(1.0,0.0,0.0), m_isSmearing(false)
	{
		m_zRotation.SetAngle(m_centerPhi);

		m_frontStripCoords.resize(s_nStrips);
		m_backStripCoords.resize(s_nStrips);
		m_rotFrontStripCoords.resize(s_nStrips);
		m_rotBackStripCoords.resize(s_nStrips);
		for(int i=0; i<s_nStrips; i++)
		{
			m_frontStripCoords[i].resize(s_nCorners);
			m_backStripCoords[i].resize(s_nCorners);
			m_rotFrontStripCoords[i].resize(s_nCorners);
			m_rotBackStripCoords[i].resize(s_nCorners);
		}
		CalculateCorners();

	}

	SX3Detector::~SX3Detector() {}

	void SX3Detector::CalculateCorners()
	{
		double y_min, y_max, z_min, z_max; 
		for (int s=0; s<s_nStrips; s++)
		{
			y_max = s_totalWidth/2.0 - s_frontStripWidth*s;
			y_min = s_totalWidth/2.0 - s_frontStripWidth*(s+1);
			z_max = m_centerZ + s_totalLength/2.0;
			z_min = m_centerZ - s_totalLength/2.0;
			m_frontStripCoords[s][2].SetXYZ(m_centerRho, y_max, z_max);
			m_frontStripCoords[s][3].SetXYZ(m_centerRho, y_max, z_min);
			m_frontStripCoords[s][0].SetXYZ(m_centerRho, y_min, z_max);
			m_frontStripCoords[s][1].SetXYZ(m_centerRho, y_min, z_min);

			z_max = (m_centerZ - s_totalLength/2.0) + (s+1)*s_backStripLength;
			z_min = (m_centerZ - s_totalLength/2.0) + (s)*s_backStripLength;
			y_max = s_totalWidth/2.0;
			y_min = -s_totalWidth/2.0;
			m_backStripCoords[s][2].SetXYZ(m_centerRho, y_max, z_max);
			m_backStripCoords[s][3].SetXYZ(m_centerRho, y_max, z_min);
			m_backStripCoords[s][0].SetXYZ(m_centerRho, y_min, z_max);
			m_backStripCoords[s][1].SetXYZ(m_centerRho, y_min, z_min);
		}

		for(int s=0; s<s_nStrips; s++)
		{
			m_rotFrontStripCoords[s][0] = m_zRotation*m_frontStripCoords[s][0];
			m_rotFrontStripCoords[s][1] = m_zRotation*m_frontStripCoords[s][1];
			m_rotFrontStripCoords[s][2] = m_zRotation*m_frontStripCoords[s][2];
			m_rotFrontStripCoords[s][3] = m_zRotation*m_frontStripCoords[s][3];

			m_rotBackStripCoords[s][0] = m_zRotation*m_backStripCoords[s][0];
			m_rotBackStripCoords[s][1] = m_zRotation*m_backStripCoords[s][1];
			m_rotBackStripCoords[s][2] = m_zRotation*m_backStripCoords[s][2];
			m_rotBackStripCoords[s][3] = m_zRotation*m_backStripCoords[s][3];
		}
	}

	ROOT::Math::XYZPoint SX3Detector::GetHitCoordinates(int front_stripch, double front_strip_ratio)
	{

		if (!ValidChannel(front_stripch) || !ValidRatio(front_strip_ratio))
			return ROOT::Math::XYZPoint(0,0,0);

		double y;
		if(m_isSmearing)
			y = -s_totalWidth/2.0 + (front_stripch + RandomGenerator::GetUniformFraction())*s_frontStripWidth;
		else
			y = -s_totalWidth/2.0 + (front_stripch+0.5)*s_frontStripWidth;

		//recall we're still assuming phi=0 det:
		ROOT::Math::XYZPoint coords(m_centerRho, y, front_strip_ratio*(s_totalLength/2) + m_centerZ);
	
		//NOW rotate by appropriate phi
		return m_zRotation*coords;

	}

	//Modified for gas target
	SX3Hit SX3Detector::GetChannelRatio(const ROOT::Math::XYZPoint& rxnPoint, double theta, double phi)
	{														
		static ROOT::Math::XYZPoint& corner = m_rotFrontStripCoords[0][0]; //Top left
		//Plane normal in rotated frame
		static ROOT::Math::XYZVector normPlane = GetNormRotated();

		ROOT::Math::XYZVector direction(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta));
		SX3Hit hit;
		//Scale factor
		double t = (corner.Dot(normPlane) - normPlane.Dot(rxnPoint))/(normPlane.Dot(direction));
		ROOT::Math::XYZPoint hitCoords = rxnPoint + t*direction;

		//Find strip channels in un-rotated frame
		hitCoords = m_zRotation.Inverse() * hitCoords;
		for (int s=0; s<s_nStrips; s++) {
			if (hitCoords.X() >=m_frontStripCoords[s][0].X() && hitCoords.X() <=m_frontStripCoords[s][0].X() && //Check min and max x (constant in flat)
				hitCoords.Y() >=m_frontStripCoords[s][1].Y() && hitCoords.Y() <=m_frontStripCoords[s][2].Y() && //Check min and max y
				hitCoords.Z() >=m_frontStripCoords[s][1].Z() && hitCoords.Z() <=m_frontStripCoords[s][0].Z()) //Check min and max z
			{
				hit.front_strip_index = s;
				hit.front_ratio = (hitCoords.Z()-m_centerZ)/(s_totalLength/2);
				break;
			}
		}

		for (int s=0; s<s_nStrips; s++) {
			if (hitCoords.X() >= m_backStripCoords[s][0].X() && hitCoords.X() <= m_backStripCoords[s][0].X() && //Check min and max x (constant in flat)
				hitCoords.Y() >= m_backStripCoords[s][1].Y() && hitCoords.Y() <= m_backStripCoords[s][2].Y() && //Check min and max y
				hitCoords.Z() >= m_backStripCoords[s][1].Z() && hitCoords.Z() <= m_backStripCoords[s][0].Z()) //Check min and max z
			{
				hit.back_strip_index = s;
				break;
			}
		}

		return hit;
	}

}