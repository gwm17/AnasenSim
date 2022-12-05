/*
	QQQDetector.h
	Class implementing geometry for QQQDetector where the detector is perpendicular to the beam axis.
	Detector is first generated centered on the x-axis (phi=0)

	Coordinate convention : +z is downstream, -z is upstream. +y is vertically down in the lab.
*/
#ifndef QQQDETECTOR_H
#define QQQDETECTOR_H

#include <cmath>
#include <vector>

#include "Sim/RandomGenerator.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/RotationZ.h"
#include "Math/Translation3D.h"
#include "IsEqual.h"

namespace AnasenSim {
	class QQQDetector
	{
	public:
		QQQDetector(double phiCentral, double zOffset, double xOffset=0, double yOffset=0);
		~QQQDetector();
		const ROOT::Math::XYZPoint& GetRingCoordinates(int ringch, int corner) { return m_ringCoords[ringch][corner]; }
		const ROOT::Math::XYZPoint& GetWedgeCoordinates(int wedgech, int corner) { return m_wedgeCoords[wedgech][corner]; }
		const ROOT::Math::XYZVector& GetNorm() { return m_norm; }
		ROOT::Math::XYZPoint GetTrajectoryCoordinates(const ROOT::Math::XYZPoint& rxnPoint, double theta, double phi);
		std::pair<int, int> GetTrajectoryRingWedge(const ROOT::Math::XYZPoint& rxnPoint, double theta, double phi);
		ROOT::Math::XYZPoint GetHitCoordinates(int ringch, int wedgech);

		void SetSmearing(bool isSmearing) { m_isSmearing = isSmearing; }

		int GetNumberOfRings() { return s_nRings; }
		int GetNumberOfWedges() { return s_nWedges; }

	private:

		bool CheckChannel(int ch) { return (ch >=0 && ch < s_nRings); }
		bool CheckCorner(int corner) { return (corner >=0 && corner < 4); }

		void CalculateCorners();
		ROOT::Math::XYZPoint TransformCoordinates(ROOT::Math::XYZPoint& vector) { return m_translation * (m_zRotation * vector) ; }

		double m_centralPhi;

		std::vector<std::vector<ROOT::Math::XYZPoint>> m_ringCoords, m_wedgeCoords;
		ROOT::Math::Translation3D m_translation;
		ROOT::Math::XYZVector m_norm;
		ROOT::Math::RotationZ m_zRotation;

		bool m_isSmearing;

		static constexpr double s_epsilon = 1.0e-6; //accuracy
		static constexpr int s_nRings = 16;
		static constexpr int s_nWedges = 16;
		static constexpr double s_deg2rad  = M_PI/180.0;
		static constexpr double s_innerR = 0.0501;
		static constexpr double s_outerR = 0.0990;
		static constexpr double s_deltaR = (s_outerR - s_innerR) / s_nRings;
		static constexpr double s_deltaPhiTotal = 1.52119;
		static constexpr double s_deltaPhi = s_deltaPhiTotal / s_nWedges;
	};

}
#endif