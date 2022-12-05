#ifndef ANASEN_ARRAY_H
#define ANASEN_ARRAY_H

#include <string>

#include "SX3Detector.h"
#include "QQQDetector.h"
#include "Sim/Target.h"
#include "Dict/Nucleus.h"

namespace AnasenSim {

	class AnasenArray
	{
	public:
		AnasenArray(const Target& gas);
		~AnasenArray();
		void IsDetected(Nucleus& nucleus);
		void DrawDetectorSystem(const std::string& filename);
		double RunConsistencyCheck();
		//void SetDeadChannelMap(const std::string& filename) { dmap.LoadMapfile(filename); }

	private:
		void IsBarrel1(Nucleus& nucleus);
		void IsBarrel2(Nucleus& nucleus);
		void IsQQQ(Nucleus& nucleus);

		std::vector<SX3Detector> m_barrel1;
		std::vector<SX3Detector> m_barrel2;
		std::vector<QQQDetector> m_qqq;

		Target m_detectorEloss;
		Target m_gasEloss;

		ROOT::Math::XYZPoint m_nullPoint;

		//AnasenDeadChannelMap dmap;

		/**** ANASEN geometry constants *****/
		static constexpr double s_epsilon = 1.0e-6; //accuracy
		static constexpr int s_nSX3PerBarrel = 12;
		static constexpr int s_nQQQ = 4;
		static constexpr double s_sx3Length = 0.075;
		static constexpr double s_barrelGap = 0.0254;  //Space between edge of frames of each SX3 barrel
		static constexpr double s_sx3FrameGap = 0.049; //0.049 is empty space due to width of SX3 barrel frame
		static constexpr double s_totalLength = 0.554; //total length of the ANASEN chamber
		static constexpr double s_qqqZ = s_totalLength;
		static constexpr double s_barrel1Z = s_qqqZ - (0.0125 + s_sx3Length * 0.5);
		static constexpr double s_barrel2Z = s_barrel1Z - (s_sx3Length + s_barrelGap + s_sx3FrameGap * 2.0);
		static constexpr double s_qqqZList[4] = {s_qqqZ, s_qqqZ, s_qqqZ, s_qqqZ};
		static constexpr double s_qqqPhiList[4] = {5.49779, 0.785398, 2.35619, 3.92699};
		static constexpr double s_barrelRhoList[12] = {0.0890601, 0.0889871, 0.0890354, 0.0890247, 0.0890354, 0.0890354, 0.0890247,
													 0.0890354, 0.0890354, 0.0890247, 0.0890354, 0.0890354};
		static constexpr double s_barrelPhiList[12] = {4.97426, 5.49739, 6.02132, 0.261868, 0.785398, 1.30893, 1.83266, 2.35619, 2.87972, 
													   3.40346, 3.92699, 4.45052};
		/*************************/

		static constexpr double s_energyThreshold = 0.6; //MeV
		static constexpr double s_deg2rad = M_PI/180.0;
		static constexpr double s_detectorDensity = 2.33; //g/cm^3, Si crystal
		static constexpr double s_detectorThickness = 0.001; //m
	};

}
#endif