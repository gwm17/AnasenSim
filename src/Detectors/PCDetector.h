#ifndef PC_DETECTOR_H
#define PC_DETECTOR_H

#include <cmath>
#include <vector>

#include "Sim/RandomGenerator.h"
#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/RotationZ.h"
#include "Math/Translation3D.h"

namespace AnasenSim {

    class PCDetector
    {
    public:
        static ROOT::Math::XYZPoint AssignPCOriginal(const ROOT::Math::XYZPoint& rxnPoint, double theta, double phi, uint32_t zp);
        static ROOT::Math::XYZPoint AssignPC(const ROOT::Math::XYZPoint& rxnPoint, double theta, double phi, uint32_t zp);

    private:
        static constexpr int s_nWires = 24;
        static constexpr double s_wireRadius = 0.03846264509; //m
        static constexpr double s_wireDeltaPhi = 2.0*M_PI/s_nWires;

        //Z-position uncertainty for protons and all others (2 cm for proton, 1cm for all else, here as sigma)
        static constexpr double s_zSigmaProton = 0.066666;
        static constexpr double s_zSigmaNonProton = 0.033333333;
    };
}

#endif