#include "PCDetector.h"
#include "Sim/SimBase.h"

namespace AnasenSim {

    //Works by taking rxn point and assigning PC wire which is closest to rxnPoint. Then z-position is calculated using 
    //projectile theta and phi from wire.
    ROOT::Math::XYZPoint PCDetector::AssignPCOriginal(const ROOT::Math::XYZPoint& rxnPoint, double theta, double phi, uint32_t zp)
    {
        ASIM_ASSERT(rxnPoint.Rho() < s_wireRadius, "Reaction did not occur within the PC");

        double wirePhi;
        double wirePhiRxn;
        double wireXRxn;
        double wireYRxn;
        double minPhiDiff = 100.0;
        double curPhiDiff = 0.0;
        int minWireID = 0;
        double minPhi = 0.0;
        double minPhiRxn = 0.0;
        double minRho = 0.0;
        for(int i=0; i<s_nWires; i++)
        {
            wirePhi = i*s_wireDeltaPhi;

            wireXRxn = s_wireRadius * std::cos(wirePhi) - rxnPoint.X();
            wireYRxn = s_wireRadius * std::sin(wirePhi) - rxnPoint.Y();

            wirePhiRxn = std::atan2(wireYRxn, wireXRxn);
            if(wirePhiRxn < 0.0)
                wirePhiRxn += M_PI*2.0;

            curPhiDiff = std::fabs(wirePhiRxn - phi);
            if(curPhiDiff < minPhiDiff)
            {
                minWireID = i;
                minPhi = wirePhi;
                minPhiRxn = wirePhiRxn;
                minRho = std::sqrt(wireXRxn * wireXRxn - wireYRxn * wireYRxn);
            }
        }

        double zhit = 0.0;
        if(theta == M_PI/2.0)
            zhit = rxnPoint.Z();
        else
            zhit = minRho/std::tan(theta) - rxnPoint.Z();

        if(zp == 1)
            zhit = RandomGenerator::GetNormal(zhit, s_zSigmaProton);
        else
            zhit = RandomGenerator::GetNormal(zhit, s_zSigmaNonProton);

        return ROOT::Math::XYZPoint(s_wireRadius*std::cos(wirePhi), s_wireRadius*std::sin(wirePhi), zhit);
    }

    //Find point of nearest approach by checking intersection with cylinder formed by pc wires. Then take nearest PC to intersection, keeping z constant
    //This method I think better captures the idea of the PC, plus no for loops
    ROOT::Math::XYZPoint PCDetector::AssignPC(const ROOT::Math::XYZPoint& rxnPoint, double theta, double phi, uint32_t zp)
    {
        ASIM_ASSERT(rxnPoint.Rho() < s_wireRadius, "Reaction did not occur within the PC");
        //Will neeeever intersect
        if(theta == 0.0)
        {
            return ROOT::Math::XYZPoint();
        }

        ROOT::Math::XYZVector particleTraj(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta));

        double b  = 2.0 * (rxnPoint.X() * particleTraj.X() + rxnPoint.Y() * particleTraj.Y());
        double a  = particleTraj.Rho() * particleTraj.Rho();
        double c = rxnPoint.Rho() * rxnPoint.Rho() - s_wireRadius * s_wireRadius;
        double radicand = b*b - 4.0 * a * c;
        if(radicand < 0.0)
            return ROOT::Math::XYZPoint();
        
        double t1 = (-1.0*b + std::sqrt(radicand))/(2.0 * a);
        double t2 = (-1.0*b - std::sqrt(radicand))/(2.0 * a);
        ROOT::Math::XYZPoint intersectionPoint;
        if(t1 < 0.0 && t2 < 0.0)
            return ROOT::Math::XYZPoint();
        if(t1 < 0.0)
            intersectionPoint = rxnPoint + t2 * particleTraj;
        else
            intersectionPoint = rxnPoint + t1 * particleTraj;

        double intersectPhi = intersectionPoint.Phi();
        if(intersectPhi < 0.0)
            intersectPhi += 2.0 * M_PI;
        
        double wireFrac = intersectPhi/s_wireDeltaPhi;
        double nearestWirePhi = std::round(wireFrac) * s_wireDeltaPhi;

        if(zp == 1)
            return ROOT::Math::XYZPoint(s_wireRadius * std::cos(nearestWirePhi), s_wireRadius * std::sin(nearestWirePhi), RandomGenerator::GetNormal(intersectionPoint.Z(), s_zSigmaProton));
        else
            return ROOT::Math::XYZPoint(s_wireRadius * std::cos(nearestWirePhi), s_wireRadius * std::sin(nearestWirePhi), RandomGenerator::GetNormal(intersectionPoint.Z(), s_zSigmaNonProton));

    }
}