#ifndef SIM_BASE_H
#define SIM_BASE_H

#include <cassert>
#include <string>

#define ASIM_ASSERT(expr, message) assert(expr && message)

#include "MassLookup.h"
#include "Dict/Nucleus.h"

namespace AnasenSim {

    static Nucleus CreateNucleus(uint32_t z, uint32_t a, Nucleus::ReactionRole role)
    {
        Nucleus nuc;
        nuc.Z = z;
        nuc.A = a;
        nuc.groundStateMass = MassLookup::GetInstance().FindMass(z, a);
        nuc.isotopicSymbol = MassLookup::GetInstance().FindSymbol(z, a);
        nuc.vec4 = ROOT::Math::PxPyPzEVector(0., 0., 0., nuc.groundStateMass);
        nuc.role = role;
        return nuc;
    }

    static std::string ReactionRoleToString(Nucleus::ReactionRole role)
    {
        switch(role)
        {
            case Nucleus::ReactionRole::Target: return "Target";
            case Nucleus::ReactionRole::Projectile: return "Projectile";
            case Nucleus::ReactionRole::Ejectile: return "Ejectile";
            case Nucleus::ReactionRole::Residual: return "Residual";
            case Nucleus::ReactionRole::Breakup1: return "Breakup1";
            case Nucleus::ReactionRole::Breakup2: return "Breakup2";
            case Nucleus::ReactionRole::None: return "None";
        }

        return "None";
    }
}

#endif