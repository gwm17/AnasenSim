#ifndef IS_EQUAL_H
#define IS_EQUAL_H

#include <cmath>
#include <algorithm>

namespace AnasenSim {

    namespace Precision {

        static constexpr bool IsFloatAlmostEqual(double val1, double val2, double epsilon)
        {
            return std::fabs(val1 - val2) < epsilon;
        }

        static constexpr bool IsFloatLessOrAlmostEqual(double lhs, double rhs, double epsilon)
        {
            return IsFloatAlmostEqual(lhs, rhs, epsilon) || lhs < rhs;
        }

        static constexpr bool IsFloatGreaterOrAlmostEqual(double lhs, double rhs, double epsilon)
        {
            return IsFloatAlmostEqual(lhs, rhs, epsilon) || lhs > rhs;
        }

        static constexpr double ClampFloat(double min, double max, double value)
        {
            return std::clamp(min, max, value);
        }
    }
}

#endif