#ifndef PHARE_ALGORITHM_H
#define PHARE_ALGORITHM_H

#include "utilities/types.h"

namespace PHARE
{
template<uint32 lhs, uint32 rhs>
constexpr uint32 max()
{
    if constexpr (lhs < rhs)
    {
        return rhs;
    }
    else if constexpr (lhs >= rhs)
    {
        return lhs;
    }
}

} // namespace PHARE

#endif
