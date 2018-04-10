

#ifndef PHARE_CORE_DATA_ELECTROMAG_ELECTROMAG_H
#define PHARE_CORE_DATA_ELECTROMAG_ELECTROMAG_H

#include <string>

#include <hybrid/hybrid_quantities.h>

namespace PHARE
{
template<typename VecFieldT>
struct Electromag
{
    Electromag(std::string name)
        : E{name + "_E", HybridQuantity::Quantity::E}
        , B{name + "_B", HybridQuantity::Quantity::B}
    {
    }
    VecFieldT E;
    VecFieldT B;
};
} // namespace PHARE
#endif
