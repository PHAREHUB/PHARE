

#ifndef PHARE_CORE_DATA_ELECTROMAG_ELECTROMAG_H
#define PHARE_CORE_DATA_ELECTROMAG_ELECTROMAG_H

#include <string>

#include <hybrid/hybrid_quantities.h>

namespace PHARE
{
template<typename VecFieldT>
struct Electromag
{
    explicit Electromag(std::string name)
        : E{name + "_E", HybridQuantity::Vector::E}
        , B{name + "_B", HybridQuantity::Vector::B}
    {
    }



    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    bool isUsable() const { return E.isUsable() && B.isUsable(); }



    bool isSettable() const { return E.isSettable() && B.isSettable(); }


    auto getCompileTimeResourcesUserList() { return std::forward_as_tuple(E, B); }


    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------



    VecFieldT E;
    VecFieldT B;
};
} // namespace PHARE
#endif
