

#ifndef PHARE_CORE_DATA_ELECTROMAG_ELECTROMAG_H
#define PHARE_CORE_DATA_ELECTROMAG_ELECTROMAG_H

#include <string>
#include <tuple>

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

    using vecfield_type = VecFieldT;


    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    bool isUsable() const { return E.isUsable() && B.isUsable(); }



    bool isSettable() const { return E.isSettable() && B.isSettable(); }


    auto getCompileTimeResourcesUserList() const { return std::forward_as_tuple(E, B); }

    auto getCompileTimeResourcesUserList() { return std::forward_as_tuple(E, B); }


    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------



    VecFieldT E;
    VecFieldT B;

    static constexpr std::size_t dimension = VecFieldT::dimension;
};
} // namespace PHARE
#endif
