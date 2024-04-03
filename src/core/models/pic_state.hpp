#ifndef PHARE_PIC_STATE_HPP
#define PHARE_PIC_STATE_HPP

#include "core/pic/pic_quantities.hpp"
#include "core/models/physical_state.hpp"
#include "initializer/data_provider.hpp"
#include "core/def.hpp"
#include "core/utilities/algorithm.hpp"

#include <cstddef>
#include <sstream>
#include <string>
#include <utility>

namespace PHARE::core
{

/**
 * @brief The PICState class is a concrete implementation of a IPhysicalState.
 * It holds an Electromag and an Fermion object manipulated by a  PIC concrete type of
 * ISolver
 */
template<typename Electromag, typename Fermions>
class PICState : public IPhysicalState
{

using VecField = typename Electromag::vecfield_type;
    
public:

    static constexpr auto dimension = Fermions::dimension;

    PICState(PHARE::initializer::PHAREDict const& dict)
        : electromag{dict["electromag"]}
        , fermions{} // TODO
        , J{"J", PICQuantity::Vector::J}
    {
    }

    Electromag electromag;
    Fermions fermions;
    VecField J;

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

        NO_DISCARD bool isUsable() const
        {
            return electromag.isUsable() && fermions.isUsable() && J.isUsable();
        }



        NO_DISCARD bool isSettable() const
        {
            return electromag.isSettable() && fermions.isSettable() && J.isSettable();
        }


        NO_DISCARD auto getCompileTimeResourcesUserList() const
        {
            return std::forward_as_tuple(electromag, fermions);
        }

        NO_DISCARD auto getCompileTimeResourcesUserList()
        {
            return std::forward_as_tuple(electromag, fermions);
        }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

};// end PICState 

} // namespace PHARE::core


#endif // PHARE_PIC_STATE_HPP
