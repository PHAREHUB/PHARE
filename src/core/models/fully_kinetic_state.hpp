#ifndef PHARE_CORE_MODELS_FULLY_KINETIC_STATE_HPP
#define PHARE_CORE_MODELS_FULLY_KINETIC_STATE_HPP


#include "core/def.hpp"
#include "core/utilities/algorithm.hpp"
#include "core/models/physical_state.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"

#include "initializer/data_provider.hpp"


#include <string>
#include <sstream>


namespace PHARE::core
{
/**
 * @brief FullyKineticState holds the state of a fully kinetic model: the electromagnetic field,
 * all particle species (ions and electrons alike, in the Ions container) and the particle
 * current J.
 * There are no fluid electrons.
 */
template<typename Electromag, typename Ions>
class FullyKineticState : public IPhysicalState
{
    using VecField = Electromag::vecfield_type;

public:
    static constexpr auto dimension = Ions::dimension;

    FullyKineticState(PHARE::initializer::PHAREDict const& dict)
        : electromag{dict["electromag"]}
        , ions{dict["ions"]}
        , J{"J", HybridQuantity::Vector::J}
    {
    }

    Electromag electromag;
    Ions ions;
    VecField J;

    NO_DISCARD std::string to_str()
    {
        std::stringstream ss;
        ss << "Fully Kinetic State\n";
        ss << "------------------------------------\n";
        ss << core::to_str(ions);
        return ss.str();
    }


    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const
    {
        return electromag.isUsable() and ions.isUsable() and J.isUsable();
    }

    NO_DISCARD bool isSettable() const
    {
        return electromag.isSettable() and ions.isSettable() and J.isSettable();
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(electromag, ions, J);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(electromag, ions, J);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------
};

} // namespace PHARE::core


#endif // PHARE_CORE_MODELS_FULLY_KINETIC_STATE_HPP
