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
 * all particle species (ions and electrons alike, in the Ions container), the particle
 * current J, and the charge density rho of the particles half a step back, the one Gauss's law
 * pairs with the E before the Ampere update (divergence cleaning only, primal).
 * There are no fluid electrons.
 */
template<typename Electromag, typename Ions>
class FullyKineticState : public IPhysicalState
{
    using VecField = Electromag::vecfield_type;
    using Field    = VecField::field_type;

public:
    static constexpr auto dimension = Ions::dimension;

    FullyKineticState(PHARE::initializer::PHAREDict const& dict)
        : electromag{dict["electromag"]}
        , ions{dict["ions"]}
        , J{"J", HybridQuantity::Vector::J}
        , rho{"FK_rho", HybridQuantity::Scalar::rho}
    {
    }

    Electromag electromag;
    Ions ions;
    VecField J;
    Field rho;

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
        return electromag.isUsable() and ions.isUsable() and J.isUsable() and rho.isUsable();
    }

    NO_DISCARD bool isSettable() const
    {
        return electromag.isSettable() and ions.isSettable() and J.isSettable()
               and rho.isSettable();
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(electromag, ions, J, rho);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(electromag, ions, J, rho);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------
};

} // namespace PHARE::core


#endif // PHARE_CORE_MODELS_FULLY_KINETIC_STATE_HPP
