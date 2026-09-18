#ifndef PHARE_CORE_MODELS_EXTERNAL_FIELD_HPP
#define PHARE_CORE_MODELS_EXTERNAL_FIELD_HPP

#include "core/def.hpp"

#include <string>
#include <tuple>

namespace PHARE::core
{
/**
 * @brief The external magnetic field and its time derivative.
 *
 * @tparam VecFieldT vecfield implementation
 *
 */
template<typename VecFieldT>
class ExternalField
{
public:
    using vecfield_type = VecFieldT;
    using tensor_type   = VecFieldT::tensor_t;

    static constexpr std::size_t dimension = VecFieldT::dimension;

    explicit ExternalField(std::string const& name)
        : B0{name + "_B0", tensor_type::B}
        , dB0dt{name + "_dB0dt", tensor_type::B}
    {
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const { return core::isUsable(B0, dB0dt); }

    NO_DISCARD bool isSettable() const { return core::isSettable(B0, dB0dt); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(B0, dB0dt);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(B0, dB0dt); }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

    VecFieldT B0;
    VecFieldT dB0dt;
};

} // namespace PHARE::core

#endif // PHARE_CORE_MODELS_EXTERNAL_FIELD_HPP
