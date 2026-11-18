#ifndef PHARE_CORE_MODELS_EXTERNAL_FIELD_STATE_HPP
#define PHARE_CORE_MODELS_EXTERNAL_FIELD_STATE_HPP

#include "core/def.hpp"
#include "core/utilities/types.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace PHARE::core
{
/**
 * @brief The external magnetic field and its time derivative.
 *
 * B0 is held only when an external field is declared, and dB0dt only when it also depends on
 * time: only those are registered and allocated, B0() and dB0dt() throw otherwise.
 *
 * @tparam VecFieldT vecfield implementation
 *
 */
template<typename VecFieldT>
class ExternalFieldState
{
public:
    using vecfield_type = VecFieldT;
    using tensor_type   = VecFieldT::tensor_t;

    static constexpr std::size_t dimension = VecFieldT::dimension;

    ExternalFieldState(std::string const& name, bool const enabled,
                       bool const timeDependent = false)
    {
        if (timeDependent and not enabled)
            throw std::runtime_error("a time dependent external field must be declared");
        if (enabled)
            fields_.emplace_back(name + "_B0", tensor_type::B);
        if (timeDependent)
            fields_.emplace_back(name + "_dB0dt", tensor_type::B);
    }

    NO_DISCARD bool isEnabled() const { return fields_.size() > 0; }
    NO_DISCARD bool isTimeDependent() const { return fields_.size() > 1; }

    NO_DISCARD VecFieldT& B0() { return field_(0, "no external field was declared"); }
    NO_DISCARD VecFieldT const& B0() const { return field_(0, "no external field was declared"); }

    NO_DISCARD VecFieldT& dB0dt() { return field_(1, "the external field is not time dependent"); }
    NO_DISCARD VecFieldT const& dB0dt() const
    {
        return field_(1, "the external field is not time dependent");
    }

    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const
    {
        return core::all(fields_, [](auto const& field) { return field.isUsable(); });
    }

    NO_DISCARD bool isSettable() const
    {
        return core::all(fields_, [](auto const& field) { return field.isSettable(); });
    }

    NO_DISCARD std::vector<VecFieldT>& getRunTimeResourcesViewList() { return fields_; }
    NO_DISCARD std::vector<VecFieldT> const& getRunTimeResourcesViewList() const { return fields_; }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

private:
    VecFieldT& field_(std::size_t const i, std::string const& absent)
    {
        throwIfAbsent_(i, absent);
        return fields_[i];
    }

    VecFieldT const& field_(std::size_t const i, std::string const& absent) const
    {
        throwIfAbsent_(i, absent);
        return fields_[i];
    }

    void throwIfAbsent_(std::size_t const i, std::string const& absent) const
    {
        if (i >= fields_.size())
            throw std::runtime_error(absent);
    }

    std::vector<VecFieldT> fields_;
};

} // namespace PHARE::core

#endif // PHARE_CORE_MODELS_EXTERNAL_FIELD_STATE_HPP
