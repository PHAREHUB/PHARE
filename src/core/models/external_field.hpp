#ifndef PHARE_CORE_MODELS_EXTERNAL_FIELD_HPP
#define PHARE_CORE_MODELS_EXTERNAL_FIELD_HPP

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
 * Holds no vecfield, and so registers and allocates nothing, when no external field is
 * declared: B0() and dB0dt() then throw.
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

    ExternalField(std::string const& name, bool const enabled)
    {
        if (enabled)
        {
            fields_.emplace_back(name + "_B0", tensor_type::B);
            fields_.emplace_back(name + "_dB0dt", tensor_type::B);
        }
    }

    NO_DISCARD bool isEnabled() const { return !fields_.empty(); }

    NO_DISCARD VecFieldT& B0() { return field_(0); }
    NO_DISCARD VecFieldT const& B0() const { return field_(0); }

    NO_DISCARD VecFieldT& dB0dt() { return field_(1); }
    NO_DISCARD VecFieldT const& dB0dt() const { return field_(1); }

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
    VecFieldT& field_(std::size_t const i)
    {
        throwIfDisabled_();
        return fields_[i];
    }

    VecFieldT const& field_(std::size_t const i) const
    {
        throwIfDisabled_();
        return fields_[i];
    }

    void throwIfDisabled_() const
    {
        if (!isEnabled())
            throw std::runtime_error("no external field was declared");
    }

    std::vector<VecFieldT> fields_;
};

} // namespace PHARE::core

#endif // PHARE_CORE_MODELS_EXTERNAL_FIELD_HPP
