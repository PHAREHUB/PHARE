#ifndef PHARE_CORE_MODELS_EXTERNAL_FIELD_UPDATER_HPP
#define PHARE_CORE_MODELS_EXTERNAL_FIELD_UPDATER_HPP

#include "core/def.hpp"
#include "core/models/external_field.hpp"
#include "core/numerics/differential_operators.hpp"
#include "core/utilities/point/point.hpp"

#include <cstddef>

namespace PHARE::core
{


/**
 * @brief Interface of the external field updaters.
 *
 * @tparam VecFieldT vecfield implementation
 * @tparam GridLayoutT grid layout implementation
 *
 * Implementations provide the vector potential @f$\mathbf{A}_0@f$ of the external field, and
 * its time derivative when the field is time dependent. The external field itself is obtained
 * here as @f$\mathbf{B}_0 = \nabla\times\mathbf{A}_0@f$, so that it is divergence free by
 * construction at the discrete level.
 */
template<typename VecFieldT, typename GridLayoutT>
class IExternalFieldUpdater
{
public:
    using vecfield_type       = VecFieldT;
    using field_type          = VecFieldT::field_type;
    using value_type          = field_type::value_type;
    using tensor_type         = VecFieldT::tensor_t;
    using component_type      = VecFieldT::component_type;
    using point_type          = Point<double, GridLayoutT::dimension>;
    using external_field_type = ExternalField<VecFieldT>;

    static constexpr std::size_t dimension = GridLayoutT::dimension;

    virtual ~IExternalFieldUpdater() = default;

    /**
     * @brief fill an external field at a given time
     *
     * @param externalField the external field to fill
     * @param scratch E-centered work vecfield, holding the vector potential on output
     * @param layout the current grid layout
     * @param time the current time
     */
    virtual void operator()(external_field_type& externalField, vecfield_type& scratch,
                            GridLayoutT const& layout, double time)
    {
        computePotential(scratch, time, layout);
        curlOnGhostBox_(externalField.B0, scratch, layout);
        if (isTimeDependent())
        {
            computePotentialTimeDerivative(scratch, time, layout);
            curlOnGhostBox_(externalField.dB0dt, scratch, layout);
        }
        else
        {
            externalField.dB0dt.zero();
        }
    }

    virtual void computePotential(vecfield_type& a0, double time, GridLayoutT const& layout) = 0;
    virtual void computePotentialTimeDerivative(vecfield_type& da0_dt, double time,
                                                GridLayoutT const& layout)
        = 0;

    NO_DISCARD virtual bool isTimeDependent() const = 0;

private:
    void curlOnGhostBox_(vecfield_type& out, vecfield_type const& in, GridLayoutT const& layout)
    {
        auto& outX = out(component_type::X);
        auto& outY = out(component_type::Y);
        auto& outZ = out(component_type::Z);

        auto const curlIn = curl(layout, in);

        layout.evalOnGhostBox(outX, [&](auto const&... ijk) { outX(ijk...) = curlIn.x({ijk...}); });
        layout.evalOnGhostBox(outY, [&](auto const&... ijk) { outY(ijk...) = curlIn.y({ijk...}); });
        layout.evalOnGhostBox(outZ, [&](auto const&... ijk) { outZ(ijk...) = curlIn.z({ijk...}); });
    }
};

} // namespace PHARE::core

#endif // PHARE_CORE_MODELS_EXTERNAL_FIELD_UPDATER_HPP
