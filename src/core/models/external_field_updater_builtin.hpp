#ifndef PHARE_CORE_MODELS_EXTERNAL_FIELD_UPDATER_BUILTIN_HPP
#define PHARE_CORE_MODELS_EXTERNAL_FIELD_UPDATER_BUILTIN_HPP

#include "core/models/external_field_updater.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/types.hpp"

namespace PHARE::core
{
/**
 * @brief CRTP base class for external fields directly implemented in PHARE.
 *
 * @tparam DerivedT CRTP derived class
 * @tparam VecFieldT vecfield implementation
 * @tparam GridLayoutT grid layout implementation
 *
 * To define a new external field, create a class that inherits this one and that implements a
 * point-wise formula for the vector potential, and its time derivative in case it is
 * time-dependent, as:
 *
 * @code
 * template<component_type i>
 * double potential(point_type const& coords, double time) const;
 *
 * template<component_type i> // only if time dependent
 * double potentialTimeDerivative(point_type const& coords, double time) const;
 * @endcode
 *
 */
template<typename DerivedT, typename VecFieldT, typename GridLayoutT>
class ExternalFieldUpdaterBuiltin : public IExternalFieldUpdater<VecFieldT, GridLayoutT>
{
public:
    using Super          = IExternalFieldUpdater<VecFieldT, GridLayoutT>;
    using vecfield_type  = Super::vecfield_type;
    using value_type     = Super::value_type;
    using point_type     = Super::point_type;
    using component_type = Super::component_type;

    static constexpr std::size_t dimension = Super::dimension;

    NO_DISCARD bool isTimeDependent() const final { return isTimeDependent_(); }

    void computePotential(vecfield_type& a0, double time, GridLayoutT const& layout) final
    {
        applyPointWiseFormula_(a0, time, layout, [this](auto c, auto const& coords, double t) {
            return derived().template potential<decltype(c)::value>(coords, t);
        });
    }

    void computePotentialTimeDerivative(vecfield_type& da0_dt, double time,
                                        GridLayoutT const& layout) final
    {
        if constexpr (isTimeDependent_())
        {
            applyPointWiseFormula_(
                da0_dt, time, layout, [this](auto c, auto const& coords, double t) {
                    return derived().template potentialTimeDerivative<decltype(c)::value>(coords,
                                                                                          t);
                });
        }
    }

private:
    /// downcast CRTP helpers
    DerivedT& derived() { return static_cast<DerivedT&>(*this); }
    DerivedT const& derived() const { return static_cast<DerivedT const&>(*this); }

    /**
     * @brief determine automatically if the external field is time dependent by checking if the
     * implementation defines a time derivative
     */
    static constexpr bool isTimeDependent_()
    {
        return requires(DerivedT const& d, point_type const& coords, double time) {
            {
                d.template potentialTimeDerivative<component_type::X>(coords, time)
            } -> std::convertible_to<double>;
        };
    }

    /**
     * @brief helper to apply a point-wise formula to a vector field
     *
     * @tparam Fn a callable invocable as fn(std::integral_constant<component_type, C>{}, point_type
     * const&, double time), for each C in {X, Y, Z}.
     */
    template<typename Fn>
    void applyPointWiseFormula_(vecfield_type& a0,         //!< vecfield to fill
                                double time,               //!< current time
                                GridLayoutT const& layout, //!< current grid layout
                                Fn&& f                     //!< the callable to apply
    ) const
    {
        auto fields = a0.components();
        for_N<3>([&](auto i) {
            constexpr auto component = static_cast<component_type>(decltype(i)::value);
            auto& field              = std::get<i>(fields);
            layout.evalOnGhostBox(field, [&](auto... ijk) {
                auto const coords
                    = layout.fieldNodeCoordinates(field, layout.localToAMR(Point{ijk...}));
                field(ijk...)
                    = f(std::integral_constant<component_type, component>{}, coords, time);
            });
        });
    }
};

} // namespace PHARE::core

#endif // PHARE_CORE_MODELS_EXTERNAL_FIELD_UPDATER_BUILTIN_HPP
