#ifndef PHARE_CORE_NUMERICS_CONSTRAINED_TRANSPORT_DISSIPATIVE_ELECTRIC_FIELD_HPP
#define PHARE_CORE_NUMERICS_CONSTRAINED_TRANSPORT_DISSIPATIVE_ELECTRIC_FIELD_HPP

#include "core/def.hpp"
#include "core/numerics/ohm/ohm.hpp"
#include "core/utilities/index/index.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "core/models/quantities/mhd_quantities.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

#include <cmath>
#include <vector>
#include <algorithm>

namespace PHARE::core
{

template<typename VecField>
class DissipativeElectricFieldState
{
public:
    DissipativeElectricFieldState() = default;
    explicit DissipativeElectricFieldState(bool const isDissipative)
    {
        if (isDissipative)
            E_.emplace_back("E_diss", MHDQuantity::Vector::E);
    }

    NO_DISCARD std::vector<VecField>& getRunTimeResourcesViewList() { return E_; }
    NO_DISCARD std::vector<VecField> const& getRunTimeResourcesViewList() const { return E_; }

    NO_DISCARD auto& E() { return E_[0]; }
    NO_DISCARD auto const& E() const { return E_[0]; }

private:
    std::vector<VecField> E_;
};


template<typename GridLayout>
class DissipativeElectricField : public OhmInfo
{
    using Super                     = OhmInfo;
    constexpr static auto dimension = GridLayout::dimension;

public:
    using Info_t = Super;

    DissipativeElectricField(OhmInfo const& info, GridLayout const& layout)
        : Super{info}
        , layout_{layout}
    {
    }

    void operator()(auto& dissipative_electric_state, auto const& mhd_state) const
    {
        auto& E         = dissipative_electric_state.E();
        auto const& J   = mhd_state.J;
        auto const& B   = mhd_state.B;
        auto const& rho = mhd_state.rho;

        Constexprifier{isResistive(), isHyperResistive(),
                       hyper_mode}([&]<bool isResistive, bool isHyperResistive, HyperMode hyper>() {
            for_N<3>([&](auto i) {
                constexpr auto component = static_cast<Component>(i());
                auto& Ec                 = E(component);
                auto const& Jc           = J(component);
                layout_.evalOnBox(Ec, [&](auto&... args) {
                    MeshIndex<dimension> idx{args...};
                    double e = 0.;
                    if constexpr (isResistive)
                        e += eta * Jc(idx);
                    if constexpr (isHyperResistive)
                        e -= hyper_coef_<component, hyper>(B, rho, idx)
                             * layout_.laplacian(Jc, idx);
                    Ec(idx) = e;
                });
            });
        });
    }

private:
    template<auto component, HyperMode hyper>
    double hyper_coef_(auto const& B, auto const& rho, MeshIndex<dimension> idx) const
    {
        if constexpr (hyper == HyperMode::constant)
            return nu;
        else
        {
            auto const& meshSize = layout_.meshSize();
            auto const dx        = *std::min_element(meshSize.begin(), meshSize.end());

            auto coef = [&]<auto BxProj, auto ByProj, auto BzProj, auto rhoProj>() {
                auto const bx = GridLayout::template project<BxProj>(B(Component::X), idx);
                auto const by = GridLayout::template project<ByProj>(B(Component::Y), idx);
                auto const bz = GridLayout::template project<BzProj>(B(Component::Z), idx);
                auto const n  = GridLayout::template project<rhoProj>(rho, idx);
                auto const b  = std::sqrt(bx * bx + by * by + bz * bz);
                return nu * dx * dx * (b / n + 1);
            };

            if constexpr (component == Component::X)
                return coef
                    .template operator()<GridLayout::BxToEx, GridLayout::ByToEx, GridLayout::BzToEx,
                                         GridLayout::implT::cellCenterToEdgeX>();
            else if constexpr (component == Component::Y)
                return coef
                    .template operator()<GridLayout::BxToEy, GridLayout::ByToEy, GridLayout::BzToEy,
                                         GridLayout::implT::cellCenterToEdgeY>();
            else
                return coef
                    .template operator()<GridLayout::BxToEz, GridLayout::ByToEz, GridLayout::BzToEz,
                                         GridLayout::implT::cellCenterToEdgeZ>();
        }
    }

    GridLayout layout_;
};

} // namespace PHARE::core

#endif
