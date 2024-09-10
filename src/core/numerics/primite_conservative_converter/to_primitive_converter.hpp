#ifndef PHARE_CORE_NUMERICS_TO_PRIMITIVE_CONVERTER_HPP
#define PHARE_CORE_NUMERICS_TO_PRIMITIVE_CONVERTER_HPP

#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"
#include "initializer/data_provider.hpp"
#include <tuple>

namespace PHARE::core
{
template<typename GridLayout>
class ToPrimitiveConverter : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:
    ToPrimitiveConverter(PHARE::initializer::PHAREDict const& dict)
        : gamma_{dict["heat_capacity_ratio"].template to<double>()}
    {
    }

    template<typename Field, typename VecField>
    void rhoVToV(const Field& rho, const VecField& rhoV, VecField& V) const
    {
        auto const& rhoVx = rhoV(Component::X);
        auto const& rhoVy = rhoV(Component::Y);
        auto const& rhoVz = rhoV(Component::Z);

        auto& Vx = V(Component::X);
        auto& Vy = V(Component::Y);
        auto& Vz = V(Component::Z);

        auto rhov_to_v = [&](const Field& rhoV, Field& V, MeshIndex<Field::dimension> index) {
            V(index) = rhoV(index) / rho(index);
        };

        layout_->evalOnBox(Vx, [&](auto&... args) mutable { rhov_to_v(rhoVx, Vx, {args...}); });
        layout_->evalOnBox(Vy, [&](auto&... args) mutable { rhov_to_v(rhoVy, Vy, {args...}); });
        layout_->evalOnBox(Vz, [&](auto&... args) mutable { rhov_to_v(rhoVz, Vz, {args...}); });
    }

    template<typename Field, typename VecField>
    void eosEtotToP(const Field& rho, const VecField& rhoV, const VecField& B, const Field& Etot,
                    Field& P) const
    {
        auto const& rhoVx = rhoV(Component::X);
        auto const& rhoVy = rhoV(Component::Y);
        auto const& rhoVz = rhoV(Component::Z);

        auto const& Bx = B(Component::X);
        auto const& By = B(Component::Y);
        auto const& Bz = B(Component::Z);

        auto eos_etot_to_p = [&](MeshIndex<Field::dimension> index) {
            auto vx  = rhoVx(index) / rho(index);
            auto vy  = rhoVy(index) / rho(index);
            auto vz  = rhoVz(index) / rho(index);
            auto bx  = GridLayout::project(Bx, index, GridLayout::faceXToCellCenter());
            auto by  = GridLayout::project(By, index, GridLayout::faceYToCellCenter());
            auto bz  = GridLayout::project(Bz, index, GridLayout::faceZToCellCenter());
            auto v2  = vx * vx + vy * vy + vz * vz;
            auto b2  = bx * bx + by * by + bz * bz;
            P(index) = (gamma_ - 1.0) * (Etot(index) - 0.5 * rho(index) * v2 - 0.5 * b2);
        };

        layout_->evalOnBox(P, [&](auto&... args) mutable { eos_etot_to_p({args...}); });
    }

private:
    double const gamma_;
};

} // namespace PHARE::core

#endif
