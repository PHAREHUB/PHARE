#ifndef PHARE_CORE_NUMERICS_PRIMITIVE_CONSERVATIVE_CONVERTER_HPP
#define PHARE_CORE_NUMERICS_PRIMITIVE_CONSERVATIVE_CONVERTER_HPP

#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"
#include "initializer/data_provider.hpp"
#include <tuple>

namespace PHARE::core
{
template<typename GridLayout>
class ToConservativeConverter : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:
    ToConservativeConverter(PHARE::initializer::PHAREDict const& dict)
        : gamma_{dict["heat_capacity_ratio"].template to<double>()}
    {
    }

    template<typename Field, typename VecField>
    void vToRhoV(const Field& rho, const VecField& V, VecField& rhoV) const
    {
        auto const& Vx = V(Component::X);
        auto const& Vy = V(Component::Y);
        auto const& Vz = V(Component::Z);

        auto& rhoVx = rhoV(Component::X);
        auto& rhoVy = rhoV(Component::Y);
        auto& rhoVz = rhoV(Component::Z);

        auto v_to_rhov = [&](const Field& V, Field& rhoV, MeshIndex<Field::dimension> index) {
            rhoV(index) = rho(index) * V(index);
        };

        layout_->evalOnBox(rhoVx, [&](auto&... args) mutable { v_to_rhov(Vx, rhoVx, {args...}); });
        layout_->evalOnBox(rhoVy, [&](auto&... args) mutable { v_to_rhov(Vy, rhoVy, {args...}); });
        layout_->evalOnBox(rhoVz, [&](auto&... args) mutable { v_to_rhov(Vz, rhoVz, {args...}); });
    }

    template<typename Field, typename VecField>
    void eosPToEtot(const Field& rho, const VecField& V, const VecField& B, const Field& P,
                    Field& Etot) const
    {
        auto const& Vx = V(Component::X);
        auto const& Vy = V(Component::Y);
        auto const& Vz = V(Component::Z);

        auto const& Bx = B(Component::X);
        auto const& By = B(Component::Y);
        auto const& Bz = B(Component::Z);

        auto eos_p_to_etot = [&](MeshIndex<Field::dimension> index) {
            auto bx     = GridLayout::project(Bx, index, GridLayout::faceXToCellCenter());
            auto by     = GridLayout::project(By, index, GridLayout::faceYToCellCenter());
            auto bz     = GridLayout::project(Bz, index, GridLayout::faceZToCellCenter());
            auto v2     = Vx(index) * Vx(index) + Vy(index) * Vy(index) + Vz(index) * Vz(index);
            auto b2     = bx * bx + by * by + bz * bz;
            Etot(index) = P(index) / (gamma_ - 1.0) + 0.5 * rho(index) * v2 + 0.5 * b2;
        };

        layout_->evalOnBox(Etot, [&](auto&... args) mutable { eos_p_to_etot({args...}); });
    }

private:
    double const gamma_;
};

} // namespace PHARE::core

#endif
