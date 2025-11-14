#ifndef PHARE_CORE_NUMERICS_PRIMITIVE_CONSERVATIVE_CONVERTER_HPP
#define PHARE_CORE_NUMERICS_PRIMITIVE_CONSERVATIVE_CONVERTER_HPP

#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"
#include "initializer/data_provider.hpp"

namespace PHARE::core
{
inline auto vToRhoV(auto const& rho, auto const& Vx, auto const& Vy, auto const& Vz)
{
    auto const rhoVx = rho * Vx;
    auto const rhoVy = rho * Vy;
    auto const rhoVz = rho * Vz;

    return std::make_tuple(rhoVx, rhoVy, rhoVz);
}

inline auto eosPToEtot(double const gamma, auto const& rho, auto const& vx, auto const& vy,
                       auto const& vz, auto const& bx, auto const& by, auto const& bz,
                       auto const& p)
{
    auto const v2 = vx * vx + vy * vy + vz * vz;
    auto const b2 = bx * bx + by * by + bz * bz;
    return p / (gamma - 1.0) + 0.5 * rho * v2 + 0.5 * b2;
}

template<typename GridLayout>
class ToConservativeConverter_ref;

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
    void operator()(Field const& rho, VecField const& V, VecField const& B, Field const& P,
                    VecField& rhoV, Field& Etot) const
    {
        ToConservativeConverter_ref<GridLayout>{*this->layout_, gamma_}(rho, V, B, P, rhoV, Etot);
    }

private:
    double const gamma_;
};

template<typename GridLayout>
class ToConservativeConverter_ref
{
    constexpr static auto dimension = GridLayout::dimension;

public:
    ToConservativeConverter_ref(GridLayout const& layout, double const gamma)
        : layout_{layout}
        , gamma_{gamma}
    {
    }

    template<typename Field, typename VecField>
    void operator()(Field const& rho, VecField const& V, VecField const& B, Field const& P,
                    VecField& rhoV, Field& Etot) const
    {
        layout_.evalOnGhostBox(rho,
                               [&](auto&... args) mutable { vToRhoV_(rho, V, rhoV, {args...}); });

        layout_.evalOnGhostBox(rho, [&](auto&... args) mutable {
            eosPToEtot_(gamma_, rho, V, B, P, Etot, {args...});
        });
    }

private:
    template<typename Field, typename VecField>
    static void vToRhoV_(Field const& rho, VecField const& V, VecField& rhoV,
                         MeshIndex<Field::dimension> index)
    {
        auto const& Vx = V(Component::X);
        auto const& Vy = V(Component::Y);
        auto const& Vz = V(Component::Z);

        auto& rhoVx = rhoV(Component::X);
        auto& rhoVy = rhoV(Component::Y);
        auto& rhoVz = rhoV(Component::Z);

        auto&& [x, y, z] = vToRhoV(rho(index), Vx(index), Vy(index), Vz(index));
        rhoVx(index)     = x;
        rhoVy(index)     = y;
        rhoVz(index)     = z;
    }

    template<typename Field, typename VecField>
    static void eosPToEtot_(double const gamma, Field const& rho, VecField const& V,
                            VecField const& B, Field const& P, Field& Etot,
                            MeshIndex<Field::dimension> index)
    {
        auto const& Vx = V(Component::X);
        auto const& Vy = V(Component::Y);
        auto const& Vz = V(Component::Z);

        auto const& Bx = B(Component::X);
        auto const& By = B(Component::Y);
        auto const& Bz = B(Component::Z);

        auto const bx = GridLayout::project(Bx, index, GridLayout::faceXToCellCenter());
        auto const by = GridLayout::project(By, index, GridLayout::faceYToCellCenter());
        auto const bz = GridLayout::project(Bz, index, GridLayout::faceZToCellCenter());

        Etot(index)
            = eosPToEtot(gamma, rho(index), Vx(index), Vy(index), Vz(index), bx, by, bz, P(index));
    }

private:
    GridLayout layout_;

    double const gamma_;
};

} // namespace PHARE::core

#endif
