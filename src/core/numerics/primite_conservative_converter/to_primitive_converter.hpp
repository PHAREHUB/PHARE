#ifndef PHARE_CORE_NUMERICS_TO_PRIMITIVE_CONVERTER_HPP
#define PHARE_CORE_NUMERICS_TO_PRIMITIVE_CONVERTER_HPP

#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"
#include "initializer/data_provider.hpp"

namespace PHARE::core
{
auto rhoVToV(auto const& rho, auto const& rhoVx, auto const& rhoVy, auto const& rhoVz)
{
    auto vx = rhoVx / rho;
    auto vy = rhoVy / rho;
    auto vz = rhoVz / rho;

    return std::make_tuple(vx, vy, vz);
}

auto eosEtotToP(double const gamma, auto const& rho, auto const& vx, auto const& vy, auto const& vz,
                auto const& bx, auto const& by, auto const& bz, auto const& etot)
{
    auto v2 = vx * vx + vy * vy + vz * vz;
    auto b2 = bx * bx + by * by + bz * bz;

    return (gamma - 1.0) * (etot - 0.5 * rho * v2 - 0.5 * b2);
}

template<typename GridLayout>
class ToPrimitiveConverter_ref;

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
    void operator()(const Field& rho, const VecField& rhoV, const VecField& B, const Field& Etot,
                    VecField& V, Field& P) const
    {
        ToPrimitiveConverter_ref<GridLayout>{*this->layout_, gamma_}(rho, rhoV, B, Etot, V, P);
    }

private:
    double const gamma_;
};

template<typename GridLayout>
class ToPrimitiveConverter_ref
{
    constexpr static auto dimension = GridLayout::dimension;

public:
    ToPrimitiveConverter_ref(GridLayout const& layout, double const gamma)
        : layout_{layout}
        , gamma_{gamma}
    {
    }

    template<typename Field, typename VecField>
    void operator()(const Field& rho, const VecField& rhoV, const VecField& B, const Field& Etot,
                    VecField& V, Field& P) const
    {
        layout_.evalOnBox(rho, [&](auto&... args) mutable { rhoVToV_(rho, rhoV, V, {args...}); });

        layout_.evalOnBox(rho, [&](auto&... args) mutable {
            eosEtotToP_(gamma_, rho, rhoV, B, Etot, P, {args...});
        });
    }

private:
    template<typename Field, typename VecField>
    static void rhoVToV_(const Field& rho, const VecField& rhoV, VecField& V,
                         MeshIndex<Field::dimension> index)
    {
        auto const& rhoVx = rhoV(Component::X);
        auto const& rhoVy = rhoV(Component::Y);
        auto const& rhoVz = rhoV(Component::Z);

        auto& Vx = V(Component::X);
        auto& Vy = V(Component::Y);
        auto& Vz = V(Component::Z);

        auto&& [x, y, z] = rhoVToV(rho(index), rhoVx(index), rhoVy(index), rhoVz(index));
        Vx(index)        = x;
        Vy(index)        = y;
        Vz(index)        = z;
    }

    template<typename Field, typename VecField>
    static void eosEtotToP_(double const gamma, const Field& rho, const VecField& rhoV,
                            const VecField& B, const Field& Etot, Field& P,
                            MeshIndex<Field::dimension> index)
    {
        auto const& rhoVx = rhoV(Component::X);
        auto const& rhoVy = rhoV(Component::Y);
        auto const& rhoVz = rhoV(Component::Z);

        auto const& Bx = B(Component::X);
        auto const& By = B(Component::Y);
        auto const& Bz = B(Component::Z);

        auto vx  = rhoVx(index) / rho(index);
        auto vy  = rhoVy(index) / rho(index);
        auto vz  = rhoVz(index) / rho(index);
        auto bx  = GridLayout::project(Bx, index, GridLayout::faceXToCellCenter());
        auto by  = GridLayout::project(By, index, GridLayout::faceYToCellCenter());
        auto bz  = GridLayout::project(Bz, index, GridLayout::faceZToCellCenter());
        P(index) = eosEtotToP(gamma, rho(index), vx, vy, vz, bx, by, bz, Etot(index));
    }


private:
    GridLayout layout_;

    double const gamma_;
};

} // namespace PHARE::core

#endif
