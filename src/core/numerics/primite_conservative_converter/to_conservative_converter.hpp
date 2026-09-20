#ifndef PHARE_CORE_NUMERICS_PRIMITIVE_CONSERVATIVE_CONVERTER_HPP
#define PHARE_CORE_NUMERICS_PRIMITIVE_CONSERVATIVE_CONVERTER_HPP


#include "core/data/field/field_box.hpp"
#include "core/utilities/index/index.hpp"
#include "core/utilities/algorithm.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

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
class ToConservativeConverter
{
    constexpr static auto dimension = GridLayout::dimension;

public:
    ToConservativeConverter(GridLayout const& layout, double const gamma)
        : layout_{layout}
        , gamma_{gamma}
    {
    }

    template<typename Field, typename VecField>
    void operator()(Field const& rho, VecField const& V, VecField const& B, Field const& P,
                    VecField& rhoV, Field& Etot) const
    {
        vToRhoV_(rho, V, rhoV);

        eosPToEtotOnGhostBox_(layout_, gamma_, rho, V(Component::X), V(Component::Y),
                              V(Component::Z), B(Component::X), B(Component::Y), B(Component::Z), P,
                              Etot);
    }

private:
    template<typename Field, typename VecField>
    static void vToRhoV_(Field const& rho, VecField const& V, VecField& rhoV)
    {
        // no neighbor access here, so unlike eosPToEtot_ below this can be a flat per-cell loop
        core::operate_on_fields(
            [](auto& rhoVx, auto& rhoVy, auto& rhoVz, auto const& r, auto const& vx, auto const& vy,
               auto const& vz) {
                auto&& [x, y, z] = vToRhoV(r, vx, vy, vz);
                rhoVx            = x;
                rhoVy            = y;
                rhoVz            = z;
            },
            rhoV(Component::X), rhoV(Component::Y), rhoV(Component::Z), rho, V(Component::X),
            V(Component::Y), V(Component::Z));
    }

    template<typename Field>
    static void eosPToEtot_(double const gamma, Field const& rho, Field const& Vx, Field const& Vy,
                            Field const& Vz, Field const& Bx, Field const& By, Field const& Bz,
                            Field const& P, Field& Etot, MeshIndex<Field::dimension> index)
    {
        auto const bx
            = GridLayout::template project<GridLayout::implT::faceXToCellCenter>(Bx, index);
        auto const by
            = GridLayout::template project<GridLayout::implT::faceYToCellCenter>(By, index);
        auto const bz
            = GridLayout::template project<GridLayout::implT::faceZToCellCenter>(Bz, index);

        Etot(index)
            = eosPToEtot(gamma, rho(index), Vx(index), Vy(index), Vz(index), bx, by, bz, P(index));
    }

    // project<faceXToCellCenter> reads Bx's own neighboring cells (see GridLayout::project),
    // so unlike vToRhoV_ this genuinely needs a per-tile GridLayout/ghost-box, not a flat loop.
    template<typename Field>
    static void eosPToEtotOnGhostBox_(GridLayout const& layout, double const gamma,
                                      Field const& rho, Field const& Vx, Field const& Vy,
                                      Field const& Vz, Field const& Bx, Field const& By,
                                      Field const& Bz, Field const& P, Field& Etot)
        requires(not is_field_tile_set_v<Field>)
    {
        layout.evalOnGhostBox(rho, [&](auto&... args) mutable {
            eosPToEtot_(gamma, rho, Vx, Vy, Vz, Bx, By, Bz, P, Etot, {args...});
        });
    }

    template<typename Field>
    static void eosPToEtotOnGhostBox_(GridLayout const&, double const gamma, Field const& rho,
                                      Field const& Vx, Field const& Vy, Field const& Vz,
                                      Field const& Bx, Field const& By, Field const& Bz,
                                      Field const& P, Field& Etot)
        requires(is_field_tile_set_v<Field>)
    {
        auto& rho_tiles  = rho();
        auto& Vx_tiles   = Vx();
        auto& Vy_tiles   = Vy();
        auto& Vz_tiles   = Vz();
        auto& Bx_tiles   = Bx();
        auto& By_tiles   = By();
        auto& Bz_tiles   = Bz();
        auto& P_tiles    = P();
        auto& Etot_tiles = Etot();

        // rho, V, B, P, Etot are all defined on the same patch, so they share the same tiling
        for (std::size_t t = 0; t < rho_tiles.size(); ++t)
            eosPToEtotOnGhostBox_(rho_tiles[t].layout(), gamma, rho_tiles[t](), Vx_tiles[t](),
                                  Vy_tiles[t](), Vz_tiles[t](), Bx_tiles[t](), By_tiles[t](),
                                  Bz_tiles[t](), P_tiles[t](), Etot_tiles[t]());
    }

private:
    GridLayout layout_;

    double const gamma_;
};

} // namespace PHARE::core

#endif
