#ifndef PHARE_CORE_NUMERICS_TO_PRIMITIVE_CONVERTER_HPP
#define PHARE_CORE_NUMERICS_TO_PRIMITIVE_CONVERTER_HPP


#include "core/utilities/index/index.hpp"
#include "core/utilities/algorithm.hpp"
#include "core/data/vecfield/vecfield_component.hpp"

namespace PHARE::core
{
static auto const min_value = std::sqrt(1024 * std::numeric_limits<double>::min());

auto rhoVToV(auto& rho, auto const& rhoVx, auto const& rhoVy, auto const& rhoVz)
{
    auto const vx = rhoVx / rho;
    auto const vy = rhoVy / rho;
    auto const vz = rhoVz / rho;

    return std::make_tuple(vx, vy, vz);
}

auto eosEtotToP(double const gamma, auto const& rho, auto const& vx, auto const& vy, auto const& vz,
                auto const& bx, auto const& by, auto const& bz, auto& etot)
{
    auto const v2 = vx * vx + vy * vy + vz * vz;
    auto const b2 = bx * bx + by * by + bz * bz;

    auto p = (gamma - 1.0) * (etot - 0.5 * rho * v2 - 0.5 * b2);
    // p      = (p < 0.) ? 1.0e-5 : p; //tbd maybe not needed
    // etot = p / (gamma - 1.0) + 0.5 * rho * v2 + 0.5 * b2;

    return p;
}



template<typename GridLayout>
class ToPrimitiveConverter
{
    constexpr static auto dimension = GridLayout::dimension;

public:
    ToPrimitiveConverter(GridLayout const& layout)
        : layout_{layout}
    {
    }

    template<typename Field, typename VecField>
    void operator()(double const gamma, Field& rho, VecField const& rhoV, VecField const& B,
                    Field& Etot, VecField& V, Field& P) const
    {
        rhoVToVOnGhostBox(rho, rhoV, V);

        eosEtotToPOnGhostBox(gamma, rho, rhoV, B, Etot, P);
    }

    // used for diagnostics
    template<typename Field, typename VecField>
    void rhoVToVOnGhostBox(Field& rho, VecField const& rhoV, VecField& V) const
    {
        rhoVToV_(rho, rhoV, V);
    }

    template<typename Field, typename VecField>
    void eosEtotToPOnGhostBox(double const gamma, Field const& rho, VecField const& rhoV,
                              VecField const& B, Field& Etot, Field& P) const
    {
        eosEtotToPOnGhostBox_(layout_, gamma, rho, rhoV(Component::X), rhoV(Component::Y),
                              rhoV(Component::Z), B(Component::X), B(Component::Y), B(Component::Z),
                              Etot, P);
    }

private:
    template<typename Field, typename VecField>
    static void rhoVToV_(Field& rho, VecField const& rhoV, VecField& V)
    {
        // no neighbor access here, so unlike eosEtotToP_ below this can be a flat per-cell loop
        core::operate_on_fields(
            [](auto& vx, auto& vy, auto& vz, auto const& r, auto const& rhovx, auto const& rhovy,
               auto const& rhovz) {
                auto&& [x, y, z] = rhoVToV(r, rhovx, rhovy, rhovz);
                vx               = x;
                vy               = y;
                vz               = z;
            },
            V(Component::X), V(Component::Y), V(Component::Z), rho, rhoV(Component::X),
            rhoV(Component::Y), rhoV(Component::Z));
    }

    template<typename Field>
    static void eosEtotToP_(double const gamma, Field const& rho, Field const& rhoVx,
                            Field const& rhoVy, Field const& rhoVz, Field const& Bx,
                            Field const& By, Field const& Bz, Field& Etot, Field& P,
                            MeshIndex<Field::dimension> index)
    {
        auto const vx = rhoVx(index) / rho(index);
        auto const vy = rhoVy(index) / rho(index);
        auto const vz = rhoVz(index) / rho(index);
        auto const bx
            = GridLayout::template project<GridLayout::implT::faceXToCellCenter>(Bx, index);
        auto const by
            = GridLayout::template project<GridLayout::implT::faceYToCellCenter>(By, index);
        auto const bz
            = GridLayout::template project<GridLayout::implT::faceZToCellCenter>(Bz, index);
        P(index) = eosEtotToP(gamma, rho(index), vx, vy, vz, bx, by, bz, Etot(index));
    }

    // project<faceXToCellCenter> reads Bx's own neighboring cells (see GridLayout::project),
    // so unlike rhoVToV_ this genuinely needs a per-tile GridLayout/ghost-box, not a flat loop.
    template<typename Field>
    static void eosEtotToPOnGhostBox_(GridLayout const& layout, double const gamma,
                                      Field const& rho, Field const& rhoVx, Field const& rhoVy,
                                      Field const& rhoVz, Field const& Bx, Field const& By,
                                      Field const& Bz, Field& Etot, Field& P)
        requires(not is_field_tile_set_v<Field>)
    {
        layout.evalOnGhostBox(rho, [&](auto&... args) mutable {
            eosEtotToP_(gamma, rho, rhoVx, rhoVy, rhoVz, Bx, By, Bz, Etot, P, {args...});
        });
    }

    template<typename Field>
    static void eosEtotToPOnGhostBox_(GridLayout const&, double const gamma, Field const& rho,
                                      Field const& rhoVx, Field const& rhoVy, Field const& rhoVz,
                                      Field const& Bx, Field const& By, Field const& Bz,
                                      Field& Etot, Field& P)
        requires(is_field_tile_set_v<Field>)
    {
        auto& rho_tiles   = rho();
        auto& rhoVx_tiles = rhoVx();
        auto& rhoVy_tiles = rhoVy();
        auto& rhoVz_tiles = rhoVz();
        auto& Bx_tiles    = Bx();
        auto& By_tiles    = By();
        auto& Bz_tiles    = Bz();
        auto& Etot_tiles  = Etot();
        auto& P_tiles     = P();

        // rho, rhoV, B, Etot, P are all defined on the same patch, so they share the same tiling
        for (std::size_t t = 0; t < rho_tiles.size(); ++t)
            eosEtotToPOnGhostBox_(rho_tiles[t].layout(), gamma, rho_tiles[t](), rhoVx_tiles[t](),
                                  rhoVy_tiles[t](), rhoVz_tiles[t](), Bx_tiles[t](), By_tiles[t](),
                                  Bz_tiles[t](), Etot_tiles[t](), P_tiles[t]());
    }


private:
    GridLayout layout_;
};

} // namespace PHARE::core

#endif
