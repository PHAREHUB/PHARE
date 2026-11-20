#ifndef PHARE_CORE_NUMERICS_FINITE_VOLUME_EULER_HPP
#define PHARE_CORE_NUMERICS_FINITE_VOLUME_EULER_HPP


#include "core/def.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/index/index.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"

#include <array>
#include <tuple>


namespace PHARE::core
{


template<typename GridLayout>
class FiniteVolumeEulerPerField
{
    constexpr static auto dimension = GridLayout::dimension;

public:
    FiniteVolumeEulerPerField(GridLayout const& layout, double const dt)
        : layout_{layout}
        , dt_{dt}
    {
    }

    template<typename Field, typename... Fluxes>
    void operator()(Field const& U, Field& Unew, Fluxes const&... fluxes) const
    {
        constexpr auto nFluxes = sizeof...(Fluxes);
        std::array<std::array<QtyCentering, dimension>, nFluxes> const centerings{
            layout_.centering(fluxes.physicalQuantity())...};
        std::array<double, nFluxes> coefs;
        for (std::size_t d = 0; d < nFluxes; ++d)
            coefs[d] = dt_ * layout_.inverseMeshSize(static_cast<Direction>(d));

        // each flux F_d is read at the cell's index and at its next index along d
        layout_.evalOnLines(Unew, [&](MeshIndex<dimension> const& first, std::uint32_t const n) {
            std::array<double const*, nFluxes> here, next;
            for_N<nFluxes>([&](auto d) {
                auto const& F = std::get<d>(std::forward_as_tuple(fluxes...));
                auto shifted  = first;
                shifted[d]    = layout_.nextIndex(centerings[d][d], first[d]);
                here[d]       = &F(first);
                next[d]       = &F(shifted);
            });
            double const* u = &U(first);
            double* unew    = &Unew(first);

            PHARE_VECTORIZE_IVDEP
            for (std::uint32_t k = 0; k < n; ++k)
            {
                double v = u[k];
                for (std::size_t d = 0; d < nFluxes; ++d)
                    v = v - coefs[d] * (next[d][k] - here[d][k]);
                unew[k] = v;
            }
        });
    }

    // as operator(), with the fluxes A_d + w B_d formed face by face (pairs {A_d, B_d})
    template<typename Field, typename... Pairs>
    void weighted(Field const& U, Field& Unew, double const w, Pairs const&... pairs) const
    {
        constexpr auto nFluxes = sizeof...(Pairs);
        std::array<std::array<QtyCentering, dimension>, nFluxes> const centerings{
            layout_.centering(pairs.first.physicalQuantity())...};
        std::array<double, nFluxes> coefs;
        for (std::size_t d = 0; d < nFluxes; ++d)
            coefs[d] = dt_ * layout_.inverseMeshSize(static_cast<Direction>(d));

        layout_.evalOnLines(Unew, [&](MeshIndex<dimension> const& first, std::uint32_t const n) {
            std::array<double const*, nFluxes> hereA, nextA, hereB, nextB;
            for_N<nFluxes>([&](auto d) {
                auto const& [A, B] = std::get<d>(std::forward_as_tuple(pairs...));
                auto shifted       = first;
                shifted[d]         = layout_.nextIndex(centerings[d][d], first[d]);
                hereA[d]           = &A(first);
                nextA[d]           = &A(shifted);
                hereB[d]           = &B(first);
                nextB[d]           = &B(shifted);
            });
            double const* u = &U(first);
            double* unew    = &Unew(first);

            PHARE_VECTORIZE_IVDEP
            for (std::uint32_t k = 0; k < n; ++k)
            {
                double v = u[k];
                for (std::size_t d = 0; d < nFluxes; ++d)
                    v = v
                        - coefs[d]
                              * ((nextA[d][k] + nextB[d][k] * w) - (hereA[d][k] + hereB[d][k] * w));
                unew[k] = v;
            }
        });
    }

private:
    GridLayout layout_;
    double dt_;
};

} // namespace PHARE::core

#endif
