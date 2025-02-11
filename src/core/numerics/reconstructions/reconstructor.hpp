#ifndef PHARE_CORE_NUMERICS_RECONSTRUCTIONS_RECONSTRUCTOR_HPP
#define PHARE_CORE_NUMERICS_RECONSTRUCTIONS_RECONSTRUCTOR_HPP

#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"
#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include <utility>

namespace PHARE::core
{
template<typename Reconstruction>
struct Reconstructor
{
public:
    using GridLayout = Reconstruction::GridLayout_t;

    template<auto direction, typename State>
    static auto reconstruct(State const& S, MeshIndex<GridLayout::dimension> index)
    {
        auto [rhoL, rhoR] = Reconstruction::template reconstruct<direction>(S.rho, index);
        auto [VxL, VxR] = Reconstruction::template reconstruct<direction>(S.V(Component::X), index);
        auto [VyL, VyR] = Reconstruction::template reconstruct<direction>(S.V(Component::Y), index);
        auto [VzL, VzR] = Reconstruction::template reconstruct<direction>(S.V(Component::Z), index);
        auto [PL, PR]   = Reconstruction::template reconstruct<direction>(S.P, index);

        auto [BxL, ByL, BzL, BxR, ByR, BzR] = center_reconstruct<direction>(
            S.B, GridLayout::faceXToCellCenter(), GridLayout::faceYToCellCenter(),
            GridLayout::faceZToCellCenter(), index);

        PerIndex uL{rhoL, VxL, VyL, VzL, BxL, ByL, BzL, PL};
        PerIndex uR{rhoR, VxR, VyR, VzR, BxR, ByR, BzR, PR};

        return std::make_pair(uL, uR);
    }

    template<auto direction, typename VecField>
    static auto center_reconstruct(const VecField& U, auto projectionX, auto projectionY,
                                   auto projectionZ, MeshIndex<VecField::dimension> index)
    {
        auto const& Ux = U(Component::X);
        auto const& Uy = U(Component::Y);
        auto const& Uz = U(Component::Z);

        auto [UxL, UxR]
            = Reconstruction::template center_reconstruct<direction>(Ux, index, projectionX);
        auto [UyL, UyR]
            = Reconstruction::template center_reconstruct<direction>(Uy, index, projectionY);
        auto [UzL, UzR]
            = Reconstruction::template center_reconstruct<direction>(Uz, index, projectionZ);

        return std::make_tuple(UxL, UyL, UzL, UxR, UyR, UzR);
    }

    template<auto direction, typename VecField>
    static auto reconstructed_laplacian(auto inverseMeshSize, VecField const& J,
                                        MeshIndex<VecField::dimension> index)
    {
        auto const& Jx = J(Component::X);
        auto const& Jy = J(Component::Y);
        auto const& Jz = J(Component::Z);

        auto&& [laplJxL, laplJxR] = reconstructed_laplacian_component_<direction>(
            inverseMeshSize, Jx, index, GridLayout::edgeXToCellCenter());

        auto&& [laplJyL, laplJyR] = reconstructed_laplacian_component_<direction>(
            inverseMeshSize, Jy, index, GridLayout::edgeYToCellCenter());

        auto&& [laplJzL, laplJzR] = reconstructed_laplacian_component_<direction>(
            inverseMeshSize, Jz, index, GridLayout::edgeZToCellCenter());

        return std::make_tuple(laplJxL, laplJyL, laplJzL, laplJxR, laplJyR, laplJzR);
    }

private:
    template<auto direction, typename Field>
    static auto reconstructed_laplacian_component_(auto inverseMeshSize, Field const& J,
                                                   MeshIndex<Field::dimension> index,
                                                   auto projection)
    {
        auto d2 = [&](auto dir, auto const& prevValue, auto const& Value, auto const& nextValue) {
            return (inverseMeshSize[dir]) * (inverseMeshSize[dir])
                   * (prevValue - 2.0 * Value + nextValue);
        };

        auto [JL, JR]
            = Reconstruction::template center_reconstruct<direction>(J, index, projection);

        MeshIndex<Field::dimension> prevX = GridLayout::template previous<Direction::X>(index);
        MeshIndex<Field::dimension> nextX = GridLayout::template next<Direction::X>(index);

        auto [JL_X_1, JR_X_1]
            = Reconstruction::template center_reconstruct<direction>(J, prevX, projection);
        auto [JL_X1, JR_X1]
            = Reconstruction::template center_reconstruct<direction>(J, nextX, projection);

        std::uint32_t dirX = static_cast<std::uint32_t>(Direction::X);

        if constexpr (Field::dimension == 1)
        {
            auto LaplJL = d2(dirX, JL_X_1, JL, JL_X1);
            auto LaplJR = d2(dirX, JR_X_1, JR, JR_X1);

            return std::make_tuple(LaplJL, LaplJR);
        }
        else if (Field::dimension >= 2)
        {
            MeshIndex<Field::dimension> prevY = GridLayout::template previous<Direction::Y>(index);
            MeshIndex<Field::dimension> nextY = GridLayout::template next<Direction::Y>(index);

            auto [JL_Y_1, JR_Y_1]
                = Reconstruction::template center_reconstruct<direction>(J, prevY, projection);
            auto [JL_Y1, JR_Y1]
                = Reconstruction::template center_reconstruct<Direction::Y>(J, nextY, projection);

            std::uint32_t dirY = static_cast<std::uint32_t>(Direction::Y);

            if constexpr (Field::dimension == 2)
            {
                auto LaplJL = d2(dirX, JL_X_1, JL, JL_X1) + d2(dirY, JL_Y_1, JL, JL_Y1);
                auto LaplJR = d2(dirX, JR_X_1, JR, JR_X1) + d2(dirY, JR_Y_1, JR, JR_Y1);

                return std::make_tuple(LaplJL, LaplJR);
            }
            if constexpr (Field::dimension == 3)
            {
                MeshIndex<Field::dimension> prevZ
                    = GridLayout::template previous<Direction::Z>(index);
                MeshIndex<Field::dimension> nextZ = GridLayout::template next<Direction::Z>(index);

                auto [JL_Z_1, JR_Z_1]
                    = Reconstruction::template center_reconstruct<direction>(J, prevZ, projection);
                auto [JL_Z1, JR_Z1]
                    = Reconstruction::template center_reconstruct<direction>(J, nextZ, projection);

                std::uint32_t dirZ = static_cast<std::uint32_t>(Direction::Z);

                auto LaplJL = d2(dirX, JL_X_1, JL, JL_X1) + d2(dirY, JL_Y_1, JL, JL_Y1)
                              + d2(dirZ, JL_Z_1, JL, JL_Z1);
                auto LaplJR = d2(dirX, JR_X_1, JR, JR_X1) + d2(dirY, JR_Y_1, JR, JR_Y1)
                              + d2(dirZ, JR_Z_1, JR, JR_Z1);

                return std::make_tuple(LaplJL, LaplJR);
            }
        }
    }
};
} // namespace PHARE::core

#endif
