#ifndef PHARE_CORE_NUMERICS_RECONSTRUCTIONS_RECONSTRUCTOR_HPP
#define PHARE_CORE_NUMERICS_RECONSTRUCTIONS_RECONSTRUCTOR_HPP

#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"
#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include <type_traits>
#include <utility>

namespace PHARE::core
{
template<typename Reconstruction>
struct Reconstructor
{
public:
    using GridLayout = Reconstruction::GridLayout_t;

    template<auto direction, typename State, typename ExternalB>
    static auto reconstruct(State const& S, ExternalB const& b0,
                            MeshIndex<GridLayout::dimension> index)
    {
        auto [rhoL, rhoR] = Reconstruction::template reconstruct<direction>(S.rho, index);
        auto [VxL, VxR] = Reconstruction::template reconstruct<direction>(S.V(Component::X), index);
        auto [VyL, VyR] = Reconstruction::template reconstruct<direction>(S.V(Component::Y), index);
        auto [VzL, VzR] = Reconstruction::template reconstruct<direction>(S.V(Component::Z), index);
        auto [PL, PR]   = Reconstruction::template reconstruct<direction>(S.P, index);

        // Only the perturbation B1 is reconstructed.
        auto [B1L, B1R] = transverse_reconstruct<direction>(S.B1, index);

        // B0 is not reconstructed: it is read once at the Riemann face (single value) and stored
        // alongside the reconstructed perturbation B1. The total field is formed by addition only,
        // where needed. This keeps B0 out of the Riemann jump (B1R - B1L) and makes a static
        // equilibrium well-balanced.
        auto const B0f = B0_at_face<direction>(b0, index);

        using Float = std::decay_t<decltype(rhoL)>;

        PerIndex<Float> uL{rhoL, PerIndexVector<Float>{VxL, VyL, VzL}, B1L, PL, B0f};
        PerIndex<Float> uR{rhoR, PerIndexVector<Float>{VxR, VyR, VzR}, B1R, PR, B0f};

        return std::make_pair(uL, uR);
    }

    template<auto direction, auto ProjectionX, auto ProjectionY, auto ProjectionZ, typename VecField>
    static auto center_reconstruct(VecField const& U, MeshIndex<VecField::dimension> index)
    {
        auto const& Ux = U(Component::X);
        auto const& Uy = U(Component::Y);
        auto const& Uz = U(Component::Z);

        auto [UxL, UxR]
            = Reconstruction::template center_reconstruct<direction, ProjectionX>(Ux, index);
        auto [UyL, UyR]
            = Reconstruction::template center_reconstruct<direction, ProjectionY>(Uy, index);
        auto [UzL, UzR]
            = Reconstruction::template center_reconstruct<direction, ProjectionZ>(Uz, index);

        return std::make_tuple(PerIndexVector{UxL, UyL, UzL}, PerIndexVector{UxR, UyR, UzR});
    }

    template<auto direction>
    static constexpr auto projection()
    {
        if constexpr (direction == Direction::X)
            return GridLayout::faceXToCellCenter();
        else if constexpr (direction == Direction::Y)
            return GridLayout::faceYToCellCenter();
        else if constexpr (direction == Direction::Z)
            return GridLayout::faceZToCellCenter();
    }

    // The normal direction for B is already face centered, so we only reconstruct the transverse
    template<auto direction, typename VecField>
    static auto transverse_reconstruct(VecField const& B, MeshIndex<VecField::dimension> index)
    {
        auto constexpr transverse = []() { // probably should be a util function somewhere
            if constexpr (direction == Direction::X)
                return std::array{Direction::Y, Direction::Z};
            else if constexpr (direction == Direction::Y)
                return std::array{Direction::X, Direction::Z};
            else if constexpr (direction == Direction::Z)
                return std::array{Direction::X, Direction::Y};
        }();


        auto const Bn  = B(static_cast<Component>(direction));
        auto const Bt0 = B(static_cast<Component>(transverse[0]));
        auto const Bt1 = B(static_cast<Component>(transverse[1]));

        auto [Bt0L, Bt0R] = Reconstruction::template center_reconstruct<
            direction, Reconstructor::template projection<transverse[0]>>(Bt0, index);
        auto [Bt1L, Bt1R] = Reconstruction::template center_reconstruct<
            direction, Reconstructor::template projection<transverse[1]>>(Bt1, index);

        PerIndexVector<typename VecField::value_type> BL, BR;
        BL(direction)     = Bn(index);
        BR(direction)     = Bn(index);
        BL(transverse[0]) = Bt0L;
        BR(transverse[0]) = Bt0R;
        BL(transverse[1]) = Bt1L;
        BR(transverse[1]) = Bt1R;

        return std::make_pair(BL, BR);
    }

    // B0 is not reconstructed: it is read once at the Riemann face (single value). The normal
    // component is naturally face-centered (read directly from B0); the transverse components are
    // linearly interpolated from the native face-centered B0 onto this face. B0 stays out of the
    // Riemann jump (single value, same on L/R), keeping a static equilibrium well-balanced.
    template<auto direction, typename VecField>
    static auto B0_at_face(VecField const& b0, MeshIndex<GridLayout::dimension> index)
    {
        using value_type = typename std::decay_t<decltype(b0(Component::X))>::value_type;
        PerIndexVector<value_type> B0f;
        B0f(direction) = b0(static_cast<Component>(direction))(index);
        if constexpr (direction == Direction::X)
        {
            B0f(Component::Y)
                = GridLayout::template project<GridLayout::ByToFaceX>(b0(Component::Y), index);
            B0f(Component::Z)
                = GridLayout::template project<GridLayout::BzToFaceX>(b0(Component::Z), index);
        }
        else if constexpr (direction == Direction::Y)
        {
            B0f(Component::X)
                = GridLayout::template project<GridLayout::BxToFaceY>(b0(Component::X), index);
            B0f(Component::Z)
                = GridLayout::template project<GridLayout::BzToFaceY>(b0(Component::Z), index);
        }
        else // Direction::Z
        {
            B0f(Component::X)
                = GridLayout::template project<GridLayout::BxToFaceZ>(b0(Component::X), index);
            B0f(Component::Y)
                = GridLayout::template project<GridLayout::ByToFaceZ>(b0(Component::Y), index);
        }
        return B0f;
    }

    template<auto direction, typename VecField>
    static auto reconstructed_laplacian(auto inverseMeshSize, VecField const& J,
                                        MeshIndex<VecField::dimension> index)
    {
        auto const& Jx = J(Component::X);
        auto const& Jy = J(Component::Y);
        auto const& Jz = J(Component::Z);

        auto const& [laplJxL, laplJxR]
            = reconstructed_laplacian_component_<direction, GridLayout::edgeXToCellCenter>(
                inverseMeshSize, Jx, index);

        auto const& [laplJyL, laplJyR]
            = reconstructed_laplacian_component_<direction, GridLayout::edgeYToCellCenter>(
                inverseMeshSize, Jy, index);

        auto const& [laplJzL, laplJzR]
            = reconstructed_laplacian_component_<direction, GridLayout::edgeZToCellCenter>(
                inverseMeshSize, Jz, index);

        return std::make_tuple(PerIndexVector{laplJxL, laplJyL, laplJzL},
                               PerIndexVector{laplJxR, laplJyR, laplJzR});
    }

private:
    template<auto direction, auto Projection, typename Field>
    static auto reconstructed_laplacian_component_(auto inverseMeshSize, Field const& J,
                                                   MeshIndex<Field::dimension> index)
    {
        auto d2 = [&](auto dir, auto const& prevValue, auto const& Value, auto const& nextValue) {
            return (inverseMeshSize[dir]) * (inverseMeshSize[dir])
                   * (prevValue - 2.0 * Value + nextValue);
        };

        auto const [JL, JR]
            = Reconstruction::template center_reconstruct<direction, Projection>(J, index);

        MeshIndex<Field::dimension> prevX = GridLayout::template previous<Direction::X>(index);
        MeshIndex<Field::dimension> nextX = GridLayout::template next<Direction::X>(index);

        auto const [JL_X_1, JR_X_1]
            = Reconstruction::template center_reconstruct<direction, Projection>(J, prevX);
        auto const [JL_X1, JR_X1]
            = Reconstruction::template center_reconstruct<direction, Projection>(J, nextX);

        std::uint32_t dirX = static_cast<std::uint32_t>(Direction::X);

        if constexpr (Field::dimension == 1)
        {
            auto const LaplJL = d2(dirX, JL_X_1, JL, JL_X1);
            auto const LaplJR = d2(dirX, JR_X_1, JR, JR_X1);

            return std::make_tuple(LaplJL, LaplJR);
        }
        else if constexpr (Field::dimension >= 2)
        {
            MeshIndex<Field::dimension> prevY = GridLayout::template previous<Direction::Y>(index);
            MeshIndex<Field::dimension> nextY = GridLayout::template next<Direction::Y>(index);

            auto const [JL_Y_1, JR_Y_1]
                = Reconstruction::template center_reconstruct<direction, Projection>(J, prevY);
            auto const [JL_Y1, JR_Y1]
                = Reconstruction::template center_reconstruct<direction, Projection>(J, nextY);

            std::uint32_t dirY = static_cast<std::uint32_t>(Direction::Y);

            if constexpr (Field::dimension == 2)
            {
                auto const LaplJL = d2(dirX, JL_X_1, JL, JL_X1) + d2(dirY, JL_Y_1, JL, JL_Y1);
                auto const LaplJR = d2(dirX, JR_X_1, JR, JR_X1) + d2(dirY, JR_Y_1, JR, JR_Y1);

                return std::make_tuple(LaplJL, LaplJR);
            }
            if constexpr (Field::dimension == 3)
            {
                MeshIndex<Field::dimension> prevZ
                    = GridLayout::template previous<Direction::Z>(index);
                MeshIndex<Field::dimension> nextZ = GridLayout::template next<Direction::Z>(index);

                auto const [JL_Z_1, JR_Z_1]
                    = Reconstruction::template center_reconstruct<direction, Projection>(J, prevZ);
                auto const [JL_Z1, JR_Z1]
                    = Reconstruction::template center_reconstruct<direction, Projection>(J, nextZ);

                std::uint32_t dirZ = static_cast<std::uint32_t>(Direction::Z);

                auto const LaplJL = d2(dirX, JL_X_1, JL, JL_X1) + d2(dirY, JL_Y_1, JL, JL_Y1)
                                    + d2(dirZ, JL_Z_1, JL, JL_Z1);
                auto const LaplJR = d2(dirX, JR_X_1, JR, JR_X1) + d2(dirY, JR_Y_1, JR, JR_Y1)
                                    + d2(dirZ, JR_Z_1, JR, JR_Z1);

                return std::make_tuple(LaplJL, LaplJR);
            }
        }
    }
};
} // namespace PHARE::core

#endif
