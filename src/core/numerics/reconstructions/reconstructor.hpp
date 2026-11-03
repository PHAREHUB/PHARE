#ifndef PHARE_CORE_NUMERICS_RECONSTRUCTIONS_RECONSTRUCTOR_HPP
#define PHARE_CORE_NUMERICS_RECONSTRUCTIONS_RECONSTRUCTOR_HPP

#include "core/models/quantities/mhd_quantities.hpp"
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

        auto [BL, BR] = transverse_reconstruct<direction>(S.B, index);

        PerIndex uL{rhoL, {VxL, VyL, VzL}, BL, PL};
        PerIndex uR{rhoR, {VxR, VyR, VzR}, BR, PR};

        return std::make_pair(uL, uR);
    }

    template<auto direction, auto srcX, auto srcY, auto srcZ, auto dst, typename VecField>
    static auto center_reconstruct(VecField const& U, MeshIndex<VecField::dimension> index)
    {
        auto const& Ux = U(Component::X);
        auto const& Uy = U(Component::Y);
        auto const& Uz = U(Component::Z);

        auto [UxL, UxR]
            = Reconstruction::template center_reconstruct<direction, srcX, dst>(Ux, index);
        auto [UyL, UyR]
            = Reconstruction::template center_reconstruct<direction, srcY, dst>(Uy, index);
        auto [UzL, UzR]
            = Reconstruction::template center_reconstruct<direction, srcZ, dst>(Uz, index);

        return std::make_tuple(PerIndexVector{UxL, UyL, UzL}, PerIndexVector{UxR, UyR, UzR});
    }

    template<auto direction>
    static constexpr auto faceQuantity()
    {
        if constexpr (direction == Direction::X)
            return MHDQuantity::Scalar::Bx;
        else if constexpr (direction == Direction::Y)
            return MHDQuantity::Scalar::By;
        else if constexpr (direction == Direction::Z)
            return MHDQuantity::Scalar::Bz;
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
            direction, Reconstructor::template faceQuantity<transverse[0]>(),
            MHDQuantity::Scalar::rho>(Bt0, index);
        auto [Bt1L, Bt1R] = Reconstruction::template center_reconstruct<
            direction, Reconstructor::template faceQuantity<transverse[1]>(),
            MHDQuantity::Scalar::rho>(Bt1, index);

        PerIndexVector<typename VecField::value_type> BL, BR;
        BL(direction)     = Bn(index);
        BR(direction)     = Bn(index);
        BL(transverse[0]) = Bt0L;
        BR(transverse[0]) = Bt0R;
        BL(transverse[1]) = Bt1L;
        BR(transverse[1]) = Bt1R;

        return std::make_pair(BL, BR);
    }
};
} // namespace PHARE::core

#endif
