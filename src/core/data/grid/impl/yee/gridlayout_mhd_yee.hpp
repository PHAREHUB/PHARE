#ifndef PHARE_CORE_GRID_GRIDLAYOUTYEE_MHD_HPP
#define PHARE_CORE_GRID_GRIDLAYOUTYEE_MHD_HPP

#include "core/def.hpp"
#include "core/utilities/constants.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/models/quantities/mhd_quantities.hpp"


#include <array>

namespace PHARE::core::mhd
{

/**
 * @brief GridLayoutNdArrayImplYee class is a concrete GridLayoutNdArrayImpl used a Yee
 * type grid layout is needed.
 *
 * It provides methods related to grid layout operations:
 * - physical domain start/end indexes
 * - indexes of the first and last ghost nodes
 * - allocation sizes for Field attributes of other classes
 * - partial derivative operator (Faraday)
 * - physical coordinate given a field and a primal point (ix, iy, iz)
 * - cell centered coordinate given a primal point (ix, iy, iz)
 */
template<auto options>
class GridLayoutImplYee
{
    using Options = decltype(options);
    using Scalar  = Options::Scalar;

    static constexpr std::size_t dim = options.dimension;

public:
    struct GridData;
    static constexpr std::size_t dimension = dim;
    static constexpr std::string_view type = "yee";
    // The MHD layout reserves ghosts based on the reconstruction stencil width, plus the one
    // layer ampere loses by computing J on the shrinked ghost box.
    // static constexpr std::uint32_t reconstruction_nghosts = reconstruction_nghosts_;
    // Ghost width computed directly based on reconstruction stencil
    // static constexpr std::uint32_t ghost_width
    //     = nbrGhostsFromReconstruction<field_options.reconstruction_nghosts>();

    /**
     * @brief GridLayoutImpl<Selector<Layout,Layout::Yee>,dim>::initLayoutCentering_ initialize
     * the table MHDQuantityCentering_. This is THE important array in the GridLayout module.
     * This table knows which quantity is primal/dual along each direction. It is **this** array
     * that
     * **defines** what a Yee Layout is. Once this array is defined, the rest of the GridLayout
     * needs this array OK and can go on from here... hence all other functions in the Yee
     * interface are just calling private implementation common to all layouts
     */
    constexpr auto static initLayoutCentering_()
    {
        GridData const data{};

        std::array<QtyCentering, NBR_COMPO> const Rho = {{data.dual, data.dual, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const Vx = {{data.dual, data.dual, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const Vy = {{data.dual, data.dual, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const Vz = {{data.dual, data.dual, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const Bx = {{data.primal, data.dual, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const By = {{data.dual, data.primal, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const Bz = {{data.dual, data.dual, data.primal}};

        std::array<QtyCentering, NBR_COMPO> const P = {{data.dual, data.dual, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const rhoVx = {{data.dual, data.dual, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const rhoVy = {{data.dual, data.dual, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const rhoVz = {{data.dual, data.dual, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const Etot = {{data.dual, data.dual, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const Ex = {{data.dual, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Ey = {{data.primal, data.dual, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Ez = {{data.primal, data.primal, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const Jx = {{data.dual, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Jy = {{data.primal, data.dual, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Jz = {{data.primal, data.primal, data.dual}};


        std::array<QtyCentering, NBR_COMPO> const ScalarFlux_x
            = {{data.primal, data.dual, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const ScalarFlux_y
            = {{data.dual, data.primal, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const ScalarFlux_z
            = {{data.dual, data.dual, data.primal}};

        std::array<QtyCentering, NBR_COMPO> const VecFluxX_x
            = {{data.primal, data.dual, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const VecFluxY_x
            = {{data.primal, data.dual, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const VecFluxZ_x
            = {{data.primal, data.dual, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const VecFluxX_y
            = {{data.dual, data.primal, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const VecFluxY_y
            = {{data.dual, data.primal, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const VecFluxZ_y
            = {{data.dual, data.primal, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const VecFluxX_z
            = {{data.dual, data.dual, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const VecFluxY_z
            = {{data.dual, data.dual, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const VecFluxZ_z
            = {{data.dual, data.dual, data.primal}};

        std::array<QtyCentering, NBR_COMPO> const ScalarAllPrimal
            = {{data.primal, data.primal, data.primal}};

        std::array<QtyCentering, NBR_COMPO> const VecAllPrimalX
            = {{data.primal, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const VecAllPrimalY
            = {{data.primal, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const VecAllPrimalZ
            = {{data.primal, data.primal, data.primal}};

        auto const _QtyCentering = std::array{Rho,
                                              Vx,
                                              Vy,
                                              Vz,
                                              Bx,
                                              By,
                                              Bz,
                                              P,
                                              Etot,
                                              rhoVx,
                                              rhoVy,
                                              rhoVz,
                                              Ex,
                                              Ey,
                                              Ez,
                                              Jx,
                                              Jy,
                                              Jz,
                                              ScalarFlux_x,
                                              ScalarFlux_y,
                                              ScalarFlux_z,
                                              VecFluxX_x,
                                              VecFluxY_x,
                                              VecFluxZ_x,
                                              VecFluxX_y,
                                              VecFluxY_y,
                                              VecFluxZ_y,
                                              VecFluxX_z,
                                              VecFluxY_z,
                                              VecFluxZ_z,
                                              ScalarAllPrimal,
                                              VecAllPrimalX,
                                              VecAllPrimalY,
                                              VecAllPrimalZ};
        static_assert(_QtyCentering.size() == static_cast<std::size_t>(MHDQuantity::Scalar::count));

        return _QtyCentering;
    }

    //! says for each MHDQuantity::Quantity whether it is primal or dual, in each direction
    constexpr static std::array<std::array<QtyCentering, NBR_COMPO>,
                                static_cast<std::size_t>(MHDQuantity::Scalar::count)> const
        qtyCentering_{initLayoutCentering_()};

    static std::size_t const dim_{dim};

    //! qtyCentering_ restricted to the first dim directions, indexed by MHDQuantity::Scalar
    constexpr static auto scalarCentering_ = [] {
        std::array<std::array<QtyCentering, dim>, qtyCentering_.size()> centerings{};
        for (std::size_t iQty = 0; iQty < centerings.size(); ++iQty)
            for (std::size_t iDir = 0; iDir < dim; ++iDir)
                centerings[iQty][iDir] = qtyCentering_[iQty][iDir];
        return centerings;
    }();

    //! centerings of each vector's components, indexed by MHDQuantity::Vector
    constexpr static auto vectorCentering_ = [] {
        std::array<std::array<std::array<QtyCentering, dim>, 3>,
                   MHDQuantity::vector_components.size()>
            centerings{};
        for (std::size_t iVec = 0; iVec < centerings.size(); ++iVec)
            for (std::size_t iComp = 0; iComp < 3; ++iComp)
                centerings[iVec][iComp] = scalarCentering_[static_cast<std::size_t>(
                    MHDQuantity::vector_components[iVec][iComp])];
        return centerings;
    }();

    // ------------------------------------------------------------------------
    //                          PUBLIC INTERFACE
    // ------------------------------------------------------------------------
public:
    NO_DISCARD constexpr static std::array<QtyCentering, dim>
    centering(MHDQuantity::Scalar const qty)
    {
        auto const idx = static_cast<std::size_t>(qty);
        if (idx >= scalarCentering_.size())
            throw std::runtime_error("Wrong MHDQuantity");
        return scalarCentering_[idx];
    }

    NO_DISCARD constexpr static std::array<std::array<QtyCentering, dim>, 3>
    centering(MHDQuantity::Vector const qty)
    {
        auto const idx = static_cast<std::size_t>(qty);
        if (idx >= vectorCentering_.size())
            throw std::runtime_error("Wrong MHDQuantity");
        return vectorCentering_[idx];
    }

    enum class InterpDir { DualToPrimal = 0, PrimalToDual = 1 };

    // we might want to support the same interpolation possibilities as for the derivative, and
    // centralise the parametrisation of it
    template<auto dir, InterpDir interp_dir, std::size_t order = 2>
    NO_DISCARD static consteval auto directionalInterp()
    {
        constexpr int baseidx = (interp_dir == InterpDir::PrimalToDual);

        if constexpr (dir >= dimension)
        {
            return std::array{WeightPoint{Point<int, dimension>{}, 1.0}};
        }
        else
        {
            auto make_p = [](auto offset, auto d) {
                Point<int, dimension> p{};
                p[d] = offset;
                return p;
            };

            if constexpr (order == 2)
            {
                auto w = 0.5;
                return std::array{WeightPoint{make_p(baseidx - 1, dir), w},
                                  WeightPoint{make_p(baseidx, dir), w}};
            }
            else if constexpr (order == 4)
            {
                auto w1 = -1.0 / 16.0;
                auto w2 = 9.0 / 16.0;

                return std::array{WeightPoint{make_p(baseidx - 2, dir), w1},
                                  WeightPoint{make_p(baseidx - 1, dir), w2},
                                  WeightPoint{make_p(baseidx, dir), w2},
                                  WeightPoint{make_p(baseidx + 1, dir), w1}};
            }
            else if constexpr (order == 6)
            {
                auto w1 = 3.0 / 256.0;
                auto w2 = -25.0 / 256.0;
                auto w3 = 150.0 / 256.0;

                return std::array{WeightPoint{make_p(baseidx - 3, dir), w1},
                                  WeightPoint{make_p(baseidx - 2, dir), w2},
                                  WeightPoint{make_p(baseidx - 1, dir), w3},
                                  WeightPoint{make_p(baseidx, dir), w3},
                                  WeightPoint{make_p(baseidx + 1, dir), w2},
                                  WeightPoint{make_p(baseidx + 2, dir), w1}};
            }
        }
    }

    /**
     * @brief Primitive B — dual-direction prolongation (coarse→fine, ratio 2).
     *
     * A quantity dual in `dir` is an average over its dual extent there. Refining splits coarse
     * cell I into two fine children sitting at ±¼ of the coarse spacing; parity p∈{0,1},
     * sign σ = 2p−1 (σ=+1 = right child, σ=−1 = left). Offsets are relative to coarse cell I.
     *
     * Even-term cancellation: subtracting the conservation constraint kills every even
     * reconstruction coefficient, so the children
     * depend only on ODD coefficients ⇒ parabolic ≡ linear (degree-2 skipped). Children are
     * antisymmetric in σ about ū_I ⇒ they always mean back to ū_I (conservative at every order).
     *
     *   order 2 (linear, 3-pt):  ū_I + σ·(ū_{I+1} − ū_{I−1})/8
     */
    template<auto dir, int sign, std::size_t order = 2>
    NO_DISCARD static consteval auto directionalProlongation()
    {
        static_assert(sign == 1 || sign == -1, "child sign σ must be ±1");
        static_assert(order == 2, "dual prolongation ladder is order 2 (degree-2 skipped)");

        if constexpr (dir >= dimension)
        {
            return std::array{WeightPoint{Point<int, dimension>{}, 1.0}};
        }
        else
        {
            auto make_p = [](int offset) {
                Point<int, dimension> p{};
                p[dir] = offset;
                return p;
            };

            // ū_I + σ·(ū_{I+1} − ū_{I−1})/8: the two neighbours enter with the same magnitude
            // and opposite signs, so the children mean back to ū_I.
            constexpr double slopeWeight = 1. / 8.;

            return std::array{WeightPoint{make_p(-1), -sign * slopeWeight},
                              WeightPoint{make_p(0), 1.0},
                              WeightPoint{make_p(1), sign * slopeWeight}};
        }
    }

    // we could possibly have a variadic version to factor out these 2 overload, but this might
    // complexify the code quite a bit. possibly not worth it
    template<auto dir1, auto dir2>
    static consteval auto tensorProduct(auto const& s1, auto const& s2)
    {
        constexpr auto N1 = std::tuple_size_v<std::decay_t<decltype(s1)>>;
        constexpr auto N2 = std::tuple_size_v<std::decay_t<decltype(s2)>>;
        static_assert(dir1 != dir2);

        std::array<WeightPoint<dimension>, N1 * N2> result{};
        std::size_t k = 0;
        for (auto const& p1 : s1)
            for (auto const& p2 : s2)
            {
                Point<int, dimension> pt{};
                if constexpr (dir1 < dimension)
                    pt[dir1] = p1.indexes[dir1];
                if constexpr (dir2 < dimension)
                    pt[dir2] = p2.indexes[dir2];
                result[k++] = {pt, p1.coef * p2.coef};
            }
        return result;
    }

    template<auto dir1, auto dir2, auto dir3>
    static consteval auto tensorProduct(auto const& s1, auto const& s2, auto const& s3)
    {
        constexpr auto N1 = std::tuple_size_v<std::decay_t<decltype(s1)>>;
        constexpr auto N2 = std::tuple_size_v<std::decay_t<decltype(s2)>>;
        constexpr auto N3 = std::tuple_size_v<std::decay_t<decltype(s3)>>;
        static_assert(dir1 != dir2 && dir1 != dir3 && dir2 != dir3);

        std::array<WeightPoint<dimension>, N1 * N2 * N3> result{};
        std::size_t k = 0;
        for (auto const& p1 : s1)
            for (auto const& p2 : s2)
                for (auto const& p3 : s3)
                {
                    Point<int, dimension> pt{};
                    if constexpr (dir1 < dimension)
                        pt[dir1] = p1.indexes[dir1];
                    if constexpr (dir2 < dimension)
                        pt[dir2] = p2.indexes[dir2];
                    if constexpr (dir3 < dimension)
                        pt[dir3] = p3.indexes[dir3];
                    result[k++] = {pt, p1.coef * p2.coef * p3.coef};
                }
        return result;
    }

    // These functions where none of the centerings are in common could actually be expensive in
    // 3d, because they will return a full cubic stencil. In MHD, they are only used for spatial
    // hyper-resistivity, which could thus become quite expensive especially for high order.
    // This should be tested.
    NO_DISCARD auto static consteval BxToEx()
    {
        // Bx is primal dual dual
        // Ex is dual primal primal
        // operation is pdd to dpp

        using PHARE::core::dirX;
        using PHARE::core::dirY;
        using PHARE::core::dirZ;

        // even tho the directional interp correcly handles lower dimensions at compile time, we
        // could still want to have a swich on the dim to be nicer to the compiler. tbd
        return tensorProduct<dirX, dirY, dirZ>(directionalInterp<dirX, InterpDir::PrimalToDual>(),
                                               directionalInterp<dirY, InterpDir::DualToPrimal>(),
                                               directionalInterp<dirZ, InterpDir::DualToPrimal>());
    }

    NO_DISCARD auto static consteval ByToEx()
    {
        // By is dual primal dual
        // Ex is dual primal primal
        // operation is thus dpD to dpP
        // shift only in the Z direction

        using PHARE::core::dirZ;

        return directionalInterp<dirZ, InterpDir::DualToPrimal>();
    }

    NO_DISCARD auto static consteval BzToEx()
    {
        // Bz is dual dual primal
        // Ex is dual primal primal
        // operation is thus pDp to pPp
        // shift only in the Y direction

        using PHARE::core::dirY;

        return directionalInterp<dirY, InterpDir::DualToPrimal>();
    }

    NO_DISCARD auto static consteval BxToEy()
    {
        // Bx is primal dual dual
        // Ey is primal dual primal
        // operation is thus pdD to pdP
        // shift only in the Z direction

        using PHARE::core::dirZ;

        return directionalInterp<dirZ, InterpDir::DualToPrimal>();
    }

    NO_DISCARD auto static consteval ByToEy()
    {
        // By is dual primal dual
        // Ey is primal dual primal
        // the operation is thus dpd to pdp

        using PHARE::core::dirX;
        using PHARE::core::dirY;
        using PHARE::core::dirZ;

        return tensorProduct<dirX, dirY, dirZ>(directionalInterp<dirX, InterpDir::DualToPrimal>(),
                                               directionalInterp<dirY, InterpDir::PrimalToDual>(),
                                               directionalInterp<dirZ, InterpDir::DualToPrimal>());
    }

    NO_DISCARD auto static consteval BzToEy()
    {
        // Bz is dual dual primal
        // Ey is primal dual primal
        // operation is thus Ddp to Pdp
        // shift only in the X direction

        using PHARE::core::dirX;

        return directionalInterp<dirX, InterpDir::DualToPrimal>();
    }

    NO_DISCARD auto static consteval BzToEz()
    {
        // Bz is dual dual primal
        // Ez is primal primal dual
        // operation is thus ddp to ppd

        using PHARE::core::dirX;
        using PHARE::core::dirY;
        using PHARE::core::dirZ;

        return tensorProduct<dirX, dirY, dirZ>(directionalInterp<dirX, InterpDir::DualToPrimal>(),
                                               directionalInterp<dirY, InterpDir::DualToPrimal>(),
                                               directionalInterp<dirZ, InterpDir::PrimalToDual>());
    }

    NO_DISCARD auto static consteval ByToEz()
    {
        // By is dual primal dual
        // Ez is primal primal dual
        // operation is thus Dpd to Ppd
        // shift only in the X direction

        using PHARE::core::dirX;

        return directionalInterp<dirX, InterpDir::DualToPrimal>();
    }

    NO_DISCARD auto static consteval BxToEz()
    {
        // Bx is primal dual dual
        // Ez is primal primal dual
        // operation is thus pDd to pPd
        // shift only in the Y direction

        using PHARE::core::dirY;

        return directionalInterp<dirY, InterpDir::DualToPrimal>();
    }

    NO_DISCARD auto static constexpr cellCenterToEdgeX()
    {
        // cellcenter is dual dual dual
        // edgeX is dual primal primal
        // operation is thus dDD to dPP
        // shift in the YZ direction

        using PHARE::core::dirY;
        using PHARE::core::dirZ;

        return tensorProduct<dirY, dirZ>(directionalInterp<dirY, InterpDir::DualToPrimal>(),
                                         directionalInterp<dirZ, InterpDir::DualToPrimal>());
    }

    NO_DISCARD auto static constexpr cellCenterToEdgeY()
    {
        // cellcenter is dual dual dual
        // edgeY is primal dual primal
        // operation is thus DdD to PdP
        // shift in the XZ direction

        using PHARE::core::dirX;
        using PHARE::core::dirZ;

        return tensorProduct<dirX, dirZ>(directionalInterp<dirX, InterpDir::DualToPrimal>(),
                                         directionalInterp<dirZ, InterpDir::DualToPrimal>());
    }

    NO_DISCARD auto static constexpr cellCenterToEdgeZ()
    {
        // cellcenter is dual dual dual
        // edgeZ is primal primal dual
        // operation is thus DDd to PPd
        // shift in the XY direction

        using PHARE::core::dirX;
        using PHARE::core::dirY;

        return tensorProduct<dirX, dirY>(directionalInterp<dirX, InterpDir::DualToPrimal>(),
                                         directionalInterp<dirY, InterpDir::DualToPrimal>());
    }

    NO_DISCARD auto static constexpr faceXToCellCenter()
    {
        // The X face is Pdd
        // the mhd quantities in FV are Ddd
        // operation is thus Pdd to Ddd
        // shift only in the X direction

        using PHARE::core::dirX;

        return directionalInterp<dirX, InterpDir::PrimalToDual>();
    }

    NO_DISCARD auto static constexpr faceYToCellCenter()
    {
        // The Y face is Dpd
        // the mhd quantities in FV are Ddd
        // operation is thus Dpd to Ddd
        // shift only in the Y direction

        using PHARE::core::dirY;

        return directionalInterp<dirY, InterpDir::PrimalToDual>();
    }

    NO_DISCARD auto static constexpr faceZToCellCenter()
    {
        // The Z face is Ddp
        // the mhd quantities in FV are Ddd
        // operation is thus Ddp to Ddd
        // shift only in the Z direction

        using PHARE::core::dirZ;

        return directionalInterp<dirZ, InterpDir::PrimalToDual>();
    }

    NO_DISCARD auto static constexpr edgeYToFaceX()
    {
        // edge-Y is primal dual primal
        // face-X is primal dual dual
        // operation is thus pdP to pdD, shift only in Z

        using PHARE::core::dirZ;

        return directionalInterp<dirZ, InterpDir::PrimalToDual>();
    }

    NO_DISCARD auto static constexpr edgeZToFaceX()
    {
        // edge-Z is primal primal dual
        // face-X is primal dual dual
        // operation is thus pPd to pDd, shift only in Y

        using PHARE::core::dirY;

        return directionalInterp<dirY, InterpDir::PrimalToDual>();
    }

    NO_DISCARD auto static constexpr edgeXToFaceY()
    {
        // edge-X is dual primal primal
        // face-Y is dual primal dual
        // operation is thus dpP to dpD, shift only in Z

        using PHARE::core::dirZ;

        return directionalInterp<dirZ, InterpDir::PrimalToDual>();
    }

    NO_DISCARD auto static constexpr edgeZToFaceY()
    {
        // edge-Z is primal primal dual
        // face-Y is dual primal dual
        // operation is thus Ppd to Dpd, shift only in X

        using PHARE::core::dirX;

        return directionalInterp<dirX, InterpDir::PrimalToDual>();
    }

    NO_DISCARD auto static constexpr edgeXToFaceZ()
    {
        // edge-X is dual primal primal
        // face-Z is dual dual primal
        // operation is thus dPp to dDp, shift only in Y

        using PHARE::core::dirY;

        return directionalInterp<dirY, InterpDir::PrimalToDual>();
    }

    NO_DISCARD auto static constexpr edgeYToFaceZ()
    {
        // edge-Y is primal dual primal
        // face-Z is dual dual primal
        // operation is thus Pdp to Ddp, shift only in X

        using PHARE::core::dirX;

        return directionalInterp<dirX, InterpDir::PrimalToDual>();
    }

    NO_DISCARD auto static constexpr ByToFaceX()
    {
        // By is dual primal dual
        // face-X is primal dual dual

        using PHARE::core::dirX;
        using PHARE::core::dirY;

        return tensorProduct<dirX, dirY>(directionalInterp<dirX, InterpDir::DualToPrimal>(),
                                         directionalInterp<dirY, InterpDir::PrimalToDual>());
    }

    NO_DISCARD auto static constexpr BzToFaceX()
    {
        // Bz is dual dual primal
        // face-X is primal dual dual

        using PHARE::core::dirX;
        using PHARE::core::dirZ;

        return tensorProduct<dirX, dirZ>(directionalInterp<dirX, InterpDir::DualToPrimal>(),
                                         directionalInterp<dirZ, InterpDir::PrimalToDual>());
    }

    NO_DISCARD auto static constexpr BxToFaceY()
    {
        // Bx is primal dual dual
        // face-Y is dual primal dual

        using PHARE::core::dirX;
        using PHARE::core::dirY;

        return tensorProduct<dirX, dirY>(directionalInterp<dirX, InterpDir::PrimalToDual>(),
                                         directionalInterp<dirY, InterpDir::DualToPrimal>());
    }

    NO_DISCARD auto static constexpr BzToFaceY()
    {
        // Bz is dual dual primal
        // face-Y is dual primal dual

        using PHARE::core::dirY;
        using PHARE::core::dirZ;

        return tensorProduct<dirY, dirZ>(directionalInterp<dirY, InterpDir::DualToPrimal>(),
                                         directionalInterp<dirZ, InterpDir::PrimalToDual>());
    }

    NO_DISCARD auto static constexpr BxToFaceZ()
    {
        // Bx is primal dual dual
        // face-Z is dual dual primal

        using PHARE::core::dirX;
        using PHARE::core::dirZ;

        return tensorProduct<dirX, dirZ>(directionalInterp<dirX, InterpDir::PrimalToDual>(),
                                         directionalInterp<dirZ, InterpDir::DualToPrimal>());
    }

    NO_DISCARD auto static constexpr ByToFaceZ()
    {
        // By is dual primal dual
        // face-Z is dual dual primal

        using PHARE::core::dirY;
        using PHARE::core::dirZ;

        return tensorProduct<dirY, dirZ>(directionalInterp<dirY, InterpDir::PrimalToDual>(),
                                         directionalInterp<dirZ, InterpDir::DualToPrimal>());
    }

    NO_DISCARD auto static constexpr edgeXToCellCenter()
    {
        // The X edge is dPP
        // the mhd quantities in FV are dDD
        // operation is thus dPP to dDD

        using PHARE::core::dirY;
        using PHARE::core::dirZ;

        return tensorProduct<dirY, dirZ>(directionalInterp<dirY, InterpDir::PrimalToDual>(),
                                         directionalInterp<dirZ, InterpDir::PrimalToDual>());
    }

    NO_DISCARD auto static constexpr edgeYToCellCenter()
    {
        // The Y edge is PdP
        // the mhd quantities in FV are DdD
        // operation is thus PdP to DdD

        using PHARE::core::dirX;
        using PHARE::core::dirZ;

        return tensorProduct<dirX, dirZ>(directionalInterp<dirX, InterpDir::PrimalToDual>(),
                                         directionalInterp<dirZ, InterpDir::PrimalToDual>());
    }

    NO_DISCARD auto static constexpr edgeZToCellCenter()
    {
        // The Z edge is PPd
        // the mhd quantities in FV are DDd
        // operation is thus PPd to DDd

        using PHARE::core::dirX;
        using PHARE::core::dirY;

        return tensorProduct<dirX, dirY>(directionalInterp<dirX, InterpDir::PrimalToDual>(),
                                         directionalInterp<dirY, InterpDir::PrimalToDual>());
    }

    NO_DISCARD auto static consteval BxToMoments()
    {
        // Bx is primal dual dual
        // moments are primal primal primal
        // operation is thus pDD to pPP

        using PHARE::core::dirY;
        using PHARE::core::dirZ;

        return tensorProduct<dirY, dirZ>(directionalInterp<dirY, InterpDir::DualToPrimal, 2>(),
                                         directionalInterp<dirZ, InterpDir::DualToPrimal, 2>());
    }



    NO_DISCARD auto static consteval ByToMoments()
    {
        // By is dual primal dual
        // moments are primal primal primal
        // operation is thus DpD to PpP

        using PHARE::core::dirX;
        using PHARE::core::dirZ;

        return tensorProduct<dirX, dirZ>(directionalInterp<dirX, InterpDir::DualToPrimal, 2>(),
                                         directionalInterp<dirZ, InterpDir::DualToPrimal, 2>());
    }




    NO_DISCARD auto static consteval BzToMoments()
    {
        // Bz is dual dual primal
        // moments are primal primal primal
        // operation is thus DDp to PPp

        using PHARE::core::dirX;
        using PHARE::core::dirY;

        return tensorProduct<dirX, dirY>(directionalInterp<dirX, InterpDir::DualToPrimal, 2>(),
                                         directionalInterp<dirY, InterpDir::DualToPrimal, 2>());
    }

    // We might not want too high order of a stencil for all of the data we have. Also these
    // interpolations assume point data, so we possibly need a different dump strategy for high
    // order mhd. I think defaulting to second order projections here for now is probably a good
    // call.
    NO_DISCARD auto static constexpr cellCenterToFullPrimal()
    {
        // operation is thus DDD to PPP

        using PHARE::core::dirX;
        using PHARE::core::dirY;
        using PHARE::core::dirZ;

        return tensorProduct<dirX, dirY, dirZ>(
            directionalInterp<dirX, InterpDir::DualToPrimal, 2>(),
            directionalInterp<dirY, InterpDir::DualToPrimal, 2>(),
            directionalInterp<dirZ, InterpDir::DualToPrimal, 2>());
    }
};



template<auto options>
struct GridLayoutImplYee<options>::GridData
{
    static constexpr Direction dirX = Direction::X;
    static constexpr Direction dirY = Direction::Y;
    static constexpr Direction dirZ = Direction::Z;

    static constexpr QtyCentering primal = QtyCentering::primal;
    static constexpr QtyCentering dual   = QtyCentering::dual;

    static constexpr std::uint32_t idirX = static_cast<std::uint32_t>(Direction::X);
    static constexpr std::uint32_t idirY = static_cast<std::uint32_t>(Direction::Y);
    static constexpr std::uint32_t idirZ = static_cast<std::uint32_t>(Direction::Z);

    static constexpr std::uint32_t irho = static_cast<std::uint32_t>(MHDQuantity::Scalar::rho);

    static constexpr std::uint32_t iVx = static_cast<std::uint32_t>(MHDQuantity::Scalar::Vx);
    static constexpr std::uint32_t iVy = static_cast<std::uint32_t>(MHDQuantity::Scalar::Vy);
    static constexpr std::uint32_t iVz = static_cast<std::uint32_t>(MHDQuantity::Scalar::Vz);

    static constexpr std::uint32_t iBx = static_cast<std::uint32_t>(MHDQuantity::Scalar::Bx);
    static constexpr std::uint32_t iBy = static_cast<std::uint32_t>(MHDQuantity::Scalar::By);
    static constexpr std::uint32_t iBz = static_cast<std::uint32_t>(MHDQuantity::Scalar::Bz);

    static constexpr std::uint32_t iP = static_cast<std::uint32_t>(MHDQuantity::Scalar::P);

    static constexpr std::uint32_t iEtot = static_cast<std::uint32_t>(MHDQuantity::Scalar::Etot);

    static constexpr std::uint32_t irhoVx = static_cast<std::uint32_t>(MHDQuantity::Scalar::rhoVx);
    static constexpr std::uint32_t irhoVy = static_cast<std::uint32_t>(MHDQuantity::Scalar::rhoVy);
    static constexpr std::uint32_t irhoVz = static_cast<std::uint32_t>(MHDQuantity::Scalar::rhoVz);

    static constexpr std::uint32_t iEx = static_cast<std::uint32_t>(MHDQuantity::Scalar::Ex);
    static constexpr std::uint32_t iEy = static_cast<std::uint32_t>(MHDQuantity::Scalar::Ey);
    static constexpr std::uint32_t iEz = static_cast<std::uint32_t>(MHDQuantity::Scalar::Ez);

    static constexpr std::uint32_t iJx = static_cast<std::uint32_t>(MHDQuantity::Scalar::Jx);
    static constexpr std::uint32_t iJy = static_cast<std::uint32_t>(MHDQuantity::Scalar::Jy);
    static constexpr std::uint32_t iJz = static_cast<std::uint32_t>(MHDQuantity::Scalar::Jz);

    static constexpr std::uint32_t iScalarFlux_x
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::ScalarFlux_x);
    static constexpr std::uint32_t iScalarFlux_y
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::ScalarFlux_y);
    static constexpr std::uint32_t iScalarFlux_z
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::ScalarFlux_z);

    static constexpr std::uint32_t iVecFluxX_x
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecFluxX_x);
    static constexpr std::uint32_t iVecFluxY_x
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecFluxY_x);
    static constexpr std::uint32_t iVecFluxZ_x
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecFluxZ_x);

    static constexpr std::uint32_t iVecFluxX_y
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecFluxX_y);
    static constexpr std::uint32_t iVecFluxY_y
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecFluxY_y);
    static constexpr std::uint32_t iVecFluxZ_y
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecFluxZ_y);

    static constexpr std::uint32_t iVecFluxX_z
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecFluxX_z);
    static constexpr std::uint32_t iVecFluxY_z
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecFluxY_z);
    static constexpr std::uint32_t iVecFluxZ_z
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecFluxZ_z);

    static constexpr std::uint32_t iScalarAllPrimal
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::ScalarAllPrimal);

    static constexpr std::uint32_t iVecAllPrimalX
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecAllPrimalX);
    static constexpr std::uint32_t iVecAllPrimalY
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecAllPrimalY);
    static constexpr std::uint32_t iVecAllPrimalZ
        = static_cast<std::uint32_t>(MHDQuantity::Scalar::VecAllPrimalZ);
};



} // namespace PHARE::core::mhd


#endif // PHARE_CORE_GRID_GRIDLAYOUTYEE_MHD_HPP
