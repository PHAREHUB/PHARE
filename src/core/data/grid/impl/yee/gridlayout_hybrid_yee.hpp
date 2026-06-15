#ifndef PHARE_CORE_GRID_IMPL_YEE_GRIDLAYOUT_HYBRID_YEE_HPP
#define PHARE_CORE_GRID_IMPL_YEE_GRIDLAYOUT_HYBRID_YEE_HPP


#include "core/def.hpp"
#include "core/utilities/constants.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"


#include <array>


namespace PHARE::core::hybrid
{


/**
 * @brief GridData provides constants used to initialize:
 * - hybridQuantity centerings
 * - physical start/end indexes
 * - ghost start/end indexes
 * - numbers of padding cells and physical cells
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


    /**
     * @brief GridLayoutImpl<Selector<Layout,Layout::Yee>,dim>::initLayoutCentering_ initialize
     * the table hybridQuantityCentering_. This is THE important array in the GridLayout module.
     * This table knows which quantity is primal/dual along each direction. It is **this** array
     * that
     * **defines** what a Yee Layout is. Once this array is defined, the rest of the GridLayout
     * needs this array OK and can go on from here... hence all other functions in the Yee
     * interface are just calling private implementation common to all layouts
     */
    constexpr auto static initLayoutCentering_()
    {
        GridData const data{};
        std::array<QtyCentering, NBR_COMPO> const Bx = {{data.primal, data.dual, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const By = {{data.dual, data.primal, data.dual}};
        std::array<QtyCentering, NBR_COMPO> const Bz = {{data.dual, data.dual, data.primal}};

        std::array<QtyCentering, NBR_COMPO> const Ex = {{data.dual, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Ey = {{data.primal, data.dual, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Ez = {{data.primal, data.primal, data.dual}};


        std::array<QtyCentering, NBR_COMPO> const Jx = {{data.dual, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Jy = {{data.primal, data.dual, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Jz = {{data.primal, data.primal, data.dual}};

        std::array<QtyCentering, NBR_COMPO> const Rho = {{data.primal, data.primal, data.primal}};

        std::array<QtyCentering, NBR_COMPO> const Vx = {{data.primal, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Vy = {{data.primal, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Vz = {{data.primal, data.primal, data.primal}};

        std::array<QtyCentering, NBR_COMPO> const Mxx = {{data.primal, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Mxy = {{data.primal, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Mxz = {{data.primal, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Myy = {{data.primal, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Myz = {{data.primal, data.primal, data.primal}};
        std::array<QtyCentering, NBR_COMPO> const Mzz = {{data.primal, data.primal, data.primal}};

        std::array<QtyCentering, NBR_COMPO> const P = {{data.primal, data.primal, data.primal}};

        std::array<std::array<QtyCentering, NBR_COMPO>,
                   static_cast<std::size_t>(HybridQuantity::Scalar::count)> const
            hybridQtyCentering{Bx, By, Bz, Ex, Ey,  Ez,  Jx,  Jy,  Jz,  Rho,
                               Vx, Vy, Vz, P,  Mxx, Mxy, Mxz, Myy, Myz, Mzz};


        return hybridQtyCentering;
    }



    //! says for each HybridQuantity::Quantity whether it is primal or dual, in each direction
    constexpr static std::array<std::array<QtyCentering, NBR_COMPO>,
                                static_cast<std::size_t>(HybridQuantity::Scalar::count)> const
        qtyCentering_{initLayoutCentering_()};

    static std::size_t const dim_{dim};

    // ------------------------------------------------------------------------
    //                          PUBLIC INTERFACE
    // ------------------------------------------------------------------------
public:
    NO_DISCARD constexpr static std::array<QtyCentering, dim>
    centering(HybridQuantity::Scalar hybridQuantity)
    {
        constexpr GridData gridData_{};
        if constexpr (dim == 1)
        {
            switch (hybridQuantity)
            {
                case HybridQuantity::Scalar::Bx:
                    return {{qtyCentering_[gridData_.iBx][gridData_.idirX]}};
                case HybridQuantity::Scalar::By:
                    return {{qtyCentering_[gridData_.iBy][gridData_.idirX]}};
                case HybridQuantity::Scalar::Bz:
                    return {{qtyCentering_[gridData_.iBz][gridData_.idirX]}};
                case HybridQuantity::Scalar::Ex:
                    return {{qtyCentering_[gridData_.iEx][gridData_.idirX]}};
                case HybridQuantity::Scalar::Ey:
                    return {{qtyCentering_[gridData_.iEy][gridData_.idirX]}};
                case HybridQuantity::Scalar::Ez:
                    return {{qtyCentering_[gridData_.iEz][gridData_.idirX]}};
                case HybridQuantity::Scalar::Jx:
                    return {{qtyCentering_[gridData_.iJx][gridData_.idirX]}};
                case HybridQuantity::Scalar::Jy:
                    return {{qtyCentering_[gridData_.iJy][gridData_.idirX]}};
                case HybridQuantity::Scalar::Jz:
                    return {{qtyCentering_[gridData_.iJz][gridData_.idirX]}};
                case HybridQuantity::Scalar::rho:
                    return {{qtyCentering_[gridData_.irho][gridData_.idirX]}};
                case HybridQuantity::Scalar::Vx:
                    return {{qtyCentering_[gridData_.iVx][gridData_.idirX]}};
                case HybridQuantity::Scalar::Vy:
                    return {{qtyCentering_[gridData_.iVy][gridData_.idirX]}};
                case HybridQuantity::Scalar::Vz:
                    return {{qtyCentering_[gridData_.iVz][gridData_.idirX]}};
                case HybridQuantity::Scalar::P:
                    return {{qtyCentering_[gridData_.iP][gridData_.idirX]}};
                case HybridQuantity::Scalar::Mxx:
                    return {{qtyCentering_[gridData_.iMxx][gridData_.idirX]}};
                case HybridQuantity::Scalar::Mxy:
                    return {{qtyCentering_[gridData_.iMxy][gridData_.idirX]}};
                case HybridQuantity::Scalar::Mxz:
                    return {{qtyCentering_[gridData_.iMxz][gridData_.idirX]}};
                case HybridQuantity::Scalar::Myy:
                    return {{qtyCentering_[gridData_.iMyy][gridData_.idirX]}};
                case HybridQuantity::Scalar::Myz:
                    return {{qtyCentering_[gridData_.iMyz][gridData_.idirX]}};
                case HybridQuantity::Scalar::Mzz:
                    return {{qtyCentering_[gridData_.iMzz][gridData_.idirX]}};
                default: throw std::runtime_error("Wrong hybridQuantity");
            }
        }

        else if constexpr (dim == 2)
        {
            switch (hybridQuantity)
            {
                case HybridQuantity::Scalar::Bx:
                    return {{qtyCentering_[gridData_.iBx][gridData_.idirX],
                             qtyCentering_[gridData_.iBx][gridData_.idirY]}};
                case HybridQuantity::Scalar::By:
                    return {{qtyCentering_[gridData_.iBy][gridData_.idirX],
                             qtyCentering_[gridData_.iBy][gridData_.idirY]}};
                case HybridQuantity::Scalar::Bz:
                    return {{qtyCentering_[gridData_.iBz][gridData_.idirX],
                             qtyCentering_[gridData_.iBz][gridData_.idirY]}};
                case HybridQuantity::Scalar::Ex:
                    return {{qtyCentering_[gridData_.iEx][gridData_.idirX],
                             qtyCentering_[gridData_.iEx][gridData_.idirY]}};
                case HybridQuantity::Scalar::Ey:
                    return {{qtyCentering_[gridData_.iEy][gridData_.idirX],
                             qtyCentering_[gridData_.iEy][gridData_.idirY]}};
                case HybridQuantity::Scalar::Ez:
                    return {{qtyCentering_[gridData_.iEz][gridData_.idirX],
                             qtyCentering_[gridData_.iEz][gridData_.idirY]}};
                case HybridQuantity::Scalar::Jx:
                    return {{qtyCentering_[gridData_.iJx][gridData_.idirX],
                             qtyCentering_[gridData_.iJx][gridData_.idirY]}};
                case HybridQuantity::Scalar::Jy:
                    return {{qtyCentering_[gridData_.iJy][gridData_.idirX],
                             qtyCentering_[gridData_.iJy][gridData_.idirY]}};
                case HybridQuantity::Scalar::Jz:
                    return {{qtyCentering_[gridData_.iJz][gridData_.idirX],
                             qtyCentering_[gridData_.iJz][gridData_.idirY]}};
                case HybridQuantity::Scalar::rho:
                    return {{qtyCentering_[gridData_.irho][gridData_.idirX],
                             qtyCentering_[gridData_.irho][gridData_.idirY]}};
                case HybridQuantity::Scalar::Vx:
                    return {{qtyCentering_[gridData_.iVx][gridData_.idirX],
                             qtyCentering_[gridData_.iVx][gridData_.idirY]}};
                case HybridQuantity::Scalar::Vy:
                    return {{qtyCentering_[gridData_.iVy][gridData_.idirX],
                             qtyCentering_[gridData_.iVy][gridData_.idirY]}};
                case HybridQuantity::Scalar::Vz:
                    return {{qtyCentering_[gridData_.iVz][gridData_.idirX],
                             qtyCentering_[gridData_.iVz][gridData_.idirY]}};
                case HybridQuantity::Scalar::P:
                    return {{qtyCentering_[gridData_.iP][gridData_.idirX],
                             qtyCentering_[gridData_.iP][gridData_.idirY]}};
                case HybridQuantity::Scalar::Mxx:
                    return {{qtyCentering_[gridData_.iMxx][gridData_.idirX],
                             qtyCentering_[gridData_.iMxx][gridData_.idirY]}};
                case HybridQuantity::Scalar::Mxy:
                    return {{qtyCentering_[gridData_.iMxy][gridData_.idirX],
                             qtyCentering_[gridData_.iMxy][gridData_.idirY]}};
                case HybridQuantity::Scalar::Mxz:
                    return {{qtyCentering_[gridData_.iMxz][gridData_.idirX],
                             qtyCentering_[gridData_.iMxz][gridData_.idirY]}};
                case HybridQuantity::Scalar::Myy:
                    return {{qtyCentering_[gridData_.iMyy][gridData_.idirX],
                             qtyCentering_[gridData_.iMyy][gridData_.idirY]}};
                case HybridQuantity::Scalar::Myz:
                    return {{qtyCentering_[gridData_.iMyz][gridData_.idirX],
                             qtyCentering_[gridData_.iMyz][gridData_.idirY]}};
                case HybridQuantity::Scalar::Mzz:
                    return {{qtyCentering_[gridData_.iMzz][gridData_.idirX],
                             qtyCentering_[gridData_.iMzz][gridData_.idirY]}};
                default: throw std::runtime_error("Wrong hybridQuantity");
            }
        }

        else if constexpr (dim == 3)
        {
            switch (hybridQuantity)
            {
                case HybridQuantity::Scalar::Bx:
                    return {{qtyCentering_[gridData_.iBx][gridData_.idirX],
                             qtyCentering_[gridData_.iBx][gridData_.idirY],
                             qtyCentering_[gridData_.iBx][gridData_.idirZ]}};
                case HybridQuantity::Scalar::By:
                    return {{qtyCentering_[gridData_.iBy][gridData_.idirX],
                             qtyCentering_[gridData_.iBy][gridData_.idirY],
                             qtyCentering_[gridData_.iBy][gridData_.idirZ]}};
                case HybridQuantity::Scalar::Bz:
                    return {{qtyCentering_[gridData_.iBz][gridData_.idirX],
                             qtyCentering_[gridData_.iBz][gridData_.idirY],
                             qtyCentering_[gridData_.iBz][gridData_.idirZ]}};
                case HybridQuantity::Scalar::Ex:
                    return {{qtyCentering_[gridData_.iEx][gridData_.idirX],
                             qtyCentering_[gridData_.iEx][gridData_.idirY],
                             qtyCentering_[gridData_.iEx][gridData_.idirZ]}};
                case HybridQuantity::Scalar::Ey:
                    return {{qtyCentering_[gridData_.iEy][gridData_.idirX],
                             qtyCentering_[gridData_.iEy][gridData_.idirY],
                             qtyCentering_[gridData_.iEy][gridData_.idirZ]}};
                case HybridQuantity::Scalar::Ez:
                    return {{qtyCentering_[gridData_.iEz][gridData_.idirX],
                             qtyCentering_[gridData_.iEz][gridData_.idirY],
                             qtyCentering_[gridData_.iEz][gridData_.idirZ]}};
                case HybridQuantity::Scalar::Jx:
                    return {{qtyCentering_[gridData_.iJx][gridData_.idirX],
                             qtyCentering_[gridData_.iJx][gridData_.idirY],
                             qtyCentering_[gridData_.iJx][gridData_.idirZ]}};
                case HybridQuantity::Scalar::Jy:
                    return {{qtyCentering_[gridData_.iJy][gridData_.idirX],
                             qtyCentering_[gridData_.iJy][gridData_.idirY],
                             qtyCentering_[gridData_.iJy][gridData_.idirZ]}};
                case HybridQuantity::Scalar::Jz:
                    return {{qtyCentering_[gridData_.iJz][gridData_.idirX],
                             qtyCentering_[gridData_.iJz][gridData_.idirY],
                             qtyCentering_[gridData_.iJz][gridData_.idirZ]}};
                case HybridQuantity::Scalar::rho:
                    return {{qtyCentering_[gridData_.irho][gridData_.idirX],
                             qtyCentering_[gridData_.irho][gridData_.idirY],
                             qtyCentering_[gridData_.irho][gridData_.idirZ]}};
                case HybridQuantity::Scalar::Vx:
                    return {{qtyCentering_[gridData_.iVx][gridData_.idirX],
                             qtyCentering_[gridData_.iVx][gridData_.idirY],
                             qtyCentering_[gridData_.iVx][gridData_.idirZ]}};
                case HybridQuantity::Scalar::Vy:
                    return {{qtyCentering_[gridData_.iVy][gridData_.idirX],
                             qtyCentering_[gridData_.iVy][gridData_.idirY],
                             qtyCentering_[gridData_.iVy][gridData_.idirZ]}};
                case HybridQuantity::Scalar::Vz:
                    return {{qtyCentering_[gridData_.iVz][gridData_.idirX],
                             qtyCentering_[gridData_.iVz][gridData_.idirY],
                             qtyCentering_[gridData_.iVz][gridData_.idirZ]}};
                case HybridQuantity::Scalar::P:
                    return {{qtyCentering_[gridData_.iP][gridData_.idirX],
                             qtyCentering_[gridData_.iP][gridData_.idirY],
                             qtyCentering_[gridData_.iP][gridData_.idirZ]}};
                case HybridQuantity::Scalar::Mxx:
                    return {{qtyCentering_[gridData_.iMxx][gridData_.idirX],
                             qtyCentering_[gridData_.iMxx][gridData_.idirY],
                             qtyCentering_[gridData_.iMxx][gridData_.idirZ]}};
                case HybridQuantity::Scalar::Mxy:
                    return {{qtyCentering_[gridData_.iMxy][gridData_.idirX],
                             qtyCentering_[gridData_.iMxy][gridData_.idirY],
                             qtyCentering_[gridData_.iMxy][gridData_.idirZ]}};
                case HybridQuantity::Scalar::Mxz:
                    return {{qtyCentering_[gridData_.iMxz][gridData_.idirX],
                             qtyCentering_[gridData_.iMxz][gridData_.idirY],
                             qtyCentering_[gridData_.iMxz][gridData_.idirZ]}};
                case HybridQuantity::Scalar::Myy:
                    return {{qtyCentering_[gridData_.iMyy][gridData_.idirX],
                             qtyCentering_[gridData_.iMyy][gridData_.idirY],
                             qtyCentering_[gridData_.iMyy][gridData_.idirZ]}};
                case HybridQuantity::Scalar::Myz:
                    return {{qtyCentering_[gridData_.iMyz][gridData_.idirX],
                             qtyCentering_[gridData_.iMyz][gridData_.idirY],
                             qtyCentering_[gridData_.iMyz][gridData_.idirZ]}};
                case HybridQuantity::Scalar::Mzz:
                    return {{qtyCentering_[gridData_.iMzz][gridData_.idirX],
                             qtyCentering_[gridData_.iMzz][gridData_.idirY],
                             qtyCentering_[gridData_.iMzz][gridData_.idirZ]}};
                default: throw std::runtime_error("Wrong hybridQuantity");
            }
        }
    }



    NO_DISCARD constexpr static std::array<std::array<QtyCentering, dim>, 3>
    centering(HybridQuantity::Vector hybridQuantity)
    {
        switch (hybridQuantity)
        {
            case HybridQuantity::Vector::B:
                return {{centering(HybridQuantity::Scalar::Bx),
                         centering(HybridQuantity::Scalar::By),
                         centering(HybridQuantity::Scalar::Bz)}};

            case HybridQuantity::Vector::V:
                return {{centering(HybridQuantity::Scalar::Vx),
                         centering(HybridQuantity::Scalar::Vy),
                         centering(HybridQuantity::Scalar::Vz)}};

            case HybridQuantity::Vector::J:
                return {{centering(HybridQuantity::Scalar::Jx),
                         centering(HybridQuantity::Scalar::Jy),
                         centering(HybridQuantity::Scalar::Jz)}};

            case HybridQuantity::Vector::E:
                return {{centering(HybridQuantity::Scalar::Ex),
                         centering(HybridQuantity::Scalar::Ey),
                         centering(HybridQuantity::Scalar::Ez)}};


            default: throw std::runtime_error("Wrong hybridQuantity");
        }
    }


    enum class InterpDir { DualToPrimal = 0, PrimalToDual = 1 };

    // we might want to support the same interpolation possibilities as for the derivative, and
    // centralise the parametrisation of it
    template<auto dir, InterpDir interp_dir, std::size_t order = 2>
    NO_DISCARD static consteval auto directionalInterp()
    {
        if constexpr (dir >= dimension)
        {
            return std::array{WeightPoint{Point<int, dimension>{}, 1.0}};
        }
        else
        {
            constexpr int baseidx = (interp_dir == InterpDir::PrimalToDual);

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
    //

    NO_DISCARD auto static consteval momentsToEx()
    {
        // Ex is dual primal primal
        // moments are primal primal primal
        // operation is thus Ppp to Dpp
        // shift only in the X direction

        using PHARE::core::dirX;

        return directionalInterp<dirX, InterpDir::PrimalToDual>();
    }

    NO_DISCARD auto static consteval momentsToEy()
    {
        // Ey is primal dual primal
        // moments are primal primal primal
        // operation is thus pPp to pDp
        // shift only in the Y direction

        using PHARE::core::dirY;

        return directionalInterp<dirY, InterpDir::PrimalToDual>();
    }

    NO_DISCARD auto static consteval momentsToEz()
    {
        // Ez is primal  primal dual
        // moments are primal primal primal
        // operation is thus ppP to ppD
        // shift only in the Z direction

        using PHARE::core::dirZ;

        return directionalInterp<dirZ, InterpDir::PrimalToDual>();
    }



    NO_DISCARD auto static consteval BxToMoments()
    {
        // Bx is primal dual dual
        // moments are primal primal primal
        // operation is thus Pdd to Ppp

        using PHARE::core::dirY;
        using PHARE::core::dirZ;

        return tensorProduct<dirY, dirZ>(directionalInterp<dirY, InterpDir::DualToPrimal>(),
                                         directionalInterp<dirZ, InterpDir::DualToPrimal>());
    }



    NO_DISCARD auto static consteval ByToMoments()
    {
        // By is dual primal dual
        // moments are primal primal primal
        // operation is thus Dpd to Ppp

        using PHARE::core::dirX;
        using PHARE::core::dirZ;

        return tensorProduct<dirX, dirZ>(directionalInterp<dirX, InterpDir::DualToPrimal>(),
                                         directionalInterp<dirZ, InterpDir::DualToPrimal>());
    }




    NO_DISCARD auto static consteval BzToMoments()
    {
        // Bz is dual dual primal
        // moments are primal primal primal
        // operation is thus Ddp to Ppp

        using PHARE::core::dirX;
        using PHARE::core::dirY;

        return tensorProduct<dirX, dirY>(directionalInterp<dirX, InterpDir::DualToPrimal>(),
                                         directionalInterp<dirY, InterpDir::DualToPrimal>());
    }

    NO_DISCARD auto static consteval ExToMoments()
    {
        // Ex is dual primal primal
        // moments are primal primal primal
        // operation is thus Dpp to Ppp
        // shift only in the X direction

        using PHARE::core::dirX;

        return directionalInterp<dirX, InterpDir::DualToPrimal>();
    }

    NO_DISCARD auto static consteval EyToMoments()
    {
        // Ey is       primal dual   primal
        // moments are primal primal primal
        // operation is thus pDp to pPp
        // shift only in the Y direction

        using PHARE::core::dirY;

        return directionalInterp<dirY, InterpDir::DualToPrimal>();
    }

    NO_DISCARD auto static consteval EzToMoments()
    {
        // Ez is       primal primal dual
        // moments are primal primal primal
        // operation is thus ppD to ppP
        // shift only in the Z direction

        using PHARE::core::dirZ;

        return directionalInterp<dirZ, InterpDir::DualToPrimal>();
    }

    NO_DISCARD auto static consteval JxToMoments()
    {
        // Jx is dual primal primal
        // moments are primal primal primal
        // operation is thus Dpp to Ppp
        // shift only in the X direction

        using PHARE::core::dirX;

        return directionalInterp<dirX, InterpDir::DualToPrimal>();
    }

    NO_DISCARD auto static consteval JyToMoments()
    {
        // Jy is primal dual primal
        // moments are primal primal primal
        // operation is thus pDp to pPp
        // shift only in the Y direction

        using PHARE::core::dirY;

        return directionalInterp<dirY, InterpDir::DualToPrimal>();
    }

    NO_DISCARD auto static consteval JzToMoments()
    {
        // Jz is primal primal dual
        // moments are primal primal primal
        // operation is thus ppD to ppP
        // shift only in the Z direction

        using PHARE::core::dirZ;

        return directionalInterp<dirZ, InterpDir::DualToPrimal>();
    }

    NO_DISCARD auto static consteval BxToEx()
    {
        // Bx is primal dual dual
        // Ex is dual primal primal
        // operation is pdd to dpp

        using PHARE::core::dirX;
        using PHARE::core::dirY;
        using PHARE::core::dirZ;

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

    NO_DISCARD auto static consteval JxToEx()
    {
        // Jx is dual primal primal
        // Ex is dual primal primal
        // operation is thus dpp to dpp
        // no shift for a yee grid

        return std::array{WeightPoint{Point<int, dimension>{}, 1.0}};
    }

    NO_DISCARD auto static consteval JyToEy()
    {
        // Jy is primal dual primal
        // Ey is primal dual primal
        // operation is thus pdp to pdp
        // no shift for a yee grid

        return std::array{WeightPoint{Point<int, dimension>{}, 1.0}};
    }

    NO_DISCARD auto static consteval JzToEz()
    {
        // Jz is primal primal dual
        // Ez is primal primal dual
        // operation is thus ppd to ppd
        // no shift for a yee grid

        return std::array{WeightPoint{Point<int, dimension>{}, 1.0}};
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

    static constexpr std::uint32_t iBx = static_cast<std::uint32_t>(HybridQuantity::Scalar::Bx);
    static constexpr std::uint32_t iBy = static_cast<std::uint32_t>(HybridQuantity::Scalar::By);
    static constexpr std::uint32_t iBz = static_cast<std::uint32_t>(HybridQuantity::Scalar::Bz);

    static constexpr std::uint32_t iEx = static_cast<std::uint32_t>(HybridQuantity::Scalar::Ex);
    static constexpr std::uint32_t iEy = static_cast<std::uint32_t>(HybridQuantity::Scalar::Ey);
    static constexpr std::uint32_t iEz = static_cast<std::uint32_t>(HybridQuantity::Scalar::Ez);

    static constexpr std::uint32_t iJx = static_cast<std::uint32_t>(HybridQuantity::Scalar::Jx);
    static constexpr std::uint32_t iJy = static_cast<std::uint32_t>(HybridQuantity::Scalar::Jy);
    static constexpr std::uint32_t iJz = static_cast<std::uint32_t>(HybridQuantity::Scalar::Jz);

    static constexpr std::uint32_t irho = static_cast<std::uint32_t>(HybridQuantity::Scalar::rho);

    static constexpr std::uint32_t iVx = static_cast<std::uint32_t>(HybridQuantity::Scalar::Vx);
    static constexpr std::uint32_t iVy = static_cast<std::uint32_t>(HybridQuantity::Scalar::Vy);
    static constexpr std::uint32_t iVz = static_cast<std::uint32_t>(HybridQuantity::Scalar::Vz);

    static constexpr std::uint32_t iMxx = static_cast<std::uint32_t>(HybridQuantity::Scalar::Mxx);
    static constexpr std::uint32_t iMxy = static_cast<std::uint32_t>(HybridQuantity::Scalar::Mxy);
    static constexpr std::uint32_t iMxz = static_cast<std::uint32_t>(HybridQuantity::Scalar::Mxz);
    static constexpr std::uint32_t iMyy = static_cast<std::uint32_t>(HybridQuantity::Scalar::Myy);
    static constexpr std::uint32_t iMyz = static_cast<std::uint32_t>(HybridQuantity::Scalar::Myz);
    static constexpr std::uint32_t iMzz = static_cast<std::uint32_t>(HybridQuantity::Scalar::Mzz);

    static constexpr std::uint32_t iP = static_cast<std::uint32_t>(HybridQuantity::Scalar::P);
};



} // namespace PHARE::core::hybrid


#endif // PHARE_CORE_GRID_IMPL_YEE_GRIDLAYOUT_HYBRID_YEE_HPP
