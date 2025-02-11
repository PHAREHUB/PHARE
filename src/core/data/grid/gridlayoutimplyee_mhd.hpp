#ifndef PHARE_CORE_GRID_GRIDLAYOUTYEE_MHD_HPP
#define PHARE_CORE_GRID_GRIDLAYOUTYEE_MHD_HPP

#include <array>
#include <cstddef>
#include <vector>

#include "core/def.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/types.hpp"
#include "gridlayoutdefs.hpp"

namespace PHARE
{
namespace core
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
    template<std::size_t dim, std::size_t interpOrder>
    class GridLayoutImplYeeMHD
    {
        // ------------------------------------------------------------------------
        //                              PRIVATE
        // ------------------------------------------------------------------------
    public:
        static constexpr std::size_t dimension    = dim;
        static constexpr std::size_t interp_order = interpOrder;
        static constexpr std::string_view type    = "yee";
        using quantity_type                       = MHDQuantity;

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
            const gridDataT data{};

            const std::array<QtyCentering, NBR_COMPO> Rho = {{data.dual, data.dual, data.dual}};

            const std::array<QtyCentering, NBR_COMPO> Vx = {{data.dual, data.dual, data.dual}};
            const std::array<QtyCentering, NBR_COMPO> Vy = {{data.dual, data.dual, data.dual}};
            const std::array<QtyCentering, NBR_COMPO> Vz = {{data.dual, data.dual, data.dual}};

            const std::array<QtyCentering, NBR_COMPO> Bx = {{data.primal, data.dual, data.dual}};
            const std::array<QtyCentering, NBR_COMPO> By = {{data.dual, data.primal, data.dual}};
            const std::array<QtyCentering, NBR_COMPO> Bz = {{data.dual, data.dual, data.primal}};

            const std::array<QtyCentering, NBR_COMPO> P = {{data.dual, data.dual, data.dual}};

            const std::array<QtyCentering, NBR_COMPO> rhoVx = {{data.dual, data.dual, data.dual}};
            const std::array<QtyCentering, NBR_COMPO> rhoVy = {{data.dual, data.dual, data.dual}};
            const std::array<QtyCentering, NBR_COMPO> rhoVz = {{data.dual, data.dual, data.dual}};

            const std::array<QtyCentering, NBR_COMPO> Etot = {{data.dual, data.dual, data.dual}};

            const std::array<QtyCentering, NBR_COMPO> Ex = {{data.dual, data.primal, data.primal}};
            const std::array<QtyCentering, NBR_COMPO> Ey = {{data.primal, data.dual, data.primal}};
            const std::array<QtyCentering, NBR_COMPO> Ez = {{data.primal, data.primal, data.dual}};

            const std::array<QtyCentering, NBR_COMPO> Jx = {{data.dual, data.primal, data.primal}};
            const std::array<QtyCentering, NBR_COMPO> Jy = {{data.primal, data.dual, data.primal}};
            const std::array<QtyCentering, NBR_COMPO> Jz = {{data.primal, data.primal, data.dual}};


            const std::array<QtyCentering, NBR_COMPO> ScalarFlux_x
                = {{data.primal, data.dual, data.dual}};

            const std::array<QtyCentering, NBR_COMPO> ScalarFlux_y
                = {{data.dual, data.primal, data.dual}};

            const std::array<QtyCentering, NBR_COMPO> ScalarFlux_z
                = {{data.dual, data.dual, data.primal}};

            const std::array<QtyCentering, NBR_COMPO> VecFluxX_x
                = {{data.primal, data.dual, data.dual}};
            const std::array<QtyCentering, NBR_COMPO> VecFluxY_x
                = {{data.primal, data.dual, data.dual}};
            const std::array<QtyCentering, NBR_COMPO> VecFluxZ_x
                = {{data.primal, data.dual, data.dual}};

            const std::array<QtyCentering, NBR_COMPO> VecFluxX_y
                = {{data.dual, data.primal, data.dual}};
            const std::array<QtyCentering, NBR_COMPO> VecFluxY_y
                = {{data.dual, data.primal, data.dual}};
            const std::array<QtyCentering, NBR_COMPO> VecFluxZ_y
                = {{data.dual, data.primal, data.dual}};

            const std::array<QtyCentering, NBR_COMPO> VecFluxX_z
                = {{data.dual, data.dual, data.primal}};
            const std::array<QtyCentering, NBR_COMPO> VecFluxY_z
                = {{data.dual, data.dual, data.primal}};
            const std::array<QtyCentering, NBR_COMPO> VecFluxZ_z
                = {{data.dual, data.dual, data.primal}};


            const std::array<std::array<QtyCentering, NBR_COMPO>,
                             static_cast<std::size_t>(MHDQuantity::Scalar::count)>
                _QtyCentering{Rho,          Vx,         Vy,         Vz,           Bx,
                              By,           Bz,         P,          rhoVx,        rhoVy,
                              rhoVz,        Etot,       Ex,         Ey,           Ez,
                              Jx,           Jy,         Jz,         ScalarFlux_x, ScalarFlux_y,
                              ScalarFlux_z, VecFluxX_x, VecFluxY_x, VecFluxZ_x,   VecFluxX_y,
                              VecFluxY_y,   VecFluxZ_y, VecFluxX_z, VecFluxY_z,   VecFluxZ_z};

            return _QtyCentering;
        }

        //! says for each MHDQuantity::Quantity whether it is primal or dual, in each direction
        constexpr const static std::array<std::array<QtyCentering, NBR_COMPO>,
                                          static_cast<std::size_t>(MHDQuantity::Scalar::count)>
            _QtyCentering_{initLayoutCentering_()};

        static const std::size_t dim_{dim};

        // ------------------------------------------------------------------------
        //                          PUBLIC INTERFACE
        // ------------------------------------------------------------------------
    public:
        NO_DISCARD constexpr static std::array<QtyCentering, dim>
        centering(MHDQuantity::Scalar MHDQuantity)
        {
            constexpr gridDataT_mhd gridData_{};
            if constexpr (dim == 1)
            {
                switch (MHDQuantity)
                {
                    case MHDQuantity::Scalar::rho:
                        return {{_QtyCentering_[gridData_.irho][gridData_.idirX]}};
                    case MHDQuantity::Scalar::Vx:
                        return {{_QtyCentering_[gridData_.iVx][gridData_.idirX]}};
                    case MHDQuantity::Scalar::Vy:
                        return {{_QtyCentering_[gridData_.iVy][gridData_.idirX]}};
                    case MHDQuantity::Scalar::Vz:
                        return {{_QtyCentering_[gridData_.iVz][gridData_.idirX]}};
                    case MHDQuantity::Scalar::Bx:
                        return {{_QtyCentering_[gridData_.iBx][gridData_.idirX]}};
                    case MHDQuantity::Scalar::By:
                        return {{_QtyCentering_[gridData_.iBy][gridData_.idirX]}};
                    case MHDQuantity::Scalar::Bz:
                        return {{_QtyCentering_[gridData_.iBz][gridData_.idirX]}};
                    case MHDQuantity::Scalar::P:
                        return {{_QtyCentering_[gridData_.iP][gridData_.idirX]}};
                    case MHDQuantity::Scalar::rhoVx:
                        return {{_QtyCentering_[gridData_.irhoVx][gridData_.idirX]}};
                    case MHDQuantity::Scalar::rhoVy:
                        return {{_QtyCentering_[gridData_.irhoVy][gridData_.idirX]}};
                    case MHDQuantity::Scalar::rhoVz:
                        return {{_QtyCentering_[gridData_.irhoVz][gridData_.idirX]}};
                    case MHDQuantity::Scalar::Etot:
                        return {{_QtyCentering_[gridData_.iEtot][gridData_.idirX]}};
                    case MHDQuantity::Scalar::Ex:
                        return {{_QtyCentering_[gridData_.iEx][gridData_.idirX]}};
                    case MHDQuantity::Scalar::Ey:
                        return {{_QtyCentering_[gridData_.iEy][gridData_.idirX]}};
                    case MHDQuantity::Scalar::Ez:
                        return {{_QtyCentering_[gridData_.iEz][gridData_.idirX]}};
                    case MHDQuantity::Scalar::Jx:
                        return {{_QtyCentering_[gridData_.iJx][gridData_.idirX]}};
                    case MHDQuantity::Scalar::Jy:
                        return {{_QtyCentering_[gridData_.iJy][gridData_.idirX]}};
                    case MHDQuantity::Scalar::Jz:
                        return {{_QtyCentering_[gridData_.iJz][gridData_.idirX]}};
                    case MHDQuantity::Scalar::ScalarFlux_x:
                        return {{_QtyCentering_[gridData_.iScalarFlux_x][gridData_.idirX]}};
                    case MHDQuantity::Scalar::VecFluxX_x:
                        return {{_QtyCentering_[gridData_.iVecFluxX_x][gridData_.idirX]}};
                    case MHDQuantity::Scalar::VecFluxY_x:
                        return {{_QtyCentering_[gridData_.iVecFluxY_x][gridData_.idirX]}};
                    case MHDQuantity::Scalar::VecFluxZ_x:
                        return {{_QtyCentering_[gridData_.iVecFluxZ_x][gridData_.idirX]}};
                    default: throw std::runtime_error("Wrong MHDQuantity");
                }
            }

            else if constexpr (dim == 2)
            {
                switch (MHDQuantity)
                {
                    case MHDQuantity::Scalar::rho:
                        return {{_QtyCentering_[gridData_.irho][gridData_.idirX],
                                 _QtyCentering_[gridData_.irho][gridData_.idirY]}};
                    case MHDQuantity::Scalar::Vx:
                        return {{_QtyCentering_[gridData_.iVx][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVx][gridData_.idirY]}};
                    case MHDQuantity::Scalar::Vy:
                        return {{_QtyCentering_[gridData_.iVy][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVy][gridData_.idirY]}};
                    case MHDQuantity::Scalar::Vz:
                        return {{_QtyCentering_[gridData_.iVz][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVz][gridData_.idirY]}};
                    case MHDQuantity::Scalar::Bx:
                        return {{_QtyCentering_[gridData_.iBx][gridData_.idirX],
                                 _QtyCentering_[gridData_.iBx][gridData_.idirY]}};
                    case MHDQuantity::Scalar::By:
                        return {{_QtyCentering_[gridData_.iBy][gridData_.idirX],
                                 _QtyCentering_[gridData_.iBy][gridData_.idirY]}};
                    case MHDQuantity::Scalar::Bz:
                        return {{_QtyCentering_[gridData_.iBz][gridData_.idirX],
                                 _QtyCentering_[gridData_.iBz][gridData_.idirY]}};
                    case MHDQuantity::Scalar::P:
                        return {{_QtyCentering_[gridData_.iP][gridData_.idirX],
                                 _QtyCentering_[gridData_.iP][gridData_.idirY]}};
                    case MHDQuantity::Scalar::rhoVx:
                        return {{_QtyCentering_[gridData_.irhoVx][gridData_.idirX],
                                 _QtyCentering_[gridData_.irhoVx][gridData_.idirY]}};
                    case MHDQuantity::Scalar::rhoVy:
                        return {{_QtyCentering_[gridData_.irhoVy][gridData_.idirX],
                                 _QtyCentering_[gridData_.irhoVy][gridData_.idirY]}};
                    case MHDQuantity::Scalar::rhoVz:
                        return {{_QtyCentering_[gridData_.irhoVz][gridData_.idirX],
                                 _QtyCentering_[gridData_.irhoVz][gridData_.idirY]}};
                    case MHDQuantity::Scalar::Etot:
                        return {{_QtyCentering_[gridData_.iEtot][gridData_.idirX],
                                 _QtyCentering_[gridData_.iEtot][gridData_.idirY]}};
                    case MHDQuantity::Scalar::Ex:
                        return {{_QtyCentering_[gridData_.iEx][gridData_.idirX],
                                 _QtyCentering_[gridData_.iEx][gridData_.idirY]}};
                    case MHDQuantity::Scalar::Ey:
                        return {{_QtyCentering_[gridData_.iEy][gridData_.idirX],
                                 _QtyCentering_[gridData_.iEy][gridData_.idirY]}};
                    case MHDQuantity::Scalar::Ez:
                        return {{_QtyCentering_[gridData_.iEz][gridData_.idirX],
                                 _QtyCentering_[gridData_.iEz][gridData_.idirY]}};
                    case MHDQuantity::Scalar::Jx:
                        return {{_QtyCentering_[gridData_.iJx][gridData_.idirX],
                                 _QtyCentering_[gridData_.iJx][gridData_.idirY]}};
                    case MHDQuantity::Scalar::Jy:
                        return {{_QtyCentering_[gridData_.iJy][gridData_.idirX],
                                 _QtyCentering_[gridData_.iJy][gridData_.idirY]}};
                    case MHDQuantity::Scalar::Jz:
                        return {{_QtyCentering_[gridData_.iJz][gridData_.idirX],
                                 _QtyCentering_[gridData_.iJz][gridData_.idirY]}};
                    case MHDQuantity::Scalar::ScalarFlux_x:
                        return {{_QtyCentering_[gridData_.iScalarFlux_x][gridData_.idirX],
                                 _QtyCentering_[gridData_.iScalarFlux_x][gridData_.idirY]}};
                    case MHDQuantity::Scalar::ScalarFlux_y:
                        return {{_QtyCentering_[gridData_.iScalarFlux_y][gridData_.idirX],
                                 _QtyCentering_[gridData_.iScalarFlux_y][gridData_.idirY]}};
                    case MHDQuantity::Scalar::VecFluxX_x:
                        return {{_QtyCentering_[gridData_.iVecFluxX_x][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecFluxX_x][gridData_.idirY]}};
                    case MHDQuantity::Scalar::VecFluxY_x:
                        return {{_QtyCentering_[gridData_.iVecFluxY_x][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecFluxY_x][gridData_.idirY]}};
                    case MHDQuantity::Scalar::VecFluxZ_x:
                        return {{_QtyCentering_[gridData_.iVecFluxZ_x][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecFluxZ_x][gridData_.idirY]}};
                    case MHDQuantity::Scalar::VecFluxX_y:
                        return {{_QtyCentering_[gridData_.iVecFluxX_y][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecFluxX_y][gridData_.idirY]}};
                    case MHDQuantity::Scalar::VecFluxY_y:
                        return {{_QtyCentering_[gridData_.iVecFluxY_y][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecFluxY_y][gridData_.idirY]}};
                    case MHDQuantity::Scalar::VecFluxZ_y:
                        return {{_QtyCentering_[gridData_.iVecFluxZ_y][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecFluxZ_y][gridData_.idirY]}};
                    default: throw std::runtime_error("Wrong MHDQuantity");
                }
            }

            else if constexpr (dim == 3)
            {
                switch (MHDQuantity)
                {
                    case MHDQuantity::Scalar::rho:
                        return {{_QtyCentering_[gridData_.irho][gridData_.idirX],
                                 _QtyCentering_[gridData_.irho][gridData_.idirY],
                                 _QtyCentering_[gridData_.irho][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::Vx:
                        return {{_QtyCentering_[gridData_.iVx][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVx][gridData_.idirY],
                                 _QtyCentering_[gridData_.iVx][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::Vy:
                        return {{_QtyCentering_[gridData_.iVy][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVy][gridData_.idirY],
                                 _QtyCentering_[gridData_.iVy][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::Vz:
                        return {{_QtyCentering_[gridData_.iVz][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVz][gridData_.idirY],
                                 _QtyCentering_[gridData_.iVz][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::Bx:
                        return {{_QtyCentering_[gridData_.iBx][gridData_.idirX],
                                 _QtyCentering_[gridData_.iBx][gridData_.idirY],
                                 _QtyCentering_[gridData_.iBx][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::By:
                        return {{_QtyCentering_[gridData_.iBy][gridData_.idirX],
                                 _QtyCentering_[gridData_.iBy][gridData_.idirY],
                                 _QtyCentering_[gridData_.iBy][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::Bz:
                        return {{_QtyCentering_[gridData_.iBz][gridData_.idirX],
                                 _QtyCentering_[gridData_.iBz][gridData_.idirY],
                                 _QtyCentering_[gridData_.iBz][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::P:
                        return {{_QtyCentering_[gridData_.iP][gridData_.idirX],
                                 _QtyCentering_[gridData_.iP][gridData_.idirY],
                                 _QtyCentering_[gridData_.iP][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::rhoVx:
                        return {{_QtyCentering_[gridData_.irhoVx][gridData_.idirX],
                                 _QtyCentering_[gridData_.irhoVx][gridData_.idirY],
                                 _QtyCentering_[gridData_.irhoVx][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::rhoVy:
                        return {{_QtyCentering_[gridData_.irhoVy][gridData_.idirX],
                                 _QtyCentering_[gridData_.irhoVy][gridData_.idirY],
                                 _QtyCentering_[gridData_.irhoVy][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::rhoVz:
                        return {{_QtyCentering_[gridData_.irhoVz][gridData_.idirX],
                                 _QtyCentering_[gridData_.irhoVz][gridData_.idirY],
                                 _QtyCentering_[gridData_.irhoVz][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::Etot:
                        return {{_QtyCentering_[gridData_.iEtot][gridData_.idirX],
                                 _QtyCentering_[gridData_.iEtot][gridData_.idirY],
                                 _QtyCentering_[gridData_.iEtot][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::Ex:
                        return {{_QtyCentering_[gridData_.iEx][gridData_.idirX],
                                 _QtyCentering_[gridData_.iEx][gridData_.idirY],
                                 _QtyCentering_[gridData_.iEx][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::Ey:
                        return {{_QtyCentering_[gridData_.iEy][gridData_.idirX],
                                 _QtyCentering_[gridData_.iEy][gridData_.idirY],
                                 _QtyCentering_[gridData_.iEy][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::Ez:
                        return {{_QtyCentering_[gridData_.iEz][gridData_.idirX],
                                 _QtyCentering_[gridData_.iEz][gridData_.idirY],
                                 _QtyCentering_[gridData_.iEz][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::Jx:
                        return {{_QtyCentering_[gridData_.iJx][gridData_.idirX],
                                 _QtyCentering_[gridData_.iJx][gridData_.idirY],
                                 _QtyCentering_[gridData_.iJx][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::Jy:
                        return {{_QtyCentering_[gridData_.iJy][gridData_.idirX],
                                 _QtyCentering_[gridData_.iJy][gridData_.idirY],
                                 _QtyCentering_[gridData_.iJy][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::Jz:
                        return {{_QtyCentering_[gridData_.iJz][gridData_.idirX],
                                 _QtyCentering_[gridData_.iJz][gridData_.idirY],
                                 _QtyCentering_[gridData_.iJz][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::ScalarFlux_x:
                        return {{_QtyCentering_[gridData_.iScalarFlux_x][gridData_.idirX],
                                 _QtyCentering_[gridData_.iScalarFlux_x][gridData_.idirY],
                                 _QtyCentering_[gridData_.iScalarFlux_x][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::ScalarFlux_y:
                        return {{_QtyCentering_[gridData_.iScalarFlux_y][gridData_.idirX],
                                 _QtyCentering_[gridData_.iScalarFlux_y][gridData_.idirY],
                                 _QtyCentering_[gridData_.iScalarFlux_y][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::ScalarFlux_z:
                        return {{_QtyCentering_[gridData_.iScalarFlux_z][gridData_.idirX],
                                 _QtyCentering_[gridData_.iScalarFlux_z][gridData_.idirY],
                                 _QtyCentering_[gridData_.iScalarFlux_z][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::VecFluxX_x:
                        return {{_QtyCentering_[gridData_.iVecFluxX_x][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecFluxX_x][gridData_.idirY],
                                 _QtyCentering_[gridData_.iVecFluxX_x][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::VecFluxY_x:
                        return {{_QtyCentering_[gridData_.iVecFluxY_x][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecFluxY_x][gridData_.idirY],
                                 _QtyCentering_[gridData_.iVecFluxY_x][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::VecFluxZ_x:
                        return {{_QtyCentering_[gridData_.iVecFluxZ_x][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecFluxZ_x][gridData_.idirY],
                                 _QtyCentering_[gridData_.iVecFluxZ_x][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::VecFluxX_y:
                        return {{_QtyCentering_[gridData_.iVecFluxX_y][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecFluxX_y][gridData_.idirY],
                                 _QtyCentering_[gridData_.iVecFluxX_y][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::VecFluxY_y:
                        return {{_QtyCentering_[gridData_.iVecFluxY_y][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecFluxY_y][gridData_.idirY],
                                 _QtyCentering_[gridData_.iVecFluxY_y][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::VecFluxZ_y:
                        return {{_QtyCentering_[gridData_.iVecFluxZ_y][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecFluxZ_y][gridData_.idirY],
                                 _QtyCentering_[gridData_.iVecFluxZ_y][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::VecFluxX_z:
                        return {{_QtyCentering_[gridData_.iVecFluxX_z][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecFluxX_z][gridData_.idirY],
                                 _QtyCentering_[gridData_.iVecFluxX_z][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::VecFluxY_z:
                        return {{_QtyCentering_[gridData_.iVecFluxY_z][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecFluxY_z][gridData_.idirY],
                                 _QtyCentering_[gridData_.iVecFluxY_z][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::VecFluxZ_z:
                        return {{_QtyCentering_[gridData_.iVecFluxZ_z][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecFluxZ_z][gridData_.idirY],
                                 _QtyCentering_[gridData_.iVecFluxZ_z][gridData_.idirZ]}};
                    default: throw std::runtime_error("Wrong MHDQuantity");
                }
            }
        }

        NO_DISCARD constexpr static std::array<std::array<QtyCentering, dim>, 3>
        centering(MHDQuantity::Vector MHDQuantity)
        {
            switch (MHDQuantity)
            {
                case MHDQuantity::Vector::V:
                    return {{centering(MHDQuantity::Scalar::Vx), centering(MHDQuantity::Scalar::Vy),
                             centering(MHDQuantity::Scalar::Vz)}};

                case MHDQuantity::Vector::B:
                    return {{centering(MHDQuantity::Scalar::Bx), centering(MHDQuantity::Scalar::By),
                             centering(MHDQuantity::Scalar::Bz)}};

                case MHDQuantity::Vector::rhoV:
                    return {{centering(MHDQuantity::Scalar::rhoVx),
                             centering(MHDQuantity::Scalar::rhoVy),
                             centering(MHDQuantity::Scalar::rhoVz)}};

                case MHDQuantity::Vector::E:
                    return {{centering(MHDQuantity::Scalar::Ex), centering(MHDQuantity::Scalar::Ey),
                             centering(MHDQuantity::Scalar::Ez)}};

                case MHDQuantity::Vector::J:
                    return {{centering(MHDQuantity::Scalar::Jx), centering(MHDQuantity::Scalar::Jy),
                             centering(MHDQuantity::Scalar::Jz)}};

                case MHDQuantity::Vector::VecFlux_x:
                    return {{centering(MHDQuantity::Scalar::VecFluxX_x),
                             centering(MHDQuantity::Scalar::VecFluxY_x),
                             centering(MHDQuantity::Scalar::VecFluxZ_x)}};

                case MHDQuantity::Vector::VecFlux_y:
                    return {{centering(MHDQuantity::Scalar::VecFluxX_y),
                             centering(MHDQuantity::Scalar::VecFluxY_y),
                             centering(MHDQuantity::Scalar::VecFluxZ_y)}};

                case MHDQuantity::Vector::VecFlux_z:
                    return {{centering(MHDQuantity::Scalar::VecFluxX_z),
                             centering(MHDQuantity::Scalar::VecFluxY_z),
                             centering(MHDQuantity::Scalar::VecFluxZ_z)}};

                default: throw std::runtime_error("Wrong MHDQuantity");
            }
        }

        NO_DISCARD auto static constexpr dualToPrimal() { return -1; }

        NO_DISCARD auto static constexpr primalToDual() { return 1; }

        NO_DISCARD auto static constexpr faceXToCellCenter()
        {
            // The X face is Pdd
            // the mhd quantities in FV are Ddd
            // operation is thus Pdd to Ddd
            // shift only in the X direction

            auto constexpr iShift = primalToDual();

            if constexpr (dimension == 1)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{iShift}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
            else if constexpr (dimension == 2)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{iShift, 0}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
            else if constexpr (dimension == 3)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0, 0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{iShift, 0, 0}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
        }

        NO_DISCARD auto static constexpr faceYToCellCenter()
        {
            // The Y face is Dpd
            // the mhd quantities in FV are Ddd
            // operation is thus Dpd to Ddd
            // shift only in the Y direction

            [[maybe_unused]] auto constexpr iShift = primalToDual();

            if constexpr (dimension == 1)
            {
                // since the linear combination is in the Y direction
                // in 1D the quantities are already on the Y face so return 1 point with no shift
                // with coef 1.
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 1.};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            else if constexpr (dimension == 2)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{0, iShift}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
            else if constexpr (dimension == 3)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0, 0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{0, iShift, 0}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
        }

        NO_DISCARD auto static constexpr faceZToCellCenter()
        {
            // The Z face is Ddp
            // the mhd quantities in FV are Ddd
            // operation is thus Ddp to Ddd
            // shift only in the Z direction

            [[maybe_unused]] auto constexpr iShift = primalToDual();

            if constexpr (dimension == 1)
            {
                // since the linear combination is in the Z direction
                // in 1D or 2D the quantities are already on the Z face so return 1 point with
                // no
                // shift with coef 1.
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 1.};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            else if constexpr (dimension == 2)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0}, 1.};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            else if constexpr (dimension == 3)
            {
                // in 3D we need two points, the second with a primalToDual shift along Z
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0, 0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{0, 0, iShift}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
        }

        NO_DISCARD auto static constexpr edgeXToCellCenter()
        {
            // The X face is Pdd
            // the mhd quantities in FV are Ddd
            // operation is thus Pdd to Ddd
            // shift only in the X direction

            auto constexpr iShift = primalToDual();

            if constexpr (dimension == 1)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 1.};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            else if constexpr (dimension == 2)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{0, iShift}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
            else if constexpr (dimension == 3)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0, 0}, 0.25};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{0, iShift, 0}, 0.25};
                constexpr WeightPoint<dimension> P3{Point<int, dimension>{0, 0, iShift}, 0.25};
                constexpr WeightPoint<dimension> P4{Point<int, dimension>{0, iShift, iShift}, 0.25};
                return std::array<WeightPoint<dimension>, 4>{P1, P2, P3, P4};
            }
        }

        NO_DISCARD auto static constexpr edgeYToCellCenter()
        {
            // The Y face is Dpd
            // the mhd quantities in FV are Ddd
            // operation is thus Dpd to Ddd
            // shift only in the Y direction

            [[maybe_unused]] auto constexpr iShift = primalToDual();

            if constexpr (dimension == 1)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{iShift}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
            else if constexpr (dimension == 2)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{iShift, 0}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
            else if constexpr (dimension == 3)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0, 0}, 0.25};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{iShift, 0, 0}, 0.25};
                constexpr WeightPoint<dimension> P3{Point<int, dimension>{0, 0, iShift}, 0.25};
                constexpr WeightPoint<dimension> P4{Point<int, dimension>{iShift, 0, iShift}, 0.25};
                return std::array<WeightPoint<dimension>, 4>{P1, P2, P3, P4};
            }
        }

        NO_DISCARD auto static constexpr edgeZToCellCenter()
        {
            // The Z face is Ddp
            // the mhd quantities in FV are Ddd
            // operation is thus Ddp to Ddd
            // shift only in the Z direction

            [[maybe_unused]] auto constexpr iShift = primalToDual();

            if constexpr (dimension == 1)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{iShift}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
            else if constexpr (dimension == 2)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0}, 0.25};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{iShift, 0}, 0.25};
                constexpr WeightPoint<dimension> P3{Point<int, dimension>{0, iShift}, 0.25};
                constexpr WeightPoint<dimension> P4{Point<int, dimension>{iShift, iShift}, 0.25};
                return std::array<WeightPoint<dimension>, 4>{P1, P2, P3, P4};
            }
            else if constexpr (dimension == 3)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0, 0}, 0.25};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{iShift, 0, 0}, 0.25};
                constexpr WeightPoint<dimension> P3{Point<int, dimension>{0, iShift, 0}, 0.25};
                constexpr WeightPoint<dimension> P4{Point<int, dimension>{iShift, iShift, 0}, 0.25};
                return std::array<WeightPoint<dimension>, 4>{P1, P2, P3, P4};
            }
        }

    }; // namespace core

} // namespace core
} // namespace PHARE

#endif // PHARE_CORE_GRID_GRIDLAYOUTYEE_MHD_HPP
