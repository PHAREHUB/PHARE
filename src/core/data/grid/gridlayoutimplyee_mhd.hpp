#ifndef PHARE_CORE_GRID_GRIDLAYOUTYEE_MHD_HPP
#define PHARE_CORE_GRID_GRIDLAYOUTYEE_MHD_HPP

#include <array>
#include <cstddef>
#include <vector>

#include "core/def.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/utilities/ghost_width_calculator.hpp"
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
    template<std::size_t dim, std::size_t interpOrder, std::uint32_t reconstruction_nghosts_ = 3>
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
        // The MHD layout reserves ghosts based on the reconstruction stencil width,
        // plus extra layers for J Laplacian and hyper-resistivity corrections.
        static constexpr std::uint32_t reconstruction_nghosts = reconstruction_nghosts_;
        // Ghost width computed directly based on reconstruction stencil
        static constexpr std::uint32_t ghost_width
            = nbrGhostsFromReconstruction<reconstruction_nghosts>();

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
            gridDataT const data{};

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

            std::array<std::array<QtyCentering, NBR_COMPO>,
                       static_cast<std::size_t>(MHDQuantity::Scalar::count)> const _QtyCentering{
                Rho,
                Vx,
                Vy,
                Vz,
                Bx,
                By,
                Bz,
                P,
                rhoVx,
                rhoVy,
                rhoVz,
                Etot,
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

            return _QtyCentering;
        }

        //! says for each MHDQuantity::Quantity whether it is primal or dual, in each direction
        constexpr static std::array<std::array<QtyCentering, NBR_COMPO>,
                                    static_cast<std::size_t>(MHDQuantity::Scalar::count)> const
            _QtyCentering_{initLayoutCentering_()};

        static std::size_t const dim_{dim};

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
                    case MHDQuantity::Scalar::ScalarAllPrimal:
                        return {{_QtyCentering_[gridData_.iScalarAllPrimal][gridData_.idirX]}};
                    case MHDQuantity::Scalar::VecAllPrimalX:
                        return {{_QtyCentering_[gridData_.iVecAllPrimalX][gridData_.idirX]}};
                    case MHDQuantity::Scalar::VecAllPrimalY:
                        return {{_QtyCentering_[gridData_.iVecAllPrimalY][gridData_.idirX]}};
                    case MHDQuantity::Scalar::VecAllPrimalZ:
                        return {{_QtyCentering_[gridData_.iVecAllPrimalZ][gridData_.idirX]}};
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
                    case MHDQuantity::Scalar::ScalarAllPrimal:
                        return {{_QtyCentering_[gridData_.iScalarAllPrimal][gridData_.idirX],
                                 _QtyCentering_[gridData_.iScalarAllPrimal][gridData_.idirY]}};
                    case MHDQuantity::Scalar::VecAllPrimalX:
                        return {{_QtyCentering_[gridData_.iVecAllPrimalX][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecAllPrimalX][gridData_.idirY]}};
                    case MHDQuantity::Scalar::VecAllPrimalY:
                        return {{_QtyCentering_[gridData_.iVecAllPrimalY][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecAllPrimalY][gridData_.idirY]}};
                    case MHDQuantity::Scalar::VecAllPrimalZ:
                        return {{_QtyCentering_[gridData_.iVecAllPrimalZ][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecAllPrimalZ][gridData_.idirY]}};
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
                    case MHDQuantity::Scalar::ScalarAllPrimal:
                        return {{_QtyCentering_[gridData_.iScalarAllPrimal][gridData_.idirX],
                                 _QtyCentering_[gridData_.iScalarAllPrimal][gridData_.idirY],
                                 _QtyCentering_[gridData_.iScalarAllPrimal][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::VecAllPrimalX:
                        return {{_QtyCentering_[gridData_.iVecAllPrimalX][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecAllPrimalX][gridData_.idirY],
                                 _QtyCentering_[gridData_.iVecAllPrimalX][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::VecAllPrimalY:
                        return {{_QtyCentering_[gridData_.iVecAllPrimalY][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecAllPrimalY][gridData_.idirY],
                                 _QtyCentering_[gridData_.iVecAllPrimalY][gridData_.idirZ]}};
                    case MHDQuantity::Scalar::VecAllPrimalZ:
                        return {{_QtyCentering_[gridData_.iVecAllPrimalZ][gridData_.idirX],
                                 _QtyCentering_[gridData_.iVecAllPrimalZ][gridData_.idirY],
                                 _QtyCentering_[gridData_.iVecAllPrimalZ][gridData_.idirZ]}};
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

                case MHDQuantity::Vector::VecAllPrimal:
                    return {{centering(MHDQuantity::Scalar::VecAllPrimalX),
                             centering(MHDQuantity::Scalar::VecAllPrimalY),
                             centering(MHDQuantity::Scalar::VecAllPrimalZ)}};

                default: throw std::runtime_error("Wrong MHDQuantity");
            }
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
            return tensorProduct<dirX, dirY, dirZ>(
                directionalInterp<dirX, InterpDir::PrimalToDual>(),
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

            return tensorProduct<dirX, dirY, dirZ>(
                directionalInterp<dirX, InterpDir::DualToPrimal>(),
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

            return tensorProduct<dirX, dirY, dirZ>(
                directionalInterp<dirX, InterpDir::DualToPrimal>(),
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

        NO_DISCARD auto static consteval cellCenterToEdgeX()
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

        NO_DISCARD auto static consteval cellCenterToEdgeY()
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

        NO_DISCARD auto static consteval cellCenterToEdgeZ()
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

        NO_DISCARD auto static consteval faceXToCellCenter()
        {
            // The X face is Pdd
            // the mhd quantities in FV are Ddd
            // operation is thus Pdd to Ddd
            // shift only in the X direction

            using PHARE::core::dirX;

            return directionalInterp<dirX, InterpDir::PrimalToDual>();
        }

        NO_DISCARD auto static consteval faceYToCellCenter()
        {
            // The Y face is Dpd
            // the mhd quantities in FV are Ddd
            // operation is thus Dpd to Ddd
            // shift only in the Y direction

            using PHARE::core::dirY;

            return directionalInterp<dirY, InterpDir::PrimalToDual>();
        }

        NO_DISCARD auto static consteval faceZToCellCenter()
        {
            // The Z face is Ddp
            // the mhd quantities in FV are Ddd
            // operation is thus Ddp to Ddd
            // shift only in the Z direction

            using PHARE::core::dirZ;

            return directionalInterp<dirZ, InterpDir::PrimalToDual>();
        }

        NO_DISCARD auto static consteval edgeXToCellCenter()
        {
            // The X edge is dPP
            // the mhd quantities in FV are dDD
            // operation is thus dPP to dDD

            using PHARE::core::dirY;
            using PHARE::core::dirZ;

            return tensorProduct<dirY, dirZ>(directionalInterp<dirY, InterpDir::PrimalToDual>(),
                                             directionalInterp<dirZ, InterpDir::PrimalToDual>());
        }

        NO_DISCARD auto static consteval edgeYToCellCenter()
        {
            // The Y edge is PdP
            // the mhd quantities in FV are DdD
            // operation is thus PdP to DdD

            using PHARE::core::dirX;
            using PHARE::core::dirZ;

            return tensorProduct<dirX, dirZ>(directionalInterp<dirX, InterpDir::PrimalToDual>(),
                                             directionalInterp<dirZ, InterpDir::PrimalToDual>());
        }

        NO_DISCARD auto static consteval edgeZToCellCenter()
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
        NO_DISCARD auto static consteval cellCenterToFullPrimal()
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
} // namespace core
} // namespace PHARE

#endif // PHARE_CORE_GRID_GRIDLAYOUTYEE_MHD_HPP
