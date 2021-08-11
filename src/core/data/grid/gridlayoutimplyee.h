#ifndef PHARE_CORE_GRID_GRIDLAYOUTYEE_H
#define PHARE_CORE_GRID_GRIDLAYOUTYEE_H


#include "core/def.h"
#include "core/hybrid/hybrid_quantities.h"
#include "core/utilities/types.h"
#include "gridlayoutdefs.h"
#include "core/utilities/constants.h"

#include <array>
#include <vector>

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
    class GridLayoutImplYee
    {
        // ------------------------------------------------------------------------
        //                              PRIVATE
        // ------------------------------------------------------------------------
    public:
        static constexpr std::size_t dimension    = dim;
        static constexpr std::size_t interp_order = interpOrder;
        static constexpr std::string_view type    = "yee";

        /*
        void constexpr initLinearCombinations_();

        LinearCombination momentsToEx_;
        LinearCombination momentsToEy_;
        LinearCombination momentsToEz_;
        LinearCombination BxToEy_;
        LinearCombination BxToEz_;
        LinearCombination ByToEx_;
        LinearCombination ByToEz_;
        LinearCombination BzToEx_;
        LinearCombination BzToEy_;
        LinearCombination ExToMoment_;
        LinearCombination EyToMoment_;
        LinearCombination EzToMoment_;
        */

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
            const gridDataT data{};
            const std::array<QtyCentering, NBR_COMPO> Bx = {{data.primal, data.dual, data.dual}};
            const std::array<QtyCentering, NBR_COMPO> By = {{data.dual, data.primal, data.dual}};
            const std::array<QtyCentering, NBR_COMPO> Bz = {{data.dual, data.dual, data.primal}};

            const std::array<QtyCentering, NBR_COMPO> Ex = {{data.dual, data.primal, data.primal}};
            const std::array<QtyCentering, NBR_COMPO> Ey = {{data.primal, data.dual, data.primal}};
            const std::array<QtyCentering, NBR_COMPO> Ez = {{data.primal, data.primal, data.dual}};


            const std::array<QtyCentering, NBR_COMPO> Jx = {{data.dual, data.primal, data.primal}};
            const std::array<QtyCentering, NBR_COMPO> Jy = {{data.primal, data.dual, data.primal}};
            const std::array<QtyCentering, NBR_COMPO> Jz = {{data.primal, data.primal, data.dual}};

            const std::array<QtyCentering, NBR_COMPO> Rho
                = {{data.primal, data.primal, data.primal}};

            const std::array<QtyCentering, NBR_COMPO> Vx
                = {{data.primal, data.primal, data.primal}};
            const std::array<QtyCentering, NBR_COMPO> Vy
                = {{data.primal, data.primal, data.primal}};
            const std::array<QtyCentering, NBR_COMPO> Vz
                = {{data.primal, data.primal, data.primal}};

            const std::array<QtyCentering, NBR_COMPO> P = {{data.primal, data.primal, data.primal}};

            const std::array<std::array<QtyCentering, NBR_COMPO>,
                             static_cast<std::size_t>(HybridQuantity::Scalar::count)>
                hybridQtyCentering{Bx, By, Bz, Ex, Ey, Ez, Jx, Jy, Jz, Rho, Vx, Vy, Vz, P};


            return hybridQtyCentering;
        }



        //! says for each HybridQuantity::Quantity whether it is primal or dual, in each direction
        constexpr const static std::array<std::array<QtyCentering, NBR_COMPO>,
                                          static_cast<std::size_t>(HybridQuantity::Scalar::count)>
            hybridQtyCentering_{initLayoutCentering_()};

        static const std::size_t dim_{dim};

        // ------------------------------------------------------------------------
        //                          PUBLIC INTERFACE
        // ------------------------------------------------------------------------
    public:
        constexpr static std::array<QtyCentering, dim>
        centering(HybridQuantity::Scalar hybridQuantity)
        {
            constexpr gridDataT gridData_{};
            if constexpr (dim == 1)
            {
                switch (hybridQuantity)
                {
                    case HybridQuantity::Scalar::Bx:
                        return {{hybridQtyCentering_[gridData_.iBx][gridData_.idirX]}};
                    case HybridQuantity::Scalar::By:
                        return {{hybridQtyCentering_[gridData_.iBy][gridData_.idirX]}};
                    case HybridQuantity::Scalar::Bz:
                        return {{hybridQtyCentering_[gridData_.iBz][gridData_.idirX]}};
                    case HybridQuantity::Scalar::Ex:
                        return {{hybridQtyCentering_[gridData_.iEx][gridData_.idirX]}};
                    case HybridQuantity::Scalar::Ey:
                        return {{hybridQtyCentering_[gridData_.iEy][gridData_.idirX]}};
                    case HybridQuantity::Scalar::Ez:
                        return {{hybridQtyCentering_[gridData_.iEz][gridData_.idirX]}};
                    case HybridQuantity::Scalar::Jx:
                        return {{hybridQtyCentering_[gridData_.iJx][gridData_.idirX]}};
                    case HybridQuantity::Scalar::Jy:
                        return {{hybridQtyCentering_[gridData_.iJy][gridData_.idirX]}};
                    case HybridQuantity::Scalar::Jz:
                        return {{hybridQtyCentering_[gridData_.iJz][gridData_.idirX]}};
                    case HybridQuantity::Scalar::rho:
                        return {{hybridQtyCentering_[gridData_.irho][gridData_.idirX]}};
                    case HybridQuantity::Scalar::Vx:
                        return {{hybridQtyCentering_[gridData_.iVx][gridData_.idirX]}};
                    case HybridQuantity::Scalar::Vy:
                        return {{hybridQtyCentering_[gridData_.iVy][gridData_.idirX]}};
                    case HybridQuantity::Scalar::Vz:
                        return {{hybridQtyCentering_[gridData_.iVz][gridData_.idirX]}};
                    case HybridQuantity::Scalar::P:
                        return {{hybridQtyCentering_[gridData_.iP][gridData_.idirX]}};
                    default: throw_runtime_error("Wrong hybridQuantity");
                }
            }

            else if constexpr (dim == 2)
            {
                switch (hybridQuantity)
                {
                    case HybridQuantity::Scalar::Bx:
                        return {{hybridQtyCentering_[gridData_.iBx][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iBx][gridData_.idirY]}};
                    case HybridQuantity::Scalar::By:
                        return {{hybridQtyCentering_[gridData_.iBy][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iBy][gridData_.idirY]}};
                    case HybridQuantity::Scalar::Bz:
                        return {{hybridQtyCentering_[gridData_.iBz][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iBz][gridData_.idirY]}};
                    case HybridQuantity::Scalar::Ex:
                        return {{hybridQtyCentering_[gridData_.iEx][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iEx][gridData_.idirY]}};
                    case HybridQuantity::Scalar::Ey:
                        return {{hybridQtyCentering_[gridData_.iEy][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iEy][gridData_.idirY]}};
                    case HybridQuantity::Scalar::Ez:
                        return {{hybridQtyCentering_[gridData_.iEz][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iEz][gridData_.idirY]}};
                    case HybridQuantity::Scalar::Jx:
                        return {{hybridQtyCentering_[gridData_.iJx][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iJx][gridData_.idirY]}};
                    case HybridQuantity::Scalar::Jy:
                        return {{hybridQtyCentering_[gridData_.iJy][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iJy][gridData_.idirY]}};
                    case HybridQuantity::Scalar::Jz:
                        return {{hybridQtyCentering_[gridData_.iJz][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iJz][gridData_.idirY]}};
                    case HybridQuantity::Scalar::rho:
                        return {{hybridQtyCentering_[gridData_.irho][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.irho][gridData_.idirY]}};
                    case HybridQuantity::Scalar::Vx:
                        return {{hybridQtyCentering_[gridData_.iVx][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iVx][gridData_.idirY]}};
                    case HybridQuantity::Scalar::Vy:
                        return {{hybridQtyCentering_[gridData_.iVy][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iVy][gridData_.idirY]}};
                    case HybridQuantity::Scalar::Vz:
                        return {{hybridQtyCentering_[gridData_.iVz][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iVz][gridData_.idirY]}};
                    case HybridQuantity::Scalar::P:
                        return {{hybridQtyCentering_[gridData_.iP][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iP][gridData_.idirY]}};
                    default: throw_runtime_error("Wrong hybridQuantity");
                }
            }

            else if constexpr (dim == 3)
            {
                switch (hybridQuantity)
                {
                    case HybridQuantity::Scalar::Bx:
                        return {{hybridQtyCentering_[gridData_.iBx][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iBx][gridData_.idirY],
                                 hybridQtyCentering_[gridData_.iBx][gridData_.idirZ]}};
                    case HybridQuantity::Scalar::By:
                        return {{hybridQtyCentering_[gridData_.iBy][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iBy][gridData_.idirY],
                                 hybridQtyCentering_[gridData_.iBy][gridData_.idirZ]}};
                    case HybridQuantity::Scalar::Bz:
                        return {{hybridQtyCentering_[gridData_.iBz][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iBz][gridData_.idirY],
                                 hybridQtyCentering_[gridData_.iBz][gridData_.idirZ]}};
                    case HybridQuantity::Scalar::Ex:
                        return {{hybridQtyCentering_[gridData_.iEx][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iEx][gridData_.idirY],
                                 hybridQtyCentering_[gridData_.iEx][gridData_.idirZ]}};
                    case HybridQuantity::Scalar::Ey:
                        return {{hybridQtyCentering_[gridData_.iEy][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iEy][gridData_.idirY],
                                 hybridQtyCentering_[gridData_.iEy][gridData_.idirZ]}};
                    case HybridQuantity::Scalar::Ez:
                        return {{hybridQtyCentering_[gridData_.iEz][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iEz][gridData_.idirY],
                                 hybridQtyCentering_[gridData_.iEz][gridData_.idirZ]}};
                    case HybridQuantity::Scalar::Jx:
                        return {{hybridQtyCentering_[gridData_.iJx][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iJx][gridData_.idirY],
                                 hybridQtyCentering_[gridData_.iJx][gridData_.idirZ]}};
                    case HybridQuantity::Scalar::Jy:
                        return {{hybridQtyCentering_[gridData_.iJy][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iJy][gridData_.idirY],
                                 hybridQtyCentering_[gridData_.iJy][gridData_.idirZ]}};
                    case HybridQuantity::Scalar::Jz:
                        return {{hybridQtyCentering_[gridData_.iJz][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iJz][gridData_.idirY],
                                 hybridQtyCentering_[gridData_.iJz][gridData_.idirZ]}};
                    case HybridQuantity::Scalar::rho:
                        return {{hybridQtyCentering_[gridData_.irho][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.irho][gridData_.idirY],
                                 hybridQtyCentering_[gridData_.irho][gridData_.idirZ]}};
                    case HybridQuantity::Scalar::Vx:
                        return {{hybridQtyCentering_[gridData_.iVx][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iVx][gridData_.idirY],
                                 hybridQtyCentering_[gridData_.iVx][gridData_.idirZ]}};
                    case HybridQuantity::Scalar::Vy:
                        return {{hybridQtyCentering_[gridData_.iVy][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iVy][gridData_.idirY],
                                 hybridQtyCentering_[gridData_.iVy][gridData_.idirZ]}};
                    case HybridQuantity::Scalar::Vz:
                        return {{hybridQtyCentering_[gridData_.iVz][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iVz][gridData_.idirY],
                                 hybridQtyCentering_[gridData_.iVz][gridData_.idirZ]}};
                    case HybridQuantity::Scalar::P:
                        return {{hybridQtyCentering_[gridData_.iP][gridData_.idirX],
                                 hybridQtyCentering_[gridData_.iP][gridData_.idirY],
                                 hybridQtyCentering_[gridData_.iP][gridData_.idirZ]}};
                    default: throw_runtime_error("Wrong hybridQuantity");
                }
            }
        }




        constexpr static std::array<std::array<QtyCentering, dim>, 3>
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


                default: throw_runtime_error("Wrong hybridQuantity");
            }
        }




        auto static constexpr dualToPrimal()
        {
            /*
             * the following is only valid when dual and primal do not have the same number of
            ghosts
             * and that depends on the interp order
             * It is commented out because ghosts are hard coded to 5 for now.
             *
            if constexpr (interp_order == 1 || interp_order == 2 || interp_order == 4)
                return -1;
            else if constexpr (interp_order == 3)
                return 1;
               */
            return -1;
        }




        auto static constexpr primalToDual()
        {
            return 1;
            /*
             * the following is only valid when dual and primal do not have the same number of
            ghosts
             * and that depends on the interp order
             * It is commented out because ghosts are hard coded to 5 for now.
             *
            if constexpr (interp_order == 1 || interp_order == 2 || interp_order == 4)
                return 1;
            else if constexpr (interp_order == 3)
                return -1;
                */
        }




        auto static constexpr momentsToEx()
        {
            // Ex is dual primal primal
            // moments are primal primal primal
            // operation is thus Ppp to Dpp
            // shift only in the X direction
            auto constexpr iShift = primalToDual();

            // P1 is always on top of Ex so no shift

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




        auto static constexpr momentsToEy()
        {
            // Ey is primal dual primal
            // moments are primal primal primal
            // operation is thus pPp to pDp
            // shift only in the Y direction
            [[maybe_unused]] auto constexpr iShift = primalToDual();

            if constexpr (dimension == 1)
            {
                // since the linear combination is in the Y direction
                // in 1D the moment is already on Ey so return 1 point with no shift
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




        auto static constexpr momentsToEz()
        {
            // Ez is primal  primal dual
            // moments are primal primal primal
            // operation is thus ppP to ppD
            // shift only in the Z direction
            [[maybe_unused]] auto constexpr iShift = primalToDual();

            if constexpr (dimension == 1)
            {
                // since the linear combination is in the Z direction
                // in 1D or 2D the moment is already on Ez so return 1 point with no shift
                // with coef 1.
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




        auto static constexpr ExToMoments()
        {
            // Ex is dual primal primal
            // moments are primal primal primal
            // operation is thus Dpp to Ppp
            // shift only in the X direction
            auto constexpr iShift = dualToPrimal();

            if constexpr (dimension == 1)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{iShift}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
            if constexpr (dimension == 2)
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




        auto static constexpr EyToMoments()
        {
            // Ey is       primal dual   primal
            // moments are primal primal primal
            // operation is thus pDp to pPp
            // shift only in the Y direction
            [[maybe_unused]] auto constexpr iShift = dualToPrimal();

            if constexpr (dimension == 1)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 1.0};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            if constexpr (dimension == 2)
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




        auto static constexpr EzToMoments()
        {
            // Ez is       primal primal dual
            // moments are primal primal primal
            // operation is thus ppD to ppP
            // shift only in the Z direction
            [[maybe_unused]] auto constexpr iShift = dualToPrimal();

            if constexpr (dimension == 1)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 1.0};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            if constexpr (dimension == 2)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0}, 1.0};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            else if constexpr (dimension == 3)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0, 0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{0, 0, iShift}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
        }


        auto static constexpr JxToMoments()
        {
            // Jx is dual primal primal
            // moments are primal primal primal
            // operation is thus Dpp to Ppp
            // shift only in the X direction

            auto constexpr iShift = dualToPrimal();

            if constexpr (dimension == 1)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{iShift}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
            if constexpr (dimension == 2)
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


        auto static constexpr JyToMoments()
        {
            // Jy is primal dual primal
            // moments are primal primal primal
            // operation is thus pDp to pPp
            // shift only in the Y direction

            [[maybe_unused]] auto constexpr iShift = dualToPrimal();

            if constexpr (dimension == 1)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 1.0};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            if constexpr (dimension == 2)
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



        auto static constexpr JzToMoments()
        {
            // Jy is primal primal dual
            // moments are primal primal primal
            // operation is thus ppD to ppP
            // shift only in the Z direction

            [[maybe_unused]] auto constexpr iShift = dualToPrimal();

            if constexpr (dimension == 1)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 1.0};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            if constexpr (dimension == 2)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0}, 1.0};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            else if constexpr (dimension == 3)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0, 0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{0, 0, iShift}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
        }



        auto static constexpr ByToEx()
        { // By is dual primal dual
            // Ex is dual primal primal
            // operation is thus dpD to dpP
            // shift only in the Z direction
            [[maybe_unused]] auto constexpr iShift = dualToPrimal();

            if constexpr (dimension == 1)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 1.0};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            if constexpr (dimension == 2)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0}, 1.0};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            else if constexpr (dimension == 3)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0, 0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{0, 0, iShift}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
        }




        auto static constexpr BzToEx()
        {
            // Bz is dual dual primal
            // Ex is dual primal primal
            // operation is thus pDp to pPp
            // shift only in the Y direction
            [[maybe_unused]] auto constexpr iShift = dualToPrimal();

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
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0, 0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{0, iShift, 0}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
        }




        auto static constexpr ByToEz()
        {
            // By is dual primal dual
            // Ez is primal primal dual
            // operation is thus Dpd to Ppd
            // shift only in the X direction
            auto constexpr iShift = dualToPrimal();

            // P1 is always on top of Ez so no shift

            if constexpr (dimension == 1)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{iShift}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
            if constexpr (dimension == 2)
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




        auto static constexpr BxToEz()
        {
            // Bx is primal dual dual
            // Ez is primal primal dual
            // operation is thus pDd to pPd
            // shift only in the Y direction
            [[maybe_unused]] auto constexpr iShift = dualToPrimal();

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
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0, 0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{0, iShift, 0}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
        }




        auto static constexpr BxToEy()
        {
            // Bx is primal dual dual
            // Ey is primal dual primal
            // operation is thus pdD to pdP
            // shift only in the Z direction
            [[maybe_unused]] auto constexpr iShift = dualToPrimal();

            if constexpr (dimension == 1)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 1};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            if constexpr (dimension == 2)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0}, 1.};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            else if constexpr (dimension == 3)
            {
                // in 3D we need two points, the second with a dualToPrimal shift along Z
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0, 0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{0, 0, iShift}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
        }




        auto static constexpr BzToEy()
        {
            // Bz is dual dual primal
            // Ey is primal dual primal
            // operation is thus Ddp to Pdp
            // shift only in the X direction
            auto constexpr iShift = dualToPrimal();

            if constexpr (dimension == 1)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 0.5};
                constexpr WeightPoint<dimension> P2{Point<int, dimension>{iShift}, 0.5};
                return std::array<WeightPoint<dimension>, 2>{P1, P2};
            }
            if constexpr (dimension == 2)
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



        auto static constexpr JxToEx()
        {
            // Jx is dual primal primal
            // Ex is dual primal primal
            // operation is thus dpp to dpp
            // no shift for a yee grid

            if constexpr (dimension == 1)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 1};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            if constexpr (dimension == 2)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0}, 1.0};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            else if constexpr (dimension == 3)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0, 0}, 1.0};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
        }



        auto static constexpr JyToEy()
        {
            // Jy is primal dual primal
            // Ey is primal dual primal
            // operation is thus pdp to pdp
            // no shift for a yee grid

            if constexpr (dimension == 1)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 1};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            if constexpr (dimension == 2)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0}, 1.0};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            else if constexpr (dimension == 3)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0, 0}, 1.0};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
        }



        auto static constexpr JzToEz()
        {
            // Jz is primal primal dual
            // Ez is primal primal dual
            // operation is thus ppd to ppd
            // no shift for a yee grid

            if constexpr (dimension == 1)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0}, 1};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            if constexpr (dimension == 2)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0}, 1.0};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
            else if constexpr (dimension == 3)
            {
                constexpr WeightPoint<dimension> P1{Point<int, dimension>{0, 0, 0}, 1.0};
                return std::array<WeightPoint<dimension>, 1>{P1};
            }
        }
    }; // namespace core


    /*

    template<std::size_t dim, std::size_t interpOrder>
    void constexpr GridLayoutImpl<dim, interpOrder>::initLinearCombinations_()
    {
        // cf https://hephaistos.lpp.polytechnique.fr/redmine/projects/hyb-par/wiki/Ohm
        // for how to calculate coefficients and shift indexes.

        int dualToPrimal = 0;
        int primalTodual = 0;

        if constexpr (interpOrder == 1 || interpOrder == 2 || interpOrder == 4)
        {
            dualToPrimal = -1;
            primalTodual = 1;
        }
        else if constexpr (interpOrder == 3)
        {
            dualToPrimal = 1;
            primalTodual = -1;
        }

        WeightPoint P1;
        WeightPoint P2;

        // moment to Ex is Ppp to Dpp
        // shift only in X
        // the average is done for all simulations
        P1.ix   = 0;
        P1.iy   = 0;
        P1.iz   = 0;
        P1.coef = 0.5;
        P2.ix   = primalTodual;
        P2.iy   = 0;
        P2.iz   = 0;
        P2.coef = 0.5;
        momentsToEx_.push_back(P1);
        momentsToEx_.push_back(P2);


        // moment to Ey is pPp to pDp
        // shift only in Y
        // the average is done only for 2D and 3D simulation
        P1.ix   = 0;
        P1.iy   = 0;
        P1.iz   = 0;
        P1.coef = (dim >= 2) ? 0.5 : 1.;
        momentsToEy_.push_back(P1);

        // in 2 and 3D, add another point and average
        if (dim >= 2)
        {
            P2.ix   = 0;
            P2.iy   = primalTodual;
            P2.iz   = 0;
            P2.coef = 0.5;
            momentsToEy_.push_back(P2);
        }


        // moment to Ez is ppP to ppD
        // shift only in Z
        // the average is done only for 3D simulation
        // hence for 1D and 2D runs coef==1
        P1.ix   = 0;
        P1.iy   = 0;
        P1.iz   = 0;
        P1.coef = (dim == 3) ? 0.5 : 1;
        momentsToEz_.push_back(P1);

        if (dim == 3)
        {
            P2.ix   = 0;
            P2.iy   = 0;
            P2.iz   = primalTodual;
            P2.coef = 0.5;
            momentsToEz_.push_back(P2);
        }




        // Bx to Ey is pdD to pdP
        // shift only in Z
        // the average is done only for 3D simulations
        P1.ix   = 0;
        P1.iy   = 0;
        P1.iz   = 0;
        P1.coef = (dim == 3) ? 0.5 : 1;
        BxToEy_.push_back(P1);

        if (dim == 3)
        {
            P2.ix   = 0;
            P2.iy   = 0;
            P2.iz   = dualToPrimal;
            P2.coef = 0.5;
            BxToEy_.push_back(P2);
        }



        // Bx to Ez is pDd to pPd
        // shift in the Y direction only
        // the average is done for 2D and 3D simulations
        // hence for 1D simulations coef is 1
        P1.ix   = 0;
        P1.iy   = 0;
        P1.iz   = 0;
        P1.coef = (dim >= 2) ? 0.5 : 1;
        BxToEz_.push_back(P1);

        if (dim >= 2)
        {
            P2.ix   = 0;
            P2.iy   = dualToPrimal;
            P2.iz   = 0;
            P2.coef = 0.5;
            BxToEz_.push_back(P2);
        }


        // By to Ex is dpD to dpP
        // shift only in the Z direction
        // averaging is done only for 3D simulations
        P1.ix   = 0;
        P1.iy   = 0;
        P1.iz   = 0;
        P1.coef = (dim == 3) ? 0.5 : 1;
        ByToEx_.push_back(P1);

        if (dim == 3)
        {
            P2.ix   = 0;
            P2.iy   = 0;
            P2.iz   = dualToPrimal;
            P2.coef = 0.5;
            ByToEx_.push_back(P2);
        }

        // By to Ez is Dpd to Ppd
        // shift only in the X direction
        // the averaging is done in all simulations
        P1.ix   = 0;
        P1.iy   = 0;
        P1.iz   = 0;
        P1.coef = 0.5;
        P2.ix   = dualToPrimal;
        P2.iy   = 0;
        P2.iz   = 0;
        P2.coef = 0.5;
        ByToEz_.push_back(P1);
        ByToEz_.push_back(P2);


        // Bz to Ex is dDp to dPp
        // shift only in the Y direction
        // the averaging is done for 2D and 3D simulations
        P1.ix   = 0;
        P1.iy   = 0;
        P1.iz   = 0;
        P1.coef = (dim >= 2) ? 0.5 : 1;
        BzToEx_.push_back(P1);

        if (dim >= 2)
        {
            P2.ix   = 0;
            P2.iy   = dualToPrimal;
            P2.iz   = 0;
            P2.coef = 0.5;
            BzToEx_.push_back(P2);
        }


        // Bz to Ey is Ddp to Pdp
        // shift only in the X direction
        // the averaging is done for all simulations
        P1.ix   = 0;
        P1.iy   = 0;
        P1.iz   = 0;
        P1.coef = 0.5;
        BzToEy_.push_back(P1);
        P2.ix   = dualToPrimal;
        P2.iy   = 0;
        P2.iz   = 0;
        P2.coef = 0.5;
        BzToEy_.push_back(P2);


        // Ex to Moment is Dpp to Ppp
        // shift only in the X direction
        // the averaging is done for all simulations
        P1.ix   = 0;
        P1.iy   = 0;
        P1.iz   = 0;
        P1.coef = 0.5;
        ExToMoment_.push_back(P1);
        P2.ix   = dualToPrimal;
        P2.iy   = 0;
        P2.iz   = 0;
        P2.coef = 0.5;
        ExToMoment_.push_back(P2);


        // Ey to Moment is pDp to PPP
        // shift tis only in the Y direction
        // the averaging is done for 2D and 3D simulations
        P1.ix   = 0;
        P1.iy   = 0;
        P1.iz   = 0;
        P1.coef = (dim >= 2) ? 0.5 : 1;
        EyToMoment_.push_back(P1);

        if (dim >= 2)
        {
            P2.ix   = 0;
            P2.iy   = dualToPrimal;
            P2.iz   = 0;
            P2.coef = 0.5;
            EyToMoment_.push_back(P2);
        }


        // Ez to Moment is ppD on ppP
        // shift only in the Z direction
        // the averaging is only for 3D simulations
        P1.ix   = 0;
        P1.iy   = 0;
        P1.iz   = 0;
        P1.coef = (dim == 3) ? 0.5 : 1;
        EzToMoment_.push_back(P1);

        if (dim == 3)
        {
            P2.ix   = 0;
            P2.iy   = 0;
            P2.iz   = dualToPrimal;
            P2.coef = 0.5;
            EzToMoment_.push_back(P2);
        }
    }




    template<std::size_t dim>
    LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::momentsToEx() const
    {
        return this->momentsToEx_;
    }


    template<std::size_t dim>
    LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::momentsToEy() const
    {
        return this->momentsToEy_;
    }


    template<std::size_t dim>
    LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::momentsToEz() const
    {
        return this->momentsToEz_;
    }

    template<std::size_t dim>
    LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::ByToEx() const
    {
        return this->ByToEx_;
    }

    template<std::size_t dim>
    LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::ByToEz() const
    {
        return this->ByToEz_;
    }

    template<std::size_t dim>
    LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::BxToEy() const
    {
        return this->BxToEy_;
    }

    template<std::size_t dim>
    LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::BxToEz() const
    {
        return this->BxToEz_;
    }

    template<std::size_t dim>
    LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::BzToEy() const
    {
        return this->BzToEy_;
    }

    template<std::size_t dim>
    LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::BzToEx() const
    {
        return this->BzToEx_;
    }



    template<std::size_t dim>
    LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::ExToMoment() const
    {
        return this->ExToMoment_;
    }


    template<std::size_t dim>
    LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::EyToMoment() const
    {
        return this->EyToMoment_;
    }


    template<std::size_t dim>
    LinearCombination const& GridLayoutImpl<Layout::Yee, dim>::EzToMoment() const
    {
        return this->EzToMoment_;
    }

    */

} // namespace core
} // namespace PHARE


#endif // PHARE_CORE_GRID_GRIDLAYOUTYEE_H
