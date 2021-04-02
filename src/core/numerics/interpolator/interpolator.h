#ifndef PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATOR_H
#define PHARE_CORE_NUMERICS_INTERPOLATOR_INTERPOLATOR_H



#include <array>
#include <cstddef>

#include "core/data/grid/gridlayout.h"
#include "core/data/vecfield/vecfield_component.h"
#include "core/utilities/point/point.h"

#include "core/logger.h"

namespace PHARE
{
namespace core
{
    //! return the size of the index and weights arrays for
    //! interpolation at a given order
    // Number of points where the interpolator of order interpOrder
    // will deposit mass/momentum. This is hence the size of the
    // index and weight arrays for interpolation at a given order
    constexpr int nbrPointsSupport(int interpOrder) { return interpOrder + 1; }


    //! return the number of points on which to interpolate at order InterpOrder
    template<std::size_t interpOrder>
    int computeStartIndex(double normalizedPos)
    {
        return static_cast<int>(normalizedPos - (static_cast<double>(interpOrder) - 1.) / 2.);
    }




    /** \brief the class Weight aims at computing the weight coefficient for
     *  interpolation at a specific order
     *
     *  This class assumes the interpolation order is known at compile-time
     *  thus there are three specialization for orders 1, 2 and 3.
     *
     *  the class has only one method called computeWeight that takes three arguments:
     *
     *  \param[in] normalized position is the particle position in the cell normalized by grid
     * spacing \param[in] startIndex first grid index where to interpolate the field \param[out]
     * weights contains the nbrPointsSupport weights calculated
     *
     *  these three parameters are given for a specific direction (x, y or z)
     */
    template<std::size_t interpOrder>
    class Weighter
    {
    };



    /** \brief Specialization of Weight for first order interpolation
     */
    template<>
    class Weighter<1>
    {
    public:
        inline void computeWeight(double normalizedPos, int startIndex,
                                  std::array<double, nbrPointsSupport(1)>& weights)
        {
            weights[1] = normalizedPos - static_cast<double>(startIndex);
            weights[0] = 1. - weights[1];
        }

        static const int interp_order = 1;
    };


    /** \brief specialization of Weighter for second order interpolation
     */
    template<>
    class Weighter<2>
    {
    public:
        inline void computeWeight(double normalizedPos, int startIndex,
                                  std::array<double, nbrPointsSupport(2)>& weights)
        {
            auto index = startIndex + 1;
            auto delta = static_cast<double>(index) - normalizedPos;
            double coef1, coef2, coef3;
            coef1 = 0.5 + delta;
            coef2 = delta;
            coef3 = 0.5 - delta;

            weights[0] = 0.5 * coef1 * coef1;
            weights[1] = 0.75 - coef2 * coef2;
            weights[2] = 0.5 * coef3 * coef3;
        }

        static const int interp_order = 2;
    };



    /** \brief specialization of Weighter for third order interpolation
     */
    template<>
    class Weighter<3>
    {
    public:
        inline void computeWeight(double normalizedPos, int startIndex,
                                  std::array<double, nbrPointsSupport(3)>& weights)
        {
            double coef1, coef2, coef3, coef4;
            auto index = static_cast<double>(startIndex) - normalizedPos;
            coef1      = 1. + 0.5 * index;
            coef2      = index + 1;
            coef3      = index + 2;
            coef4      = 1. - 0.5 * (index + 3);

            double coef2_sq  = coef2 * coef2;
            double coef2_cub = coef2_sq * coef2;
            double coef3_sq  = coef3 * coef3;
            double coef3_cub = coef3_sq * coef3;

            weights[0] = (4. / 3.) * coef1 * coef1 * coef1;
            weights[1] = 2. / 3. - coef2_sq - 0.5 * coef2_cub;
            weights[2] = 2. / 3. - coef3_sq + 0.5 * coef3_cub;
            weights[3] = (4. / 3.) * coef4 * coef4 * coef4;
        }

        static const int interp_order = 3;
    };




    //! Interpol performs the interpolation of a field using precomputed weights at
    //! indices starting at startIndex. The class is templated by the Dimensionality
    template<std::size_t dim>
    class MeshToParticle
    {
    };



    /** \brief specialization of Interpol for 1D interpolation
     */
    template<>
    class MeshToParticle<1>
    {
    public:
        /** Performs the 1D interpolation
         * \param[in] field is the field from which values are interpolated
         * \param[in] fieldCentering is the centering (dual or primal) of the field
         * \param[in] startIndex is the first of the nbrPointsSupport indices where to interpolate
         * the field \param[in] weights are the nbrPointsSupport weights used for the interpolation
         */
        template<typename Field, typename Array1, typename Array2>
        inline double operator()(Field const& field,
                                 std::array<QtyCentering, 1> const& fieldCentering,
                                 Array1 const& startIndex, Array2 const& weights)
        {
            auto fieldAtParticle    = 0.;
            auto const& xStartIndex = startIndex[static_cast<int>(fieldCentering[0])][0];
            auto const& xWeights    = weights[static_cast<int>(fieldCentering[0])][0];
            auto order_size         = xWeights.size();

            for (auto ik = 0u; ik < order_size; ++ik)
            {
                fieldAtParticle += field(xStartIndex + ik) * xWeights[ik];
            }
            return fieldAtParticle;
        }
    };


    /**\brief Specialization of Interpol for 2D interpolation
     */
    template<>
    class MeshToParticle<2>
    {
    public:
        /** Performs the 2D interpolation
         * \param[in] field is the field from which values are interpolated
         * \param[in] fieldCentering is the centering (dual or primal) of the field in each
         * direction \param[in] startIndex is the first of the nbrPointsSupport indices where to
         * interpolate the field in both directions \param[in] weights are the nbrPointsSupport
         * weights used for the interpolation in both directions
         */

        template<typename Field, typename Array1, typename Array2>
        inline double operator()(Field const& field,
                                 std::array<QtyCentering, 2> const& fieldCentering,
                                 Array1 const& startIndex, Array2 const& weights)
        {
            auto const& xStartIndex = startIndex[static_cast<int>(fieldCentering[0])][0];
            auto const& yStartIndex = startIndex[static_cast<int>(fieldCentering[1])][1];
            auto const& xWeights    = weights[static_cast<int>(fieldCentering[0])][0];
            auto const& yWeights    = weights[static_cast<int>(fieldCentering[1])][1];

            auto order_size        = xWeights.size();
            double fieldAtParticle = 0.;
            for (auto ix = 0u; ix < order_size; ++ix)
            {
                double Yinterp = 0.;
                for (auto iy = 0u; iy < order_size; ++iy)
                {
                    Yinterp += field(xStartIndex + ix, yStartIndex + iy) * yWeights[iy];
                }
                fieldAtParticle += Yinterp * xWeights[ix];
            }

            return fieldAtParticle;
        }
    };



    /** \brief Specialization of Interpol for 3D interpolation
     */
    template<>
    class MeshToParticle<3>
    {
    public:
        /** Performs the 3D interpolation
         * \param[in] field is the field from which values are interpolated
         * \param[in] fieldCentering is the centering (dual or primal) of the field in each
         * direction \param[in] startIndex is the first of the nbrPointsSupport indices where to
         * interpolate the field in the 3 directions \param[in] weights are the nbrPointsSupport
         * weights used for the interpolation in the 3 directions
         */
        template<typename Field, typename Array1, typename Array2>
        inline double operator()(Field const& field,
                                 std::array<QtyCentering, 3> const& fieldCentering,
                                 Array1 const& startIndex, Array2 const& weights)
        {
            auto const& xStartIndex = startIndex[static_cast<std::size_t>(fieldCentering[0])][0];
            auto const& yStartIndex = startIndex[static_cast<std::size_t>(fieldCentering[1])][1];
            auto const& zStartIndex = startIndex[static_cast<std::size_t>(fieldCentering[2])][2];
            auto const& xWeights    = weights[static_cast<std::size_t>(fieldCentering[0])][0];
            auto const& yWeights    = weights[static_cast<std::size_t>(fieldCentering[1])][1];
            auto const& zWeights    = weights[static_cast<std::size_t>(fieldCentering[2])][2];

            auto order_size        = xWeights.size();
            double fieldAtParticle = 0.;
            for (auto ix = 0u; ix < order_size; ++ix)
            {
                double Yinterp = 0.;
                for (auto iy = 0u; iy < order_size; ++iy)
                {
                    double Zinterp = 0.;
                    for (auto iz = 0u; iz < order_size; ++iz)
                    {
                        Zinterp += field(xStartIndex + ix, yStartIndex + iy, zStartIndex + iz)
                                   * zWeights[iz];
                    }
                    Yinterp += Zinterp * yWeights[iy];
                }
                fieldAtParticle += Yinterp * xWeights[ix];
            }
            return fieldAtParticle;
        }
    };




    //! ParticleToMesh projects a particle density and flux to given grids
    template<std::size_t dim>
    class ParticleToMesh
    {
    };



    /** \brief specialization of ParticleToMesh for 1D interpolation
     */
    template<>
    class ParticleToMesh<1>
    {
    public: /** Performs the 1D interpolation
             * \param[in] density is the field that will be interpolated from the particle Particle
             * \param[in] xFlux is the field that will be interpolated from the particle Particle
             * \param[in] yFlux is the field that will be interpolated from the particle Particle
             * \param[in] zFlux is the field that will be interpolated from the particle Particle
             * \param[in] fieldCentering is the centering (dual or primal) of the field in each
             * direction \param[in] particle is the single particle used for the interpolation of
             * density and flux \param[in] startIndex is the first index for which a particle will
             * contribute \param[in] weights is the arrays of weights for the associated index
             */
        template<typename Field, typename Array1, typename Array2, typename Particle,
                 typename VectorCenteringArray>
        inline void operator()(Field& density, Field& xFlux, Field& yFlux, Field& zFlux,
                               std::array<QtyCentering, 1> const& densityCentering,
                               VectorCenteringArray const& fluxCentering, Particle const& particle,
                               Array1 const& startIndex, Array2 const& weights, double coef = 1.)
        {
            auto const& xDenStartIndex = startIndex[static_cast<int>(densityCentering[0])][0];
            auto const& xDenWeights    = weights[static_cast<int>(densityCentering[0])][0];

            auto const& xXFluxStartIndex = startIndex[static_cast<int>(fluxCentering[0][0])][0];
            auto const& xXFluxWeights    = weights[static_cast<int>(fluxCentering[0][0])][0];

            auto const& xYFluxStartIndex = startIndex[static_cast<int>(fluxCentering[1][0])][0];
            auto const& xYFluxWeights    = weights[static_cast<int>(fluxCentering[1][0])][0];

            auto const& xZFluxStartIndex = startIndex[static_cast<int>(fluxCentering[2][0])][0];
            auto const& xZFluxWeights    = weights[static_cast<int>(fluxCentering[2][0])][0];

            auto order_size = xDenWeights.size();

            auto const partRho   = particle.weight;
            auto const xPartFlux = particle.v[0] * particle.weight;
            auto const yPartFlux = particle.v[1] * particle.weight;
            auto const zPartFlux = particle.v[2] * particle.weight;

            for (auto ik = 0u; ik < order_size; ++ik)
            {
                density(xDenStartIndex + ik) += partRho * xDenWeights[ik] * coef;

                xFlux(xXFluxStartIndex + ik) += xPartFlux * xXFluxWeights[ik] * coef;
                yFlux(xYFluxStartIndex + ik) += yPartFlux * xYFluxWeights[ik] * coef;
                zFlux(xZFluxStartIndex + ik) += zPartFlux * xZFluxWeights[ik] * coef;
            }
        }
    };




    /** \brief specialization of ParticleToMesh for 2D interpolation
     */
    template<>
    class ParticleToMesh<2>
    {
    public: /** Performs the 2D interpolation
             * \param[in] density is the field that will be interpolated from the particle Particle
             * \param[in] xFlux is the field that will be interpolated from the particle Particle
             * \param[in] yFlux is the field that will be interpolated from the particle Particle
             * \param[in] zFlux is the field that will be interpolated from the particle Particle
             * \param[in] fieldCentering is the centering (dual or primal) of the field in each
             * direction \param[in] particle is the single particle used for the interpolation of
             * density and flux \param[in] startIndex is the first index for which a particle will
             * contribute \param[in] weights is the arrays of weights for the associated index
             */
        template<typename Field, typename Array1, typename Array2, typename Particle,
                 typename VectorCenteringArray>
        inline void operator()(Field& density, Field& xFlux, Field& yFlux, Field& zFlux,
                               std::array<QtyCentering, 2> const& densityCentering,
                               VectorCenteringArray const& fluxCentering, Particle const& particle,
                               Array1 const& startIndex, Array2 const& weights, double coef = 1.)
        {
            auto const& xDenStartIndex = startIndex[static_cast<int>(densityCentering[0])][0];
            auto const& xDenWeights    = weights[static_cast<int>(densityCentering[0])][0];

            auto const& xXFluxStartIndex = startIndex[static_cast<int>(fluxCentering[0][0])][0];
            auto const& xXFluxWeights    = weights[static_cast<int>(fluxCentering[0][0])][0];

            auto const& xYFluxStartIndex = startIndex[static_cast<int>(fluxCentering[1][0])][0];
            auto const& xYFluxWeights    = weights[static_cast<int>(fluxCentering[1][0])][0];

            auto const& xZFluxStartIndex = startIndex[static_cast<int>(fluxCentering[2][0])][0];
            auto const& xZFluxWeights    = weights[static_cast<int>(fluxCentering[2][0])][0];




            auto const& yDenStartIndex = startIndex[static_cast<int>(densityCentering[0])][1];
            auto const& yDenWeights    = weights[static_cast<int>(densityCentering[0])][1];

            auto const& yXFluxStartIndex = startIndex[static_cast<int>(fluxCentering[0][1])][1];
            auto const& yXFluxWeights    = weights[static_cast<int>(fluxCentering[0][1])][1];

            auto const& yYFluxStartIndex = startIndex[static_cast<int>(fluxCentering[1][1])][1];
            auto const& yYFluxWeights    = weights[static_cast<int>(fluxCentering[1][1])][1];

            auto const& yZFluxStartIndex = startIndex[static_cast<int>(fluxCentering[2][1])][1];
            auto const& yZFluxWeights    = weights[static_cast<int>(fluxCentering[2][1])][1];



            auto const partRho   = particle.weight * coef;
            auto const xPartFlux = particle.v[0] * particle.weight * coef;
            auto const yPartFlux = particle.v[1] * particle.weight * coef;
            auto const zPartFlux = particle.v[2] * particle.weight * coef;

            auto order_size = xDenWeights.size();
            for (auto ix = 0u; ix < order_size; ++ix)
            {
                for (auto iy = 0u; iy < order_size; ++iy)
                {
                    density(xDenStartIndex + ix, yDenStartIndex + iy)
                        += partRho * xDenWeights[ix] * yDenWeights[iy];

                    xFlux(xXFluxStartIndex + ix, yXFluxStartIndex + iy)
                        += xPartFlux * xXFluxWeights[ix] * yXFluxWeights[iy];

                    yFlux(xYFluxStartIndex + ix, yYFluxStartIndex + iy)
                        += yPartFlux * xYFluxWeights[ix] * yYFluxWeights[iy];

                    zFlux(xZFluxStartIndex + ix, yZFluxStartIndex + iy)
                        += zPartFlux * xZFluxWeights[ix] * yZFluxWeights[iy];
                }
            }
        }
    };




    /** \brief specialization of ParticleToMesh for 3D interpolation
     */
    template<>
    class ParticleToMesh<3>
    {
    public: /** Performs the 3D interpolation
             * \param[in] density is the field that will be interpolated from the particle Particle
             * \param[in] xFlux is the field that will be interpolated from the particle Particle
             * \param[in] yFlux is the field that will be interpolated from the particle Particle
             * \param[in] zFlux is the field that will be interpolated from the particle Particle
             * \param[in] fieldCentering is the centering (dual or primal) of the field in each
             * direction \param[in] particle is the single particle used for the interpolation of
             * density and flux \param[in] startIndex is the first index for which a particle will
             * contribute \param[in] weights is the arrays of weights for the associated index
             */
        template<typename Field, typename Array1, typename Array2, typename Particle,
                 typename VectorCenteringArray>
        inline void operator()(Field& density, Field& xFlux, Field& yFlux, Field& zFlux,
                               std::array<QtyCentering, 3> const& densityCentering,
                               VectorCenteringArray const& fluxCentering, Particle const& particle,
                               Array1 const& startIndex, Array2 const& weights, double coef = 1.)
        {
            auto const& xDenStartIndex = startIndex[static_cast<int>(densityCentering[0])][0];
            auto const& xDenWeights    = weights[static_cast<int>(densityCentering[0])][0];

            auto const& xXFluxStartIndex = startIndex[static_cast<int>(fluxCentering[0][0])][0];
            auto const& xXFluxWeights    = weights[static_cast<int>(fluxCentering[0][0])][0];

            auto const& xYFluxStartIndex = startIndex[static_cast<int>(fluxCentering[1][0])][0];
            auto const& xYFluxWeights    = weights[static_cast<int>(fluxCentering[1][0])][0];

            auto const& xZFluxStartIndex = startIndex[static_cast<int>(fluxCentering[2][0])][0];
            auto const& xZFluxWeights    = weights[static_cast<int>(fluxCentering[2][0])][0];




            auto const& yDenStartIndex = startIndex[static_cast<int>(densityCentering[0])][1];
            auto const& yDenWeights    = weights[static_cast<int>(densityCentering[0])][1];

            auto const& yXFluxStartIndex = startIndex[static_cast<int>(fluxCentering[0][1])][1];
            auto const& yXFluxWeights    = weights[static_cast<int>(fluxCentering[0][1])][1];

            auto const& yYFluxStartIndex = startIndex[static_cast<int>(fluxCentering[1][1])][1];
            auto const& yYFluxWeights    = weights[static_cast<int>(fluxCentering[1][1])][1];

            auto const& yZFluxStartIndex = startIndex[static_cast<int>(fluxCentering[2][1])][1];
            auto const& yZFluxWeights    = weights[static_cast<int>(fluxCentering[2][1])][1];




            auto const& zDenStartIndex = startIndex[static_cast<int>(densityCentering[0])][2];
            auto const& zDenWeights    = weights[static_cast<int>(densityCentering[0])][2];

            auto const& zXFluxStartIndex = startIndex[static_cast<int>(fluxCentering[0][2])][2];
            auto const& zXFluxWeights    = weights[static_cast<int>(fluxCentering[0][2])][2];

            auto const& zYFluxStartIndex = startIndex[static_cast<int>(fluxCentering[1][2])][2];
            auto const& zYFluxWeights    = weights[static_cast<int>(fluxCentering[1][2])][2];

            auto const& zZFluxStartIndex = startIndex[static_cast<int>(fluxCentering[2][2])][2];
            auto const& zZFluxWeights    = weights[static_cast<int>(fluxCentering[2][2])][2];

            auto const partRho   = particle.weight * coef;
            auto const xPartFlux = particle.v[0] * particle.weight * coef;
            auto const yPartFlux = particle.v[1] * particle.weight * coef;
            auto const zPartFlux = particle.v[2] * particle.weight * coef;

            auto order_size = xDenWeights.size();
            for (auto ix = 0u; ix < order_size; ++ix)
            {
                for (auto iy = 0u; iy < order_size; ++iy)
                {
                    for (auto iz = 0u; iz < order_size; ++iz)
                    {
                        density(xDenStartIndex + ix, yDenStartIndex + iy, zDenStartIndex + iz)
                            += partRho * xDenWeights[ix] * yDenWeights[iy] * zDenWeights[iz];

                        xFlux(xXFluxStartIndex + ix, yXFluxStartIndex + iy, zXFluxStartIndex + iz)
                            += xPartFlux * xXFluxWeights[ix] * yXFluxWeights[iy]
                               * zXFluxWeights[iz];

                        yFlux(xYFluxStartIndex + ix, yYFluxStartIndex + iy, zYFluxStartIndex + iz)
                            += yPartFlux * xYFluxWeights[ix] * yYFluxWeights[iy]
                               * zYFluxWeights[iz];

                        zFlux(xZFluxStartIndex + ix, yZFluxStartIndex + iy, zZFluxStartIndex + iz)
                            += zPartFlux * xZFluxWeights[ix] * yZFluxWeights[iy]
                               * zZFluxWeights[iz];
                    }
                }
            }
        }
    };




    /** \brief Interpolator is used to perform particle-mesh interpolations using
     * 1st, 2nd or 3rd order interpolation in 1D, 2D or 3D, on a given layout.
     */
    template<std::size_t dim, std::size_t interpOrder>
    class Interpolator : private Weighter<interpOrder>
    {
    public:
        auto static constexpr interp_order = interpOrder;
        auto static constexpr dimension    = dim;
        /**\brief interpolate electromagnetic fields on all particles in the range
         *
         * For each particle :
         *  - The function first calculates the startIndex and weights for interpolation at
         * order InterpOrder and in dimension dim for dual and primal nodes
         *  - then it uses Interpol<> to calculate the interpolation of E and B components
         * onto the particle.
         */
        template<typename PartIterator, typename Electromag, typename GridLayout>
        inline void operator()(PartIterator begin, PartIterator end, Electromag const& Em,
                               GridLayout const& layout)
        {
            PHARE_LOG_SCOPE("Interpolator::operator()");

            // this lambda calculates the startIndex and the nbrPointsSupport() weights for
            // dual field interpolation and puts this at the corresponding location
            // in 'startIndex' and 'weights'. For dual fields, the normalizedPosition
            // is offseted compared to primal ones.
            auto indexAndWeightDual = [this, &layout](auto const& part) {
                for (auto iDim = 0u; iDim < dimension; ++iDim)
                {
                    auto iCell           = layout.AMRToLocal(Point{part.iCell});
                    double normalizedPos = iCell[iDim] + part.delta[iDim] + dualOffset(interpOrder);

                    startIndex_[centering2int(QtyCentering::dual)][iDim]
                        = computeStartIndex<interpOrder>(normalizedPos);

                    weightComputer_.computeWeight(
                        normalizedPos, startIndex_[centering2int(QtyCentering::dual)][iDim],
                        weights_[centering2int(QtyCentering::dual)][iDim]);
                }
            };

            // does the same as above but for a primal field
            auto indexAndWeightPrimal = [this, &layout](auto const& part) {
                for (auto iDim = 0u; iDim < dimension; ++iDim)
                {
                    auto iCell           = layout.AMRToLocal(Point{part.iCell});
                    double normalizedPos = iCell[iDim] + part.delta[iDim];

                    startIndex_[centering2int(QtyCentering::primal)][iDim]
                        = computeStartIndex<interpOrder>(normalizedPos);

                    weightComputer_.computeWeight(
                        normalizedPos, startIndex_[centering2int(QtyCentering::primal)][iDim],
                        weights_[centering2int(QtyCentering::primal)][iDim]);
                }
            };

            auto const& Ex = Em.E.getComponent(Component::X);
            auto const& Ey = Em.E.getComponent(Component::Y);
            auto const& Ez = Em.E.getComponent(Component::Z);
            auto const& Bx = Em.B.getComponent(Component::X);
            auto const& By = Em.B.getComponent(Component::Y);
            auto const& Bz = Em.B.getComponent(Component::Z);

            auto const ExCentering = GridLayout::centering(HybridQuantity::Scalar::Ex);
            auto const EyCentering = GridLayout::centering(HybridQuantity::Scalar::Ey);
            auto const EzCentering = GridLayout::centering(HybridQuantity::Scalar::Ez);
            auto const BxCentering = GridLayout::centering(HybridQuantity::Scalar::Bx);
            auto const ByCentering = GridLayout::centering(HybridQuantity::Scalar::By);
            auto const BzCentering = GridLayout::centering(HybridQuantity::Scalar::Bz);


            // for each particle, first calculate the startIndex and weights
            // for dual and primal quantities.
            // then, knowing the centering (primal or dual) of each electromagnetic
            // component, we use Interpol to actually perform the interpolation.
            // the trick here is that the StartIndex and weights have only been calculated
            // twice, and not for each E,B component.

            PHARE_LOG_START("MeshToParticle::operator()");
            for (auto currPart = begin; currPart != end; ++currPart)
            {
                indexAndWeightPrimal(*currPart);
                indexAndWeightDual(*currPart);

                currPart->Ex = meshToParticle_(Ex, ExCentering, startIndex_, weights_);
                currPart->Ey = meshToParticle_(Ey, EyCentering, startIndex_, weights_);
                currPart->Ez = meshToParticle_(Ez, EzCentering, startIndex_, weights_);
                currPart->Bx = meshToParticle_(Bx, BxCentering, startIndex_, weights_);
                currPart->By = meshToParticle_(By, ByCentering, startIndex_, weights_);
                currPart->Bz = meshToParticle_(Bz, BzCentering, startIndex_, weights_);
            }
            PHARE_LOG_STOP("MeshToParticle::operator()");
        }




        /**\brief interpolate electromagnetic fields on all particles in the range
         *
         * For each particle :
         *  - The function first calculates the startIndex and weights for interpolation at
         * order InterpOrder and in dimension dim for dual and primal nodes
         *  - then it uses Interpol<> to calculate the interpolation of E and B components
         * onto the particle.
         */
        template<typename PartIterator, typename VecField, typename GridLayout,
                 typename Field = typename VecField::field_type>
        inline void operator()(PartIterator begin, PartIterator end, Field& density, VecField& flux,
                               GridLayout const& layout, double coef = 1.)
        {
            // this lambda calculates the startIndex and the order+1 weights for
            // dual field interpolation and puts this at the corresponding location
            // in 'startIndex' and 'weights'. For dual fields, the normalizedPosition
            // is offseted compared to primal ones.
            auto indexAndWeightDual = [this, &layout](auto const& part) {
                for (auto iDim = 0u; iDim < dimension; ++iDim)
                {
                    auto iCell           = layout.AMRToLocal(Point{part.iCell});
                    double normalizedPos = iCell[iDim] + part.delta[iDim] + dualOffset(interpOrder);

                    startIndex_[centering2int(QtyCentering::dual)][iDim]
                        = computeStartIndex<interpOrder>(normalizedPos);

                    weightComputer_.computeWeight(
                        normalizedPos, startIndex_[centering2int(QtyCentering::dual)][iDim],
                        weights_[centering2int(QtyCentering::dual)][iDim]);
                }
            };

            // does the same as above but for a primal field
            auto indexAndWeightPrimal = [this, &layout](auto const& part) {
                for (auto iDim = 0u; iDim < dimension; ++iDim)
                {
                    auto iCell           = layout.AMRToLocal(Point{part.iCell});
                    double normalizedPos = iCell[iDim] + part.delta[iDim];

                    startIndex_[centering2int(QtyCentering::primal)][iDim]
                        = computeStartIndex<interpOrder>(normalizedPos);

                    weightComputer_.computeWeight(
                        normalizedPos, startIndex_[centering2int(QtyCentering::primal)][iDim],
                        weights_[centering2int(QtyCentering::primal)][iDim]);
                }
            };

            auto& xFlux = flux.getComponent(Component::X);
            auto& yFlux = flux.getComponent(Component::Y);
            auto& zFlux = flux.getComponent(Component::Z);

            auto constexpr densityCentering = GridLayout::centering(HybridQuantity::Scalar::rho);
            auto constexpr fluxCentering    = GridLayout::centering(HybridQuantity::Vector::V);


            // for each particle, first calculate the startIndex and weights
            // for dual and primal quantities.
            // then, knowing the centering (primal or dual) of each electromagnetic
            // component, we use Interpol to actually perform the interpolation.
            // the trick here is that the StartIndex and weights have only been calculated
            // twice, and not for each E,B component.

            PHARE_LOG_START("ParticleToMesh::operator()");
            for (auto currPart = begin; currPart != end; ++currPart)
            {
                // TODO #3375
                indexAndWeightPrimal(*currPart);
                indexAndWeightDual(*currPart);

                particleToMesh_(density, xFlux, yFlux, zFlux, densityCentering, fluxCentering,
                                *currPart, startIndex_, weights_, coef);
            }
            PHARE_LOG_STOP("ParticleToMesh::operator()");
        }




    private:
        static_assert(dimension <= 3 && dimension > 0 && interpOrder >= 1 && interpOrder <= 3,
                      "error");

        Weighter<interpOrder> weightComputer_;
        MeshToParticle<dimension> meshToParticle_;
        ParticleToMesh<dimension> particleToMesh_;

        // array[dual/primal][dim]
        std::array<std::array<int, dimension>, 2> startIndex_;
        std::array<std::array<std::array<double, nbrPointsSupport(interpOrder)>, dimension>, 2>
            weights_;

        /**
         * @brief dualOffset returns the offset by which changing the
         * startIndex for dual node interpolation. This offset depends on
         * interpolation order. */
        inline constexpr double dualOffset(int order)
        {
            std::array<double, 3> offsets = {{-0.5, -0.5, 0.5}};
            return offsets[static_cast<std::array<double, 3>::size_type>(order - 1)];
        }
    };


} // namespace core

} // namespace PHARE

#endif
