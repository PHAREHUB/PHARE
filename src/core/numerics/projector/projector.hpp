#ifndef PHARE_CORE_NUMERICS_PROJECTOR_HPP
#define PHARE_CORE_NUMERICS_PROJECTOR_HPP

#include <cmath>
#include <iostream>
#include <cstddef>

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayout_utils.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/data/particles/particle.hpp"

#include "core/numerics/interpolator/interpolator.hpp"


namespace PHARE::core
{

double one_third = 1./3.;

template<typename GridLayout>
class Projector : public LayoutHolder<GridLayout>
{
    constexpr static auto dimension = GridLayout::dimension;
    using LayoutHolder<GridLayout>::layout_;

public:
    // Return J associated with a single particle. 
    // Will need to be called in the solver as it depends on its old and new positions.
    template<typename VecField, typename GridLayout, typename ParticleRange>
    inline void operator()(VecField& J, ParticleRange& rangeOut, ParticleRange& rangeIn, double dt)
    {
        for (auto inIdx = rangeIn.ibegin(), outIdx = rangeOut.ibegin(); inIdx < rangeIn.iend();
                    ++inIdx, ++outIdx)
                {
                    test<GridLayout::dimension>(J, outParticles[outIdx], inParticles[inIdx], dt);
                }   
    }

} // END Projector



// Calculate weights for a single particle, for use in calculating Esirkepov coefficients
void E_weights_(Particle const& part)
{   
    auto const& interpOrder = PHARE::core::Weighter::interp_order;
    auto const& supportpts  = PHARE::core::nbrPointsSupport(interpOrder);
    
    std::vector<std::vector<double>> weights(dimension, std::vector<double>(supportpts));
    
    for (auto i = 0u; i < dimension; ++i)
    {
        auto position = part.iCell[i] + part.delta[i] - 0.5;
        auto startIndex = part.iCell[i] - PHARE::core::computeStartLeftShift<QtyCentering, QtyCentering::dual>(Part.delta[i]);
        PHARE::core::computeWeight(position, startIndex, weights[i]);

        // Add 0.0 to the beginning and end of the weights for calculation purposes
        weights[i].insert(weights[i].begin(), 0.0);
        weights[i].insert(weights[i].end(), 0.0);

    }

    return weights;

} // END E_weights_



// call using test(J, oldCentralNode, E_weights(partOut, partIn))
template<std::size_t dim>
class test
{
};


/** \brief specialization of ParticleToMesh for 1D interpolation of J */
template<>
class test<1>
{
public:
    template<typename VecField, typename Particle>
    inline void operator()(VecField& J, Particle const& partIn,
                           Particle const& partOut, double dt)
    {
        auto& Jx = J(Component::X);
        auto& Jy = J(Component::Y);
        auto& Jz = J(Component::Z);

        auto const& xStartIndex = oldPart.iCell ; // central dual node
        auto const& S0          = E_weights_(partOut);
        auto const& S1          = E_weights_(partIn);
        auto const& order_size  = S0.size();

        // requisite for appropriate centering (greatest weight at node considered)
        int iCorr = order_size/2;
        if Sx0[1] > Sx0[Sx0.size()-2] { 
            iCorr -= 1
        }

        double inv_cell_volume = 1.0 / layout_->cellVolume(); // TODO: see if this is necessary, ie. if the weight is relative to initial cell volume or not
        double charge_weight = inv_cell_volume * partIn.charge * partIn.weight;
            
        double cr_p = charge_weight/dt;   // current density in the evaluated dimension, assuming dx=1
        double cry_p_1D = charge_weight*partIn.v[1];   // current density in the y-direction in 1D
        double crz_p_1D2D = charge_weight*partIn.v[2]; // current density in the z-direction in 1D or 2D

        std::vector<double> Jx_p(order_size, 0.);


        for (auto i = 0u; i < order_size; ++i)
        {
            x = xStartIndex + i - iCorr; // eg, i from -2 to 2 for 3rd order B-splines.

            Wl[i] = S0[i] - S1[i];
            Wt[i] = 0.5 * ( S0[i] + S1[i] );

            Jx_p[i] = Jx_p[i-1] + cr_p * Wl[i-1];
            Jx(x) += Jx_p[i];
            Jy(x)  += cry_p_1D * Wt[i];
            Jz(x)  += crz_p_1D2D * Wt[i];
        }
    }
}; // END test<1> specialization


/** \brief specialization of ParticleToMesh for 2D interpolation of J */
template<>
class test<2>
{
public:
    template<typename VecField, typename Particle>
    inline void operator()(VecField& J, Particle const& partIn,
                           Particle const& partOut, double dt)
    {
        auto& Jx = J(Component::X);
        auto& Jy = J(Component::Y);
        auto& Jz = J(Component::Z);

        auto const& xStartIndex = oldPart.iCell[0];
        auto const& yStartIndex = oldPart.iCell[1];
        auto const& oldWeights    = E_weights_(partOut);
        auto const& newWeights    = E_weights_(partIn);

        auto const& [Sx0, Sy0]      = oldWeights;
        auto const& [Sx1, Sy1]      = newWeights;
        auto const& order_size      = Sx0.size();

        // requisite for appropriate centering
        int iCorr = order_size/2;
        if Sx0[1] > Sx0[Sx0.size()-2] { 
            iCorr -= 1
        }
        int jCorr = order_size/2;
        if Sy0[1] > Sy0[Sy0.size()-2] {
            jCorr -= 1
        }

        double inv_cell_volume = 1.0 / layout_->cellVolume();
        double charge_weight = inv_cell_volume * partIn.charge * partIn.weight;
            
        double cr_p = charge_weight/dt;  // current density in the evaluated dimension assuming dx=1
        double crz_p_1D2D = charge_weight*partIn.v[2]; // current density in the z-direction in 1D or 2D

        std::vector<std::vector<double>> Jx_p(order_size, std::vector<double>(order_size, 0.));
        std::vector<std::vector<double>> Jy_p(order_size, std::vector<double>(order_size, 0.));


        for(auto i = 0u; i < order_size; ++i)
        {
            DSx = Sx1[i] - Sx0[i];
            DSy = Sy1[i] - Sy0[i];
        }


        for(auto i = 0u; i < order_size; ++i)
        {
            for(auto j = 0u; j < order_size; ++j)
            {
                Wx(i,j) = DSx[i] * (Sy0[j] + 0.5*DSy[j]);
                Wy(i,j) = DSy[j] * (Sx0[i] + 0.5*DSx[i]);
                Wz(i,j) = Sx0[i] * Sy0[j] + 0.5*DSx[i]*Sy0[j] + 0.5*Sx0[i]*DSy[j] + one_third*DSx[i]*DSy[j];
            }
        }


        for (auto i = 0u; i < order_size; ++i)
        {
            for(auto j = 0u; j < order_size; ++j)
            {
                auto x = xStartIndex + i - iCorr; // eg, i from -2 to 2 for 3rd order B-splines.
                auto y = yStartIndex + j - jCorr;

                Jx_p(i, j) = Jx_p(i-1, j) + cr_p * Wx(i-1, j);
                Jx(x, y) += Jx_p(i, j) ;

                Jy_p(i, j) = Jy_p(i, j-1) + cr_p * Wy(i, j-1);
                Jy(x, y) += Jy_p(i, j) ;

                Jz(x, y) += crz_p_1D2D * Wz(i, j);
            }
        }

    }
}; // END test<2> specialization




/** \brief specialization of ParticleToMesh for 3D interpolation of J */
template<>
class test<3>
{
public:
    template<typename VecField, typename Particle>
    inline void operator()(VecField& J, Particle const& partIn,
                           Particle const& partOut, double dt)
    {
        auto& Jx = J(Component::X);
        auto& Jy = J(Component::Y);
        auto& Jz = J(Component::Z);

        auto const& xStartIndex = oldPart.iCell[0];
        auto const& yStartIndex = oldPart.iCell[1];
        auto const& zStartIndex = oldPart.iCell[2];

        auto const& oldWeights    = E_weights_(partOut);
        auto const& newWeights    = E_weights_(partIn);
        auto const& [Sx0, Sy0, Sz0]      = oldWeights;
        auto const& [Sx1, Sy1, Sz1]      = newWeights;

        auto const& order_size      = Sx0.size();

        // requisite for appropriate centering
        int iCorr = order_size/2;
        if Sx0[1] > Sx0[Sx0.size()-2] { 
            iCorr -= 1
        }
        int jCorr = order_size/2;
        if Sy0[1] > Sy0[Sy0.size()-2] {
            jCorr -= 1
        }
        int kCorr = order_size/2;
        if Sz0[1] > Sz0[Sz0.size()-2] {
            kCorr -= 1
        }

        double inv_cell_volume = 1.0 / layout_->cellVolume();
        double charge_weight = inv_cell_volume * partIn.charge * partIn.weight;
            
        double cr_p = charge_weight/dt;         // current density in the evaluated dimension
        double cry_p_1D = charge_weight*partIn.v[1];   // current density in the y-direction in 1D
        double crz_p_1D2D = charge_weight*partIn.v[2]; // current density in the z-direction in 1D or 2D

        std::vector<std::vector<std::vector<double>>> Jx_p(order_size, std::vector<std::vector<double>>(order_size, std::vector<double>(order_size, 0.)));
        std::vector<std::vector<std::vector<double>>> Jy_p(order_size, std::vector<std::vector<double>>(order_size, std::vector<double>(order_size, 0.)));
        std::vector<std::vector<std::vector<double>>> Jz_p(order_size, std::vector<std::vector<double>>(order_size, std::vector<double>(order_size, 0.)));

        for(auto i = 0u; i < order_size; ++i)
        {
            DSx = Sx1[i] - Sx0[i];
            DSy = Sy1[i] - Sy0[i];
            DSz = Sz1[i] - Sz0[i];
        }
        
        for(auto i = 0u; i < order_size; ++i)
        {
            for(auto j = 0u; j < order_size; ++j)
            {
                for(auto k = 0u; k < order_size; ++k)
                {
                Wx(i,j) = DSx[i] * (Sy0[j]*Sz0[k] + 0.5*DSy[j]*Sz0[k] + 0.5*Sy0[j]*DSz[k] + one_third*DSy[j]*DSz[k]);
                Wy(i,j) = DSy[j] * (Sx0[i]*Sz0[k] + 0.5*DSx[i]*Sz0[k] + 0.5*Sx0[i]*DSz[k] + one_third*DSx[i]*DSz[k]);
                Wz(i,j) = DSz[k] * (Sx0[i]*Sy0[j] + 0.5*DSx[i]*Sy0[j] + 0.5*Sx0[i]*DSy[j] + one_third*DSx[i]*DSy[j]);
                }
            }
        }

        for (auto i = 0u; i < order_size; ++i)
        {
            for (auto j = 0u; j < order_size; ++j)
            {
                for (auto k = 0u; k < order_size; ++k)
                {
                auto x = xStartIndex + i - iCorr; // eg, i from -2 to 2 for 3rd order B-splines.
                auto y = yStartIndex + j - jCorr;
                auto z = zStartIndex + k - kCorr;

                Jx_p(i, j, k) = Jx_p(i-1, j, k) + cr_p * Wx(i-1, j, k);
                Jx(x, y, z) += Jx_p(i, j, k) ;

                Jy_p(i, j, k) = Jy_p(i, j-1, k) + cr_p * Wy(i, j-1, k);
                Jy(x, y, z) += Jy_p(i, j, k) ;

                Jz_p(i, j, k) = Jz_p(i, j-1, k) + cr_p * Wz(i, j-1, k);
                Jz(x, y, z) += Jz_p(i, j, k) ;

                }
            }
        }
    }
}; // END test<3> specialization



} //namespace PHARE::core 
#endif
