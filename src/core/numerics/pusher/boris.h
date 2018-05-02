#ifndef PHARE_CORE_PUSHER_BORIS_H
#define PHARE_CORE_PUSHER_BORIS_H

#include <array>
#include <cmath>
#include <cstddef>


template<std::size_t dim>
class BorisPusher
{
public:
    template<typename ParticleRange, typename Electromag, typename Interpolator,
             typename ParticleSelector, typename BoundaryCondition>
    auto move(ParticleRange rangeIn, ParticleRange rangeOut, Electromag const& em, double mass,
              Interpolator& interpolator, ParticleSelector const& selector, BoundaryCondition&& bc)
    {
        auto firstLeaving
            = pushStep_(rangeIn.begin(), rangeIn.end(), rangeOut.begin(), rangeOut.end(), selector);

        auto newEnd = bc.applyOutgoingParticleBC(firstLeaving, rangeOut.end());
        interpolator(rangeOut.begin(), rangeOut.end(), em);
    }



    template<typename ParticleRange, typename Electromag, typename Interpolator,
             typename ParticleSelector>
    auto move(ParticleRange rangeIn, ParticleRange rangeOut, Electromag const& EM, double mass,
              Interpolator& interpolator, ParticleSelector const& selector)
    {
        auto pivot
            = pushStep_(rangeIn.begin(), rangeIn.end(), rangeOut.begin(), rangeOut.end(), selector);
    }



    void setMeshAndTimeStep(std::array<double, dim> ms, double ts)
    {
        std::transform(std::begin(ms), std::end(ms), std::begin(halfDtOverDl_),
                       [ts](double& x) { return 0.5 * x / ts; });
        dt_ = ts;
    }


private:
    /**
     *
     */
    template<typename ParticleIter>
    void advancePosition_(ParticleIter partIn, ParticleIter partOut)
    {
        // push the particle
        for (std::size_t iDim = 0; iDim < dim; ++iDim)
        {
            float delta
                = partIn.delta[iDim] + static_cast<float>(halfDtOverDl_[iDim] * partOut.v[iDim]);

            float iCell         = std::floor(delta);
            partIn.delta[iDim]  = delta - iCell;
            partOut.iCell[iDim] = static_cast<int>(iCell + partIn.iCell[iDim]);
        }
    }



    template<typename ParticleIterIn, typename ParticleIterOut, typename ParticleSelector>
    ParticleIterOut pushStep_(ParticleIterIn firstIn, ParticleIterIn lastIn,
                              ParticleIterOut firstOut, ParticleIterOut lastOut,
                              ParticleSelector const& selector)
    {
        auto pivot = lastOut;
        --pivot;
        auto pivotNext = lastOut;

        ParticleIterOut currentOut{firstOut};

        for (auto& currentIn = firstIn; currentIn != lastIn; ++currentIn)
        {
            // push the particle
            advancePosition_(*currentIn, *currentOut);

            if (!selector(*currentOut))
            {
                // if the particle satisfies the predicate
                // swap it with the pivot
                // and decrement the pivot

                std::swap(*currentOut, *pivot);
                --pivotNext;
                --pivot;
            }
            else
            {
                // we advance the output iterator
                // only if currentOut has not been
                // swapped
                ++currentOut;
            }
        }

        // now all particles have been pushed
        // those not satisfying the predicate after the push
        // are found in [pivotNext:end[
        // those for which pred is true are in [firstOut,pivotNext[
        return pivotNext;
    }




    template<typename ParticleIterIn, typename ParticleIterOut>
    void pushVelocity_(ParticleIterIn firstIn, ParticleIterIn lastIn, ParticleIterOut firstOut,
                       ParticleIterOut lastOut, double mass)
    {
        double dto2m = 0.5 * dt_ / mass;

        ParticleIterOut currentOut = firstOut;

        for (auto& currentIn = firstIn; currentIn != lastIn; ++currentIn, ++currentOut)
        {
            double coef1 = currentIn->charge * dto2m;

            // We now apply the 3 steps of the BORIS PUSHER

            // 1st half push of the electric field
            double velx1 = currentIn->v[0] + coef1 * currentIn->Ex;
            double vely1 = currentIn->v[1] + coef1 * currentIn->Ey;
            double velz1 = currentIn->v[2] + coef1 * currentIn->Ez;


            // preparing variables for magnetic rotation
            double const rx = coef1 * currentIn->Bx;
            double const ry = coef1 * currentIn->By;
            double const rz = coef1 * currentIn->Bz;

            double const rx2  = rx * rx;
            double const ry2  = ry * ry;
            double const rz2  = rz * rz;
            double const rxry = rx * ry;
            double const rxrz = rx * rz;
            double const ryrz = ry * rz;

            double const invDet = 1. / (1. + rx2 + ry2 + rz2);

            // preparing rotation matrix due to the magnetic field
            // m = invDet*(I + r*r - r x I) - I where x denotes the cross product
            double const mxx = 1. + rx2 - ry2 - rz2;
            double const mxy = 2. * (rxry + rz);
            double const mxz = 2. * (rxrz - ry);

            double const myx = 2. * (rxry - rz);
            double const myy = 1. + ry2 - rx2 - rz2;
            double const myz = 2. * (ryrz + rx);

            double const mzx = 2. * (rxrz + ry);
            double const mzy = 2. * (ryrz - rx);
            double const mzz = 1. + rz2 - rx2 - ry2;

            // magnetic rotation
            double const velx2 = (mxx * velx1 + mxy * vely1 + mxz * velz1) * invDet;
            double const vely2 = (myx * velx1 + myy * vely1 + myz * velz1) * invDet;
            double const velz2 = (mzx * velx1 + mzy * vely1 + mzz * velz1) * invDet;


            // 2nd half push of the electric field
            velx1 = velx2 + coef1 * currentIn->Ex;
            vely1 = vely2 + coef1 * currentIn->Ey;
            velz1 = velz2 + coef1 * currentIn->Ez;

            // Update particle velocity
            currentOut->v[0] = velx1;
            currentOut->v[1] = vely1;
            currentOut->v[2] = velz1;
        }
    }




    std::array<double, dim> halfDtOverDl_;
    double dt_;
};


#endif
