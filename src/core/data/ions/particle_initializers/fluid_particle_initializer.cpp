
//#include "data/ions/particle_initializers/fluid_particle_initializer.h"
#include <array>
#include <random>

#include "utilities/types.h"


namespace PHARE
{
namespace core
{
    void maxwellianVelocity(std::array<double, 3> V, std::array<double, 3> Vth,
                            std::mt19937_64 generator, std::array<double, 3>& partVelocity)
    {
        std::normal_distribution<> maxwellX(V[0], Vth[0]);
        std::normal_distribution<> maxwellY(V[1], Vth[1]);
        std::normal_distribution<> maxwellZ(V[2], Vth[2]);

        partVelocity[0] = maxwellX(generator);
        partVelocity[1] = maxwellY(generator);
        partVelocity[2] = maxwellZ(generator);
    }




    std::array<double, 3> basisTransform(const std::array<std::array<double, 3>, 3> basis,
                                         std::array<double, 3> vec)
    {
        std::array<double, 3> newVec;

        for (uint32 comp = 0; comp < 3; comp++)
        {
            newVec[comp]
                = basis[0][comp] * vec[0] + basis[1][comp] * vec[1] + basis[2][comp] * vec[2];
        }

        return newVec;
    }



    void localMagneticBasis(std::array<double, 3> B, std::array<std::array<double, 3>, 3>& basis)
    {
        auto b2 = norm(B);

        if (b2 < 1e-8)
        {
            basis[0][0] = 1.0;
            basis[0][1] = 0.0;
            basis[0][2] = 0.0;

            basis[1][0] = 0.0;
            basis[1][1] = 1.0;
            basis[1][2] = 0.0;

            basis[2][0] = 0.0;
            basis[2][1] = 0.0;
            basis[2][2] = 1.0;
        }
        else
        {
            // first basis vector is aligned with B
            basis[0][0] = B[0] / b2;
            basis[0][1] = B[1] / b2;
            basis[0][2] = B[2] / b2;

            // second vector is (1,1,1) x B
            basis[1][0] = B[2] - B[1];
            basis[1][1] = B[0] - B[2];
            basis[1][2] = B[1] - B[0];

            auto vecNorm = norm(basis[1]);
            basis[1][0] /= vecNorm;
            basis[1][1] /= vecNorm;
            basis[1][2] /= vecNorm;

            // last vector is just the cross product of the first two vectors
            basis[2][0] = basis[0][1] * basis[1][2] - basis[0][2] * basis[1][1];
            basis[2][1] = basis[0][2] * basis[1][0] - basis[0][0] * basis[1][2];
            basis[2][2] = basis[0][0] * basis[1][1] - basis[0][1] * basis[1][0];
        }
    }


} // namespace core
} // namespace PHARE
