#ifndef PHARE_MAXWELLIAN_PARTICLE_INITIALIZER_BASE_HPP
#define PHARE_MAXWELLIAN_PARTICLE_INITIALIZER_BASE_HPP

#include <array>
#include <tuple>
#include <random>

#include "core/def.hpp"
#include "core/utilities/types.hpp"


namespace PHARE
{
namespace core
{

    template<typename T, typename... Args>
    void maxwellianVelocity(std::array<T, 3>& partVelocity, std::mt19937_64& generator,
                            Args const... args)
    {
        auto const& [V0, V1, V2, Vth0, Vth1, Vth2] = std::forward_as_tuple(args...);

        std::normal_distribution<> maxwellX(V0, Vth0);
        std::normal_distribution<> maxwellY(V1, Vth1);
        std::normal_distribution<> maxwellZ(V2, Vth2);

        partVelocity[0] = maxwellX(generator);
        partVelocity[1] = maxwellY(generator);
        partVelocity[2] = maxwellZ(generator);
    }



    template<typename T>
    NO_DISCARD std::array<T, 3> basisTransform(std::array<std::array<double, 3>, 3> const& basis,
                                               std::array<T, 3> const& vec)
    {
        std::array<T, 3> newVec;

        for (std::uint32_t comp = 0; comp < 3; comp++)
        {
            newVec[comp]
                = basis[0][comp] * vec[0] + basis[1][comp] * vec[1] + basis[2][comp] * vec[2];
        }

        return newVec;
    }


    template<typename T0, typename... Bargs>
    void localMagneticBasis(std::array<std::array<T0, 3>, 3>& basis, Bargs const... bargs)
    {
        std::array const B{bargs...};

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

#endif /*PHARE_MAXWELLIAN_PARTICLE_INITIALIZER_BASE_HPP*/
