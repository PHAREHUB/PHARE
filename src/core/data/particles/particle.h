#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_H

#include "core/utilities/point/point.h"

#include <array>
#include <type_traits>




namespace PHARE
{
namespace core
{
    template<std::size_t dim>
    struct Particle
    {
    };



    template<>
    struct Particle<1>
    {
        double weight;
        double charge;

        std::array<int, 1> iCell   = {{0}};
        std::array<float, 1> delta = {{0.0f}};
        std::array<double, 3> v    = {{0., 0., 0.}};

        double Ex = 0, Ey = 0, Ez = 0;
        double Bx = 0, By = 0, Bz = 0;

        static const std::size_t dimension = 1;
    };




    template<>
    struct Particle<2>
    {
        double weight;
        double charge;

        std::array<int, 2> iCell   = {{0, 0}};
        std::array<float, 2> delta = {{0.0f, 0.0f}};
        std::array<double, 3> v    = {{0., 0., 0.}};

        double Ex = 0, Ey = 0, Ez = 0;
        double Bx = 0, By = 0, Bz = 0;

        static const std::size_t dimension = 2;
    };



    template<>
    struct Particle<3>
    {
        double weight;
        double charge;

        std::array<int, 3> iCell   = {{0, 0, 0}};
        std::array<float, 3> delta = {{0.f, 0.f, 0.f}};
        std::array<double, 3> v    = {{0., 0., 0.}};

        double Ex = 0, Ey = 0, Ez = 0;
        double Bx = 0, By = 0, Bz = 0;

        static const std::size_t dimension = 3;
    };




    template<typename Particle>
    auto cellAsPoint(Particle const& particle)
    {
        return Point<int, Particle::dimension>{particle.iCell};
    }


} // namespace core

} // namespace PHARE

#endif
