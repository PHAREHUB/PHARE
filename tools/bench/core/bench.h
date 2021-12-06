#ifndef PHARE_BENCH_CORE_BENCH_H
#define PHARE_BENCH_CORE_BENCH_H


#include "benchmark/benchmark.h"

#include "phare_core.h"
#include "core/utilities/types.h"



namespace PHARE::core::bench
{
template<std::size_t dim>
using Field = PHARE::core::Field<PHARE::core::NdArrayVector<dim>,
                                 typename PHARE::core::HybridQuantity::Scalar>;
template<std::size_t dim>
using VecField
    = PHARE::core::VecField<PHARE::core::NdArrayVector<dim>, typename PHARE::core::HybridQuantity>;

// clang-format off
template<std::size_t dim>
PHARE::core::Particle<dim> particle(int icell = 15)
{
    return {
        /*.weight = */ .00001,
        /*.charge = */ .00001,
        /*.iCell  = */ PHARE::core::ConstArray<int, dim>(icell),
        /*.delta  = */ PHARE::core::ConstArray<double, dim>(.5),
        /*.v      = */ {{.00001, .00001, .00001}},
    };
}
// clang-format on

template<std::size_t dim>
PHARE::core::ParticleArray<dim> make_particles(std::size_t ppc)
{
    return {ppc, particle<dim>()};
}

template<typename Particles, typename Point>
void disperse(Particles& particles, Point lo, Point up, std::optional<int> seed = std::nullopt)
{
    auto gen = [&]() {
        if (!seed.has_value())
        {
            std::random_device rd;
            std::seed_seq seed_seq{rd(), rd(), rd(), rd(), rd(), rd(), rd()};
            return std::mt19937_64(seed_seq);
        }
        return std::mt19937_64(*seed);
    }();
    for (std::size_t i = 0; i < Particles::dimension; i++)
    {
        std::uniform_int_distribution<> distrib(lo[i], up[i]);
        for (auto& particle : particles)
            particle.iCell[i] = distrib(gen);
    }
}
template<typename Particles>
void disperse(Particles& particles, std::size_t lo, std::size_t up,
              std::optional<int> seed = std::nullopt)
{
    auto constexpr static dim = Particles::dimension;

    disperse(particles, core::ConstArray<int, dim>(lo), core::ConstArray<int, dim>(up), seed);
}

template<std::size_t dim, typename Box>
auto make_particles(std::size_t ppc, Box disperse_in, std::optional<int> seed = std::nullopt)
{
    auto particles = make_particles<dim>(ppc * disperse_in.size());
    disperse(particles, disperse_in.lower, disperse_in.upper, seed);
    return particles;
}


template<typename GridLayout, typename Quantity, std::size_t dim = GridLayout::dimension>
Field<dim> field(std::string key, Quantity type, GridLayout const& layout)
{
    Field<dim> feeld{key, type, layout.allocSize(type)};
    std::fill(feeld.begin(), feeld.end(), 1);
    return feeld;
}


template<typename GridLayout, std::size_t dim = GridLayout::dimension>
auto rho(GridLayout const& layout)
{
    return field("rho", HybridQuantity::Scalar::rho, layout);
}


template<typename GridLayout>
class _VF_ : public VecField<GridLayout::dimension>
{
public:
    using Super = VecField<GridLayout::dimension>;

    _VF_(_VF_ const&) = delete;
    _VF_(_VF_&&)      = delete;
    auto& operator=(_VF_ const&) = delete;
    auto& operator=(_VF_&&) = delete;

    _VF_(GridLayout const& layout, HybridQuantity::Vector v_id, std::string id,
         std::array<HybridQuantity::Scalar, 3> quantities)
        : Super{id, v_id}
        , xyz{field(id + "x", quantities[0], layout), field(id + "y", quantities[1], layout),
              field(id + "z", quantities[2], layout)}
    {
        Super::setBuffer(id + "_x", &xyz[0]);
        Super::setBuffer(id + "_y", &xyz[1]);
        Super::setBuffer(id + "_z", &xyz[2]);
    }


    _VF_(GridLayout const& layout, HybridQuantity::Vector v_id, std::string id)
        : _VF_{layout,
               v_id,
               id,
               {HybridQuantity::Scalar::Vx, HybridQuantity::Scalar::Vy, HybridQuantity::Scalar::Vz}}
    {
    }

private:
    std::array<Field<GridLayout::dimension>, 3> xyz;
};

template<typename GridLayout>
class Flux : public _VF_<GridLayout>
{
public:
    using Super = _VF_<GridLayout>;

    Flux(GridLayout const& layout)
        : Super{layout, HybridQuantity::Vector::V, "F"}
    {
    }
};

template<typename GridLayout>
class J : public _VF_<GridLayout>
{
public:
    using Super = _VF_<GridLayout>;

    J(GridLayout const& layout)
        : Super{layout, HybridQuantity::Vector::J, "J"}
    {
    }
};


template<typename GridLayout>
struct Electromag
{
    Electromag(Electromag const&) = delete;

    Electromag(GridLayout const& layout)
        : E{layout,
            HybridQuantity::Vector::E,
            "E",
            {HybridQuantity::Scalar::Ex, HybridQuantity::Scalar::Ey, HybridQuantity::Scalar::Ez}}
        , B{layout,
            HybridQuantity::Vector::B,
            "B",
            {HybridQuantity::Scalar::Bx, HybridQuantity::Scalar::By, HybridQuantity::Scalar::Bz}}
    {
    }

    auto operator()() { return std::forward_as_tuple(E, B); }

    _VF_<GridLayout> E, B;
};



} // namespace PHARE::core::bench

#endif /*PHARE_BENCH_CORE_BENCH_H*/
