#ifndef PHARE_FLUID_PARTICLE_INITIALIZER_HPP
#define PHARE_FLUID_PARTICLE_INITIALIZER_HPP

#include <memory>
#include <random>
#include <cassert>
#include <functional>

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/utilities/types.hpp"
#include "core/data/ions/particle_initializers/particle_initializer.hpp"
#include "core/data/particles/particle.hpp"
#include "initializer/data_provider.hpp"
#include "core/utilities/point/point.hpp"
#include "core/def.hpp"


namespace PHARE::core
{
void maxwellianVelocity(std::array<double, 3> const& V, std::array<double, 3> const& Vth,
                        std::mt19937_64& generator, std::array<double, 3>& partVelocity);


std::array<double, 3> basisTransform(std::array<std::array<double, 3>, 3> const& basis,
                                     std::array<double, 3> const& vec);

void localMagneticBasis(std::array<double, 3> B, std::array<std::array<double, 3>, 3>& basis);


/** @brief a MaxwellianParticleInitializer is a ParticleInitializer that loads particles from a
 * local Maxwellian distribution given density, bulk velocity and thermal velocity profiles.
 */
template<typename ParticleArray, typename GridLayout>
class MaxwellianParticleInitializer : public ParticleInitializer<ParticleArray, GridLayout>
{
public:
    static constexpr auto dimension = GridLayout::dimension;
    using InputFunction             = initializer::InitFunction<dimension>;

    MaxwellianParticleInitializer(
        InputFunction density, std::array<InputFunction, 3> const& bulkVelocity,
        std::array<InputFunction, 3> const& thermalVelocity, double const particleCharge,
        std::uint32_t const& nbrParticlesPerCell, std::optional<std::size_t> seed = {},
        Basis const basis                                 = Basis::Cartesian,
        std::array<InputFunction, 3> const& magneticField = {nullptr, nullptr, nullptr},
        double densityCutOff                              = 1e-5)
        : density_{density}
        , bulkVelocity_{bulkVelocity}
        , thermalVelocity_{thermalVelocity}
        , magneticField_{magneticField}
        , particleCharge_{particleCharge}
        , densityCutOff_{densityCutOff}
        , nbrParticlePerCell_{nbrParticlesPerCell}
        , basis_{basis}
        , rngSeed_{seed}
    {
    }


    /**
     * @brief load pacore
} // namespace PHARErticles in a ParticleArray in a domain defined by the given layout
     */
    void loadParticles(ParticleArray& particles, GridLayout const& layout) const override;


    virtual ~MaxwellianParticleInitializer() = default;


    NO_DISCARD static std::mt19937_64 getRNG(std::optional<std::size_t> const& seed)
    {
        if (!seed.has_value())
        {
            std::random_device randSeed;
            std::seed_seq seed_seq{randSeed(), randSeed(), randSeed(), randSeed(),
                                   randSeed(), randSeed(), randSeed(), randSeed()};
            return std::mt19937_64(seed_seq);
        }
        return std::mt19937_64(*seed);
    }

private:
    using Particle = typename ParticleArray::value_type;
    InputFunction density_;
    std::array<InputFunction, 3> bulkVelocity_;
    std::array<InputFunction, 3> thermalVelocity_;
    std::array<InputFunction, 3> magneticField_;

    double particleCharge_;
    double densityCutOff_ = 1e-5;
    std::uint32_t nbrParticlePerCell_;
    Basis basis_;
    std::optional<std::size_t> rngSeed_;
};


class MaxwellianInitFunctions
{
public:
    template<typename Function, typename FunctionArray, typename Basis, typename... Coords>
    MaxwellianInitFunctions(Function& density, FunctionArray& bulkVelocity,
                            FunctionArray& thermalVelocity, FunctionArray& magneticField,
                            Basis const& basis, Coords const&... coords)
        : _n{density(coords...)}

    {
        static_assert(sizeof...(coords) <= 3, "can only provide up to 3 coordinates");
        for (std::uint32_t i = 0; i < 3; i++)
        {
            _V[i]   = bulkVelocity[i](coords...);
            _Vth[i] = thermalVelocity[i](coords...);
        }
        if (basis == Basis::Magnetic)
            for (std::uint32_t i = 0; i < 3; i++)
                _B[i] = magneticField[i](coords...);
    }

    NO_DISCARD std::array<double const*, 3> B() const { return ptrs(_B); }

    NO_DISCARD auto operator()() const { return std::make_tuple(_n->data(), ptrs(_V), ptrs(_Vth)); }

private:
    NO_DISCARD std::array<double const*, 3>
    ptrs(std::array<std::shared_ptr<PHARE::core::Span<double>>, 3> const& v) const
    {
        return {v[0]->data(), v[1]->data(), v[2]->data()};
    }

    std::shared_ptr<PHARE::core::Span<double>> const _n, _ppc;
    std::array<std::shared_ptr<PHARE::core::Span<double>>, 3> _B, _V, _Vth;
};




template<typename ParticleArray, typename GridLayout>
void MaxwellianParticleInitializer<ParticleArray, GridLayout>::loadParticles(
    ParticleArray& particles, GridLayout const& layout) const
{
    auto const icell = [](auto const idx, auto const& indices) {
        return for_N_make_array<dimension>([&](auto i) { return std::get<i>(indices[idx]); });
    };

    auto const deltas = [](auto& pos, auto& gen) {
        return for_N_make_array<dimension>([&](auto) { return pos(gen); });
    };

    // in the following two calls,
    // primal indexes are given here because that's what cellCenteredCoordinates takes

    // indices = std::vector<std::tuple<int, per dim>>
    auto const ndCellIndices = layout.indicis(layout.AMRBox());

    // coords = std::tuple<std::vector<double>,  per dim>
    auto cellCoords = layout.indexesToCoordVectors(
        ndCellIndices, QtyCentering::primal, [](auto const& gridLayout, auto const&... indexes) {
            return gridLayout.cellCenteredCoordinates(indexes...);
        });

    auto const fns = std::make_from_tuple<MaxwellianInitFunctions>(std::tuple_cat(
        std::forward_as_tuple(density_, bulkVelocity_, thermalVelocity_, magneticField_, basis_),
        cellCoords));

    auto const [n, V, Vth] = fns();
    auto randGen           = getRNG(rngSeed_);
    ParticleDeltaDistribution<double> deltaDistrib;

    for (std::size_t flatCellIdx = 0; flatCellIdx < ndCellIndices.size(); ++flatCellIdx)
    {
        if (n[flatCellIdx] < densityCutOff_)
            continue;

        auto const cellWeight = n[flatCellIdx] / nbrParticlePerCell_;
        auto const iCell      = icell(flatCellIdx, ndCellIndices);
        std::array<double, 3> particleVelocity;
        std::array<std::array<double, 3>, 3> basis;

        if (basis_ == Basis::Magnetic)
        {
            auto const B = fns.B();
            localMagneticBasis({B[0][flatCellIdx], B[1][flatCellIdx], B[2][flatCellIdx]}, basis);
        }

        for (std::uint32_t ipart = 0; ipart < nbrParticlePerCell_; ++ipart)
        {
            maxwellianVelocity({V[0][flatCellIdx], V[1][flatCellIdx], V[2][flatCellIdx]},
                               {Vth[0][flatCellIdx], Vth[1][flatCellIdx], Vth[2][flatCellIdx]}, //
                               randGen, particleVelocity);

            if (basis_ == Basis::Magnetic)
                particleVelocity = basisTransform(basis, particleVelocity);

            particles.emplace_back(Particle{cellWeight, particleCharge_, iCell,
                                            deltas(deltaDistrib, randGen), particleVelocity});
        }
    }
}

} // namespace PHARE::core


#endif
