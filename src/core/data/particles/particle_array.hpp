#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_HPP

#include "core/data/particles/particle.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_detail.hpp"

#include "core/utilities/box/box.hpp"
#include "core/utilities/equality.hpp"

#include <utility>
#include <sstream>
#include <stdexcept>

namespace PHARE::core
{

template<auto opts /* defaulted in details header */>
class ParticleArray : public ResolvedParticleArray_t<opts>
{
    using This      = ParticleArray<opts>;
    using internals = ParticleArrayResolver<opts>;

public:
    using Super      = ResolvedParticleArray_t<opts>;
    using value_type = ParticleDefaults<opts.dim>::Particle_t;
    using view_t     = ParticleArray<opts.with_storage(StorageMode::SPAN)>;

    auto static constexpr options      = opts;
    auto static constexpr dimension    = opts.dim;
    auto static constexpr alloc_mode   = opts.alloc_mode;
    auto static constexpr layout_mode  = opts.layout_mode;
    auto static constexpr storage_mode = opts.storage_mode;
    auto static constexpr type_id      = internals::type_id;

    std::string static id() { return std::string{type_id}; }

    ParticleArray(ParticleArray&& that)
        : Super{std::forward<Super>(that)}
    {
    }
    ParticleArray(ParticleArray const& that)
        : Super{that}
    {
    }

    ParticleArray& operator=(ParticleArray&& that)
    {
        super() = std::move(that.super());
        return *this;
    }
    ParticleArray& operator=(ParticleArray const& that)
    {
        super() = that.super();
        return *this;
    }

    template<typename... Args>
    ParticleArray(Args&&... args)
        requires(self_excluding_constructible<This, Super, Args...>())
        : Super{std::forward<Args>(args)...}
    {
    }

    auto view() { return view_t{*this}; }
    auto view() const { return view_t{*this}; }

    auto view(std::size_t i) // to take only i particles and ignore the rest
    {
        view_t v{*this};
        v.super().resize(i);
        return v;
    }

    auto view(std::size_t const start, std::size_t const size)
    {
        return view_t{*this, start, size};
    }

    auto operator*() { return view(); }
    auto operator*() const { return view(); }

    auto begin()
    {
        if constexpr (requires { super().begin(); })
            return super().begin();
        else
            static_assert(dependent_false_v<This>, "iteration not supported for this layout");
    }
    auto begin() const
    {
        if constexpr (requires { super().begin(); })
            return super().begin();
        else
            static_assert(dependent_false_v<This>, "iteration not supported for this layout");
    }
    auto end()
    {
        if constexpr (requires { super().end(); })
            return super().end();
        else
            static_assert(dependent_false_v<This>, "iteration not supported for this layout");
    }
    auto end() const
    {
        if constexpr (requires { super().end(); })
            return super().end();
        else
            static_assert(dependent_false_v<This>, "iteration not supported for this layout");
    }

    Super& super() { return *this; }
    Super const& super() const { return *this; }

    template<auto _opts>
    friend std::ostream& operator<<(std::ostream& out, ParticleArray<_opts> const&);
};

template<std::size_t dim>
using AoSParticleArray = ParticleArray<ParticleArrayOptions{dim, LayoutMode::AoS}>;


template<std::size_t dim>
using AoSMappedParticleArray = ParticleArray<ParticleArrayOptions{dim, LayoutMode::AoSMapped}>;

// internal only - see LayoutMode::SoA
template<std::size_t dim>
using SoAParticleArray = ParticleArray<ParticleArrayOptions{dim, LayoutMode::SoA}>;

template<auto opts>
std::ostream& operator<<(std::ostream& out, ParticleArray<opts> const& arr)
{
    for (auto const& p : arr)
        out << p.copy();
    return out;
}




template<auto opts>
void empty(ParticleArray<opts>& array)
{
    array.clear();
}


template<auto opts>
void swap(ParticleArray<opts>& array1, ParticleArray<opts>& array2)
{
    array1.swap(array2);
}

template<typename P0, typename P1>
EqualityReport particle_compare(P0 const& p0, P1 const& p1, std::size_t const i = 0,
                                double const atol = 1e-15)
{
    std::string idx = std::to_string(i);
    if (p0.iCell() != p1.iCell())
        return EqualityReport{false, "icell mismatch at index: " + idx, i};

    if (!float_equals(p0.v(), p1.v(), atol))
        return EqualityReport{false, "v mismatch at index: " + idx, i};
    if (!float_equals(p0.delta(), p1.delta(), atol))
        return EqualityReport{false, "delta mismatch at index: " + idx, i};


    return EqualityReport{true};
}


template<typename P0, typename P1>
EqualityReport particles_equals(P0 const& ref, P1 const& cmp, double const atol = 1e-15)
{
    if (ref.size() != cmp.size())
        return EqualityReport{false, "different sizes: " + std::to_string(ref.size()) + " vs "
                                         + std::to_string(cmp.size())};

    auto rit      = ref.begin();
    auto cit      = cmp.begin();
    std::size_t i = 0;

    for (; rit != ref.end(); ++rit, ++cit, ++i)
        if (auto const eq = particle_compare(*rit, *cit, i, atol); !eq)
            return eq;

    return EqualityReport{true};
}

template<auto o>
EqualityReport operator==(ParticleArray<o> const& p0, ParticleArray<o> const& p1)
{
    auto report = particles_equals(p0, p1);
    if (!report)
    {
        PHARE_LOG_LINE_STR(p0[report.idx].copy());
        PHARE_LOG_LINE_STR(p1[report.idx].copy());
    }
    return report;
}

template<auto o>
EqualityReport operator==(ParticleArray<o> const& p0, std::vector<Particle<o>> const& p1)
{
    return particles_equals(p0, p1);
}



template<typename ParticleArray_t>
auto constexpr base_layout_type()
{
    return LayoutMode::AoS;
}


template<typename ParticleArray_t>
void check_level_ghost_particles(ParticleArray_t const& particles)
    requires(ParticleArray_t::layout_mode == LayoutMode::AoSPCTS)
{
    // level ghost particles are expected to be duplicated per tile: any cell reachable
    // (via ghost halo) from more than one tile must hold the same particle count in
    // every tile that reaches it, or some tile is missing contributions.

    auto const& tiles = particles();

    auto const count_at = [](auto const& ps, auto const& cell) { return ps(cell).size(); };

    for (auto const& bix : particles.ghost_box())
    {
        if (isIn(bix, particles.box()))
            continue;

        for (std::size_t i = 0; i < tiles.size(); ++i)
        {
            auto const& tile_0 = tiles[i];

            if (not isIn(bix, tile_0().ghost_box()))
                continue;

            for (std::size_t j = i + 1; j < tiles.size(); ++j)
            {
                auto const& tile_1 = tiles[j];

                if (not isIn(bix, tile_1().ghost_box()))
                    continue;

                auto const lcl_0 = (bix - tile_0().ghost_box().lower).as_unsigned();
                auto const lcl_1 = (bix - tile_1().ghost_box().lower).as_unsigned();
                auto const na    = count_at(tile_0(), lcl_0);
                auto const nb    = count_at(tile_1(), lcl_1);

                if (na != nb)
                {
                    std::ostringstream oss;
                    oss << "check_level_ghost_particles: tile mismatch at cell " << bix << ": "
                        << na << " vs " << nb;
                    throw std::runtime_error(oss.str());
                }
            }
        }
    }
}


template<typename ParticleArray_t>
void check_level_ghost_particles(ParticleArray_t const& particles)
{
    // fallthrough
}


} // namespace PHARE::core


#endif
