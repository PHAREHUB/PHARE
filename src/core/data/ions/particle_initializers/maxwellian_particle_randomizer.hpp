#ifndef PHARE_FLUID_PARTICLE_RANDOMIZER_HPP
#define PHARE_FLUID_PARTICLE_RANDOMIZER_HPP

#include <random>

#include <fstream>

#include "core/utilities/mpi_utils.hpp"

#if defined(PHARE_DEBUG_PARTICLE_INIT) && PHARE_DEBUG_PARTICLE_INIT == 1
#define PHARE_WITH_DEBUG_PARTICLE_INIT(...) __VA_ARGS__
#else
#define PHARE_WITH_DEBUG_PARTICLE_INIT(...)
#endif

#if !defined(PHARE_DEBUG_PARTICLE_INIT)
#define PHARE_DEBUG_PARTICLE_INIT 0
#endif


namespace PHARE::core::detail
{
inline auto seed_file_name(std::string const& patch_id)
{
    return ".phare/debug/seeds_" + patch_id + "_" + std::to_string(mpi::rank());
}

template<typename SeedSequence>
void dump_to_rank_file(std::string const& patch_id, SeedSequence& seq)
{
    std::vector<std::uint32_t> seeds(seq.size());
    seq.generate(seeds.begin(), seeds.end());
    std::ofstream seed_file{seed_file_name(patch_id), std::ios::out | std::ios::trunc};
    for (auto const& seed : seeds)
        seed_file << seed << "\n";
    seed_file << std::flush;
    seed_file.close();
}

template<typename String>
auto load_from_rank_file(String const& patch_id)
{
    auto file_contents = [&]() {
        std::ifstream seed_file{seed_file_name(patch_id), std::ios::in};
        std::stringstream buffer;
        buffer << seed_file.rdbuf();
        return buffer;
    }();

    std::vector<std::uint32_t> seeds;
    std::string line;
    while (std::getline(file_contents, line))
        if (!line.empty())
        {
            std::stringstream to_int(line);
            to_int >> seeds.emplace_back();
            if (to_int.fail())
                throw std::runtime_error("Failed to parse string as std::uint32_t");
        }

    return std::seed_seq(seeds.begin(), seeds.end());
}
} // namespace PHARE::core::detail

namespace PHARE::core
{
struct MaxwellianParticleRandomizer
{
    constexpr static bool debug_init = PHARE_DEBUG_PARTICLE_INIT;

    static auto getRNG([[maybe_unused]] std::string const& patch_id,
                       std::optional<std::size_t> const& seed)
    {
        if constexpr (debug_init == true)
            if (auto env = get_env("PHARE_RELOAD_SEEDS"); env and *env == "1")
            {
                PHARE_WITH_DEBUG_PARTICLE_INIT(                       //
                    auto seq = detail::load_from_rank_file(patch_id); //
                    return std::mt19937_64(seq));
            }

        if (!seed.has_value())
        {
            std::random_device randSeed;
            std::seed_seq seed_seq{randSeed(), randSeed(), randSeed(), randSeed(),
                                   randSeed(), randSeed(), randSeed(), randSeed()};

            if constexpr (debug_init == true)
                if (auto env = get_env("PHARE_DUMP_SEEDS"); env and *env == "1")
                {
                    PHARE_WITH_DEBUG_PARTICLE_INIT(detail::dump_to_rank_file(patch_id, seed_seq));
                }

            return std::mt19937_64(seed_seq);
        }
        return std::mt19937_64(*seed);
    }
};
} // namespace PHARE::core


#endif /* PHARE_FLUID_PARTICLE_RANDOMIZER_HPP */
