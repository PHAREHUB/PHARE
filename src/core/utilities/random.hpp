#ifndef PHARE_CORE_UTILITIES_RANDOM_HPP
#define PHARE_CORE_UTILITIES_RANDOM_HPP

#include <array>
#include <cstdint>
#include <limits>
#include <string_view>


namespace PHARE::core
{
/** @brief splitmix64 finalizer: a bijective 64-bit mixing function */
constexpr std::uint64_t splitmix64_mix(std::uint64_t z)
{
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}


/** @brief combines value into key, bijective in value for a fixed key */
constexpr std::uint64_t hash_combine(std::uint64_t const key, std::uint64_t const value)
{
    return splitmix64_mix(key ^ splitmix64_mix(value + 0x9e3779b97f4a7c15ULL));
}


/** @brief FNV-1a 64-bit string hash, stable across platforms unlike std::hash */
constexpr std::uint64_t fnv1a64(std::string_view const str)
{
    std::uint64_t hash = 0xcbf29ce484222325ULL;
    for (char const c : str)
    {
        hash ^= static_cast<std::uint8_t>(c);
        hash *= 0x100000001b3ULL;
    }
    return hash;
}


/** @brief key of the random stream of one cell, a pure function of (seed, population, level,
 * global AMR cell index): it does not depend on how the level is cut into patches.
 * Combining order: seed, population, level, cell[0], cell[1], cell[2].
 * For a fixed prefix the last combine is bijective, so cells differing in their last index
 * never share a key; other collisions have probability ~2^-64.
 */
template<std::size_t dim>
constexpr std::uint64_t cell_key(std::uint64_t const seed, std::uint64_t const populationKey,
                                 int const level, std::array<int, dim> const& amrCell)
{
    auto key = hash_combine(hash_combine(seed, populationKey), static_cast<std::uint64_t>(level));
    for (auto const index : amrCell)
        key = hash_combine(key, static_cast<std::uint64_t>(static_cast<std::int64_t>(index)));
    return key;
}


/** @brief SplitMix64 generator (Steele, Lea & Flood 2014; Vigna's reference version),
 * a UniformRandomBitGenerator with an 8-byte state, cheap enough to build one per cell.
 */
class SplitMix64
{
public:
    using result_type = std::uint64_t;

    constexpr explicit SplitMix64(std::uint64_t const seed)
        : state_{seed}
    {
    }

    static constexpr result_type min() { return std::numeric_limits<result_type>::min(); }
    static constexpr result_type max() { return std::numeric_limits<result_type>::max(); }

    constexpr result_type operator()()
    {
        state_ += 0x9e3779b97f4a7c15ULL;
        return splitmix64_mix(state_);
    }

private:
    std::uint64_t state_;
};

} // namespace PHARE::core

#endif
