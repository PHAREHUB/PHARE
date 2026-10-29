#ifndef PHARE_CORE_UTILITIES_RANDOM_HPP
#define PHARE_CORE_UTILITIES_RANDOM_HPP

#include <array>
#include <cstdint>
#include <limits>
#include <string_view>

/* Tools to give every cell its own reproducible random stream.
 *
 * Seeded particle loading must not depend on how a level is cut into patches or spread over
 * MPI ranks. So the random numbers of a cell are derived from what identifies the cell, never
 * from the patch: (user seed, population, level, global AMR cell index) is hashed into a
 * 64-bit "cell key", and that key seeds a small generator used for that cell only.
 *
 *   seed ──┐
 *   pop  ──┼─ hash_combine chain ─> cell key ─> SplitMix64 ─> deltas, velocities of the cell
 *   level ─┤
 *   cell ──┘
 *
 * Everything here is plain 64-bit integer arithmetic, so it gives the same bits on every
 * platform and compiler (unsigned overflow wraps modulo 2^64 by definition in C++).
 */

namespace PHARE::core
{
/** @brief 2^64 / golden ratio, rounded to an odd integer ("golden gamma").
 *
 * Adding it repeatedly modulo 2^64 visits all 2^64 values before coming back (any odd step
 * does), and the golden ratio spreads consecutive values far apart. Used as the SplitMix64
 * increment and as an offset in hash_combine.
 */
inline constexpr std::uint64_t golden_gamma = 0x9e3779b97f4a7c15ULL;


/** @brief SplitMix64 finalizer: scrambles a 64-bit integer into a random-looking one.
 *
 * Goal: inputs that are close (cell 3 and cell 4, or a counter n and n+1) must give outputs
 * that look unrelated. Flipping any single input bit flips about half of the output bits
 * ("avalanche").
 *
 * How: two rounds of "xor with a shifted copy, then multiply by a large odd constant".
 * - z ^ (z >> k) folds the high bits down into the low bits;
 * - multiplying by an odd constant pushes the low bits up into the high bits.
 * Alternating both lets every input bit influence every output bit.
 *
 * Bijective: each step can be undone (an xor-shift is invertible, an odd number is invertible
 * modulo 2^64), so two different inputs never give the same output. Note mix(0) == 0, which is
 * why callers add golden_gamma before mixing.
 *
 * The shifts (30, 27, 31) and constants are D. Stafford's "Mix13" variant of the MurmurHash3
 * 64-bit finalizer, as used in Java's SplittableRandom and in S. Vigna's splitmix64.c.
 */
constexpr std::uint64_t splitmix64_mix(std::uint64_t z)
{
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}


/** @brief folds value into key, giving a new key that depends on both.
 *
 * Used to build a key from several integers one after the other:
 *   key = hash_combine(hash_combine(hash_combine(seed, a), b), c)
 * - order matters: combining (3, 7) and (7, 3) gives different keys, so cell (3,7) and cell
 *   (7,3) get different streams;
 * - for a fixed key, two different values always give two different results (each step is
 *   bijective: + golden_gamma, mix, xor key, mix). So two cells that only differ in their last
 *   index can never share a key. Any other coincidence has a probability of about 2^-64.
 * - value is mixed before the xor so that small values (0, 1, 2...) do not just flip a few low
 *   bits of the key; the final mix spreads the result again.
 */
constexpr std::uint64_t hash_combine(std::uint64_t const key, std::uint64_t const value)
{
    return splitmix64_mix(key ^ splitmix64_mix(value + golden_gamma));
}


/** @brief FNV-1a 64-bit hash of a string, to turn a population name into a number.
 *
 * std::hash cannot be used: its value is only guaranteed within one program execution and
 * differs between standard libraries, so the same seed would give different particles on two
 * machines. FNV-1a is fixed by its definition: start from the "offset basis" and, for each
 * byte, xor it in then multiply by the "FNV prime". It is not a strong hash, which is fine:
 * its output goes through hash_combine, which does the scrambling.
 */
constexpr std::uint64_t fnv1a64(std::string_view const str)
{
    std::uint64_t hash = 0xcbf29ce484222325ULL; // FNV-1a 64-bit offset basis
    for (char const c : str)
    {
        hash ^= static_cast<std::uint8_t>(c);
        hash *= 0x100000001b3ULL; // FNV-1a 64-bit prime
    }
    return hash;
}


/** @brief key of the random stream of one cell, a pure function of (seed, population, level,
 * global AMR cell index): it does not depend on how the level is cut into patches.
 *
 * - population: two populations with the same seed must not get the same particles;
 * - level: the same AMR index exists on every level but is a different physical cell;
 * - cell: global AMR index, the same whichever patch or rank holds the cell. Negative indices
 *   are sign-extended to 64 bits, so -2 and 2 give different keys.
 *
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


/** @brief SplitMix64 generator, a UniformRandomBitGenerator with an 8-byte state.
 *
 * How it works: the state is a counter that advances by golden_gamma at each call, and the
 * output is that counter scrambled by splitmix64_mix. The n-th number of a stream started from
 * key k is mix(k + n * golden_gamma): the counter guarantees a period of 2^64, the mix makes
 * consecutive outputs look independent.
 *
 * Why it, and not std::mt19937_64: particles are loaded with one generator per cell, so
 * building a generator must be cheap. SplitMix64 starts from a single 64-bit value, whereas
 * std::mt19937_64 carries ~2.5 kB of state and a costly seeding. SplitMix64 passes TestU01
 * BigCrush, the standard battery of statistical tests for random generators.
 *
 * Streams of two cells overlap only if their keys differ by fewer than m steps of
 * golden_gamma, m being the numbers drawn per cell. For N cells the chance is about
 * N^2 m / 2^64: ~1e-2 for 1.7e7 cells (256^3) drawing 600 numbers each, and even then only
 * one pair of cells would be affected.
 *
 * result_type, min() and max() are what std:: distributions (std::normal_distribution...)
 * require from a generator.
 *
 * References:
 * - G. L. Steele, D. Lea, C. H. Flood, "Fast splittable pseudorandom number generators",
 *   OOPSLA 2014, doi:10.1145/2660193.2660195 (the generator behind Java's SplittableRandom)
 * - S. Vigna, reference implementation https://prng.di.unimi.it/splitmix64.c, reproduced here
 *   (state += 0x9e3779b97f4a7c15, output = splitmix64_mix(state))
 * - P. L'Ecuyer, R. Simard, "TestU01: a C library for empirical testing of random number
 *   generators", ACM Trans. Math. Softw. 33(4), 2007 (BigCrush)
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
        state_ += golden_gamma;
        return splitmix64_mix(state_);
    }

private:
    std::uint64_t state_;
};

} // namespace PHARE::core

#endif
