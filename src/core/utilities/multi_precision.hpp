#ifndef PHARE_CORE_UTILITIES_MULTI_PRECISION_HPP
#define PHARE_CORE_UTILITIES_MULTI_PRECISION_HPP


#include "core/def.hpp"
#include "core/utilities/types.hpp"


#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <ostream>


// Reduced precision storage types.
//  Values are decoded to double on read and encoded on assignment, so all arithmetic
//  is done in double and only the stored representation loses precision.
//  Storage is a byte array (alignof == 1), so Bytes need not be a power of two.
//  Encoding rounds to nearest, and encode(decode(x)) == x bitwise.

namespace PHARE::core
{
static_assert(std::endian::native == std::endian::little,
              "multi precision types assume little endian");


// Bytes bytes to/from the low bytes of a uint64, in power of two chunks.
//  A single memcpy of e.g. 6 bytes into a uint64 goes through the stack as two narrow
//  stores and one wide load, which stalls on store forwarding on every decode.
template<std::size_t Bytes>
NO_DISCARD inline std::uint64_t load_low_bytes(std::uint8_t const* const src)
{
    static_assert(Bytes >= 1 and Bytes <= 8);
    std::uint64_t u  = 0;
    std::size_t off  = 0;
    auto const chunk = [&]<typename U>() {
        U c;
        std::memcpy(&c, src + off, sizeof(U));
        u |= static_cast<std::uint64_t>(c) << (8 * off);
        off += sizeof(U);
    };
    if constexpr (Bytes & 8)
        chunk.template operator()<std::uint64_t>();
    if constexpr (Bytes & 4)
        chunk.template operator()<std::uint32_t>();
    if constexpr (Bytes & 2)
        chunk.template operator()<std::uint16_t>();
    if constexpr (Bytes & 1)
        chunk.template operator()<std::uint8_t>();
    return u;
}

template<std::size_t Bytes>
inline void store_low_bytes(std::uint8_t* const dst, std::uint64_t const u)
{
    static_assert(Bytes >= 1 and Bytes <= 8);
    std::size_t off  = 0;
    auto const chunk = [&]<typename U>() {
        U const c = static_cast<U>(u >> (8 * off));
        std::memcpy(dst + off, &c, sizeof(U));
        off += sizeof(U);
    };
    if constexpr (Bytes & 8)
        chunk.template operator()<std::uint64_t>();
    if constexpr (Bytes & 4)
        chunk.template operator()<std::uint32_t>();
    if constexpr (Bytes & 2)
        chunk.template operator()<std::uint16_t>();
    if constexpr (Bytes & 1)
        chunk.template operator()<std::uint8_t>();
}


template<typename Derived, typename Real>
struct MultiPrecisionOps
{
    using real_type = Real;

    auto& derived() { return static_cast<Derived&>(*this); }

    Derived& operator+=(auto const& v)
    {
        return derived() = static_cast<Real>(derived()) + static_cast<Real>(v);
    }

    Derived& operator-=(auto const& v)
    {
        return derived() = static_cast<Real>(derived()) - static_cast<Real>(v);
    }

    Derived& operator*=(auto const& v)
    {
        return derived() = static_cast<Real>(derived()) * static_cast<Real>(v);
    }

    friend std::ostream& operator<<(std::ostream& out, Derived const& d)
    {
        return out << static_cast<Real>(d);
    }
};


/** The Bytes most significant bytes of an IEEE-754 double:
 *   sign, 11 exponent bits and 8 * Bytes - 12 mantissa bits.
 *   Dropped mantissa bits are rounded to nearest, ties to even.
 */
template<std::size_t Bytes>
struct TruncatedDouble : MultiPrecisionOps<TruncatedDouble<Bytes>, double>
{
    static_assert(Bytes >= 2 and Bytes <= 7);
    static constexpr std::size_t bytes          = Bytes;
    static constexpr std::uint32_t dropped_bits = 64 - 8 * Bytes;

    constexpr TruncatedDouble() = default;
    TruncatedDouble(double const v) { *this = v; }

    TruncatedDouble& operator=(double const v)
    {
        auto u                  = std::bit_cast<std::uint64_t>(v);
        std::uint64_t const lsb = (u >> dropped_bits) & 1;
        // a carry out of the mantissa increments the exponent, which is the correct rounding
        u += (std::uint64_t{1} << (dropped_bits - 1)) - 1 + lsb;
        u >>= dropped_bits;
        store_low_bytes<Bytes>(data.data(), u);
        return *this;
    }

    operator double() const
    {
        return std::bit_cast<double>(load_low_bytes<Bytes>(data.data()) << dropped_bits);
    }

    std::array<std::uint8_t, Bytes> data{};
};


/** Unsigned fixed point for values in [0, 1), resolution 2^-(8 * Bytes).
 *   Rounds to nearest and clamps to [0, 1 - 2^-(8 * Bytes)], so a value that would round
 *   up to 1 can never wrap to 0 without its cell being incremented.
 *   Bytes <= 6 keeps the decoded value exactly representable as a double.
 */
template<std::size_t Bytes>
struct FixedPointUnit : MultiPrecisionOps<FixedPointUnit<Bytes>, double>
{
    static_assert(Bytes >= 1 and Bytes <= 6);
    static constexpr std::size_t bytes  = Bytes;
    static constexpr std::uint32_t bits = 8 * Bytes;
    static constexpr std::uint64_t max  = (std::uint64_t{1} << bits) - 1;
    static constexpr double scale       = static_cast<double>(std::uint64_t{1} << bits);

    constexpr FixedPointUnit() = default;
    FixedPointUnit(double const v) { *this = v; }

    FixedPointUnit& operator=(double const v)
    {
        double const s        = std::nearbyint(v * scale);
        std::uint64_t const u = !(s > 0)                        ? 0 // also NaN
                                : s >= static_cast<double>(max) ? max
                                                                : static_cast<std::uint64_t>(s);
        store_low_bytes<Bytes>(data.data(), u);
        return *this;
    }

    operator double() const
    {
        return static_cast<double>(load_low_bytes<Bytes>(data.data())) / scale;
    }

    std::array<std::uint8_t, Bytes> data{};
};



// T::real_type for the storage types above, else T
template<typename T>
struct real_type_of
{
    using type = T;
};
template<typename T>
    requires requires { typename T::real_type; }
struct real_type_of<T>
{
    using type = typename T::real_type;
};
template<typename T>
using real_type_t = typename real_type_of<T>::type;


/** Fixed size array whose storage precision is selectable per instantiation,
 *   T may be a plain type or any of the above. Converts to and from any other precision,
 *   element-wise via the real type, so mixed precision interop is hidden from callers.
 *   Same layout as std::array<T, N>.
 */
template<typename T, std::size_t N>
struct MultiPrecisionArray
{
    using value_type     = T;
    using real_type      = real_type_t<T>;
    using iterator       = typename std::array<T, N>::iterator;
    using const_iterator = typename std::array<T, N>::const_iterator;

    MultiPrecisionArray() = default;

    explicit MultiPrecisionArray(std::array<real_type, N> const& from)
        : data{array_cast<T>(from)}
    {
    }

    MultiPrecisionArray& operator=(std::array<real_type, N> const& from)
    {
        data = array_cast<T>(from);
        return *this;
    }

    template<typename U>
    MultiPrecisionArray(MultiPrecisionArray<U, N> const& that)
        : data{array_cast<T>(that.data)}
    {
    }

    // exact, not float_equals: e.g. particle == means identity (round trips through packing,
    //  messaging, copies), not closeness. Storage types round trip exactly, but mixed precision
    //  only equals if exactly representable in both. Value compare, so +0 == -0, NaN != NaN
    template<typename U>
    NO_DISCARD bool operator==(MultiPrecisionArray<U, N> const& that) const
    {
        return array_equals(data, that.data);
    }

    NO_DISCARD auto& operator[](std::size_t const i) { return data[i]; }
    NO_DISCARD auto& operator[](std::size_t const i) const { return data[i]; }

    NO_DISCARD auto begin() { return data.begin(); }
    NO_DISCARD auto begin() const { return data.begin(); }
    NO_DISCARD auto end() { return data.end(); }
    NO_DISCARD auto end() const { return data.end(); }

    NO_DISCARD static constexpr std::size_t size() { return N; }

    std::array<T, N> data{};
};

template<std::size_t Bytes, std::size_t N>
using FixedPointArray = MultiPrecisionArray<FixedPointUnit<Bytes>, N>;

template<std::size_t Bytes, std::size_t N>
using TruncatedFPArray = MultiPrecisionArray<TruncatedDouble<Bytes>, N>;


} // namespace PHARE::core


#endif /* PHARE_CORE_UTILITIES_MULTI_PRECISION_HPP */
