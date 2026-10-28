#ifndef PHARE_CORE_UTILITIES_ENUM_TABLE_HPP
#define PHARE_CORE_UTILITIES_ENUM_TABLE_HPP

#include <array>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <utility>

namespace PHARE::core
{

//! builds a std::array indexed by Enum from {key, value} entries given in any order
//! every key in [0, Enum::count) must appear exactly once, so the table cannot silently
//! go out of order with the enum: in a constant expression a missing, duplicate or out of
//! range key is a compile error
template<typename Enum, typename Value>
constexpr auto make_enum_table(std::initializer_list<std::pair<Enum, Value>> const entries)
{
    constexpr auto size = static_cast<std::size_t>(Enum::count);

    std::array<Value, size> table{};
    std::array<bool, size> isSet{};

    for (auto const& [key, value] : entries)
    {
        auto const idx = static_cast<std::size_t>(key);
        if (idx >= size)
            throw std::runtime_error("make_enum_table: key out of range");
        if (isSet[idx])
            throw std::runtime_error("make_enum_table: duplicate key");
        table[idx] = value;
        isSet[idx] = true;
    }

    for (auto const set : isSet)
        if (!set)
            throw std::runtime_error("make_enum_table: missing key");

    return table;
}

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_ENUM_TABLE_HPP
