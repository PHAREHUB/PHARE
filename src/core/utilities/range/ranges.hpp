#ifndef PHARE_CORE_UTILITITES_RANGES_HPP
#define PHARE_CORE_UTILITITES_RANGES_HPP

#include <vector>
#include <iterator>
#include <cstddef>
#include <iterator>
#include <algorithm>

#include "core/utilities/types.hpp"
#include "core/utilities/range/range.hpp"

namespace PHARE::core
{
template<typename Container, typename Iterator = typename Container::iterator>
auto ranges(Container& container, std::size_t batch_size, std::size_t offset = 0)
{
    // using Iterator = std::conditional_t<std::is_const_v<decltype(container)>, typename
    // Container::const_iterator, typename Container::iterator>;

    std::size_t size     = container.size() - offset;
    std::size_t modulo   = size % batch_size;
    std::size_t n_ranges = size / batch_size;
    std::vector<Range<Iterator>> ranges;
    ;
    ranges.reserve(n_ranges + (modulo > 0 ? 1 : 0));
    auto begin = container.begin();
    for (std::size_t i = 0; i < n_ranges; ++i)
    {
        begin = container.begin() + offset + i * batch_size;
        ranges.emplace_back(makeRange(begin, begin + batch_size));
    }
    if (modulo > 0)
    {
        begin = container.begin() + offset + n_ranges * batch_size;
        ranges.emplace_back(makeRange(begin, begin + modulo));
    }
    return ranges;
}
template<typename Container, typename Iterator = typename Container::iterator>
auto reverse_ranges(Container& container, std::size_t batch_size, std::size_t offset = 0)
{
    auto vec = ranges(container, batch_size, offset);
    std::reverse(vec.begin(), vec.end());
    return vec;
}

namespace detail
{
    struct BalancedRangesDefaults
    {
        auto static constexpr accessor = [](auto& el) -> auto& { return el; };
        using Accessor                 = decltype(accessor);
        auto static constexpr builder  = [](auto range, auto& el) { return range; };
        using Builder                  = decltype(builder);
    };
} // namespace detail


template<typename Container, typename Defaults = detail::BalancedRangesDefaults,
         typename Accessor = typename Defaults::Accessor,
         typename Builder  = typename Defaults::Builder>
auto make_balanced_ranges(Container& container, std::size_t split = 10,
                          Accessor accessor = Defaults::accessor,
                          Builder builder   = Defaults::builder)
{
    using Array = std::decay_t<std::result_of_t<Accessor&(typename Container::value_type&)>>;

    using Range_t = Range<typename Array::iterator>;
    using BuildRange_t
        = std::decay_t<std::result_of_t<Builder&(Range_t const, typename Container::value_type&)>>;

    auto arrays = core::generate([&](auto& el) { return &accessor(el); }, container);

    std::vector<std::vector<BuildRange_t>> ranges_vec(split);
    std::size_t arr_idx = 0, offset = 0;

    auto load_ranges = [&](auto const t_idx, auto const size) {
        auto rem = size;
        while (rem > 0)
        {
            auto& arr  = *arrays[arr_idx];
            auto range = makeRange(arr.begin() + offset, arr.end());
            auto op    = rem >= range.size() ? range.size() : rem;

            ranges_vec[t_idx].emplace_back(
                builder(makeRange(range.begin(), range.begin() + op), container[arr_idx]));

            offset += op;
            if (rem >= range.size())
            {
                ++arr_idx;
                offset = 0;
            }
            rem -= op;
        }
    };

    auto const total = sum_from(container, [&](auto const& el) { return accessor(el).size(); });
    std::size_t const each   = total / split;
    std::size_t const modulo = total % split;

    // give any remaining to index 0, assuming it's the main thread
    load_ranges(0, each + modulo);
    for (std::size_t i = 1; i < ranges_vec.size(); ++i)
        load_ranges(i, each);

    return ranges_vec;
}

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITITES_RANGES_HPP
