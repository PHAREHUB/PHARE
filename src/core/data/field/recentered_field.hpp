#ifndef PHARE_CORE_DATA_FIELD_RECENTERED_FIELD_HPP
#define PHARE_CORE_DATA_FIELD_RECENTERED_FIELD_HPP

#include "core/data/grid/gridlayoutdefs.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace PHARE::core
{
template<typename GridLayout, typename Field>
class RecenteredField
{
public:
    static constexpr auto dimension = GridLayout::dimension;
    using value_type                = Field::value_type;

    RecenteredField(Field const& field, std::array<QtyCentering, dimension> const& destination)
        : field_{field}
        , source_{GridLayout::centering(field)}
        , destination_{destination}
    {
    }

    template<typename... Indexes>
    value_type operator()(Indexes... indexes) const
        requires(sizeof...(Indexes) == dimension && (std::is_integral_v<Indexes> && ...))
    {
        return (*this)(
            std::array<std::uint32_t, dimension>{static_cast<std::uint32_t>(indexes)...});
    }

    value_type operator()(std::array<std::uint32_t, dimension> const& index) const
    {
        return averageFrom_<0>(index);
    }

    auto const& centering() const { return destination_; }

private:
    template<std::size_t d>
    value_type averageFrom_(std::array<std::uint32_t, dimension> const& node) const
    {
        if constexpr (d == dimension)
            return field_(node);
        else
        {
            if (source_[d] == destination_[d])
                return averageFrom_<d + 1>(node);

            auto lower = node;
            auto upper = node;
            if (source_[d] == QtyCentering::primal)
                upper[d] += 1;
            else
                lower[d] -= 1;
            return (averageFrom_<d + 1>(lower) + averageFrom_<d + 1>(upper)) / 2;
        }
    }

    Field field_;
    std::array<QtyCentering, dimension> source_;
    std::array<QtyCentering, dimension> destination_;
};

} // namespace PHARE::core

#endif
