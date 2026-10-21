#ifndef PHARE_CORE_DATA_FIELD_RECENTERED_FIELD_HPP
#define PHARE_CORE_DATA_FIELD_RECENTERED_FIELD_HPP

#include "core/data/grid/gridlayoutdefs.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

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
    {
        static_assert(sizeof...(Indexes) == dimension);
        return (*this)(
            std::array<std::uint32_t, dimension>{static_cast<std::uint32_t>(indexes)...});
    }

    value_type operator()(std::array<std::uint32_t, dimension> const& index) const
    {
        value_type value = 0;
        for (std::uint32_t corner = 0; corner < (1u << dimension); ++corner)
        {
            auto node         = index;
            value_type weight = 1;
            bool inStencil    = true;
            for (std::size_t d = 0; d < dimension; ++d)
            {
                std::uint32_t const upper = (corner >> d) & 1u;
                if (source_[d] == destination_[d])
                {
                    if (upper)
                    {
                        inStencil = false;
                        break;
                    }
                }
                else
                {
                    weight *= 0.5;
                    node[d] = source_[d] == QtyCentering::primal ? node[d] + upper
                                                                 : node[d] + upper - 1;
                }
            }
            if (inStencil)
                value += weight * field_(node);
        }
        return value;
    }

    auto const& centering() const { return destination_; }

private:
    Field field_;
    std::array<QtyCentering, dimension> source_;
    std::array<QtyCentering, dimension> destination_;
};

} // namespace PHARE::core

#endif
