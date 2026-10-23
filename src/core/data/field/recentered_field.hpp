#ifndef PHARE_CORE_DATA_FIELD_RECENTERED_FIELD_HPP
#define PHARE_CORE_DATA_FIELD_RECENTERED_FIELD_HPP

#include "core/data/grid/gridlayoutdefs.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <type_traits>

namespace PHARE::core
{

/**
 * @brief Read-only view of a field as if it had another centering.
 *
 * Values are computed when read by linear interpolation. Nothing is allocated.
 *
 * Primal -> dual is valid on every index. Dual -> primal is not valid on the outermost
 * primal ghost node of each side.
 *
 * Holds a copy of the Field (pointer and shape): rebuild it once the field is set on
 * another patch.
 */
template<typename Field>
class RecenteredField
{
public:
    static constexpr auto dimension = Field::dimension;
    using value_type                = Field::value_type;

    template<typename GridLayout>
    RecenteredField(Field const& field, std::array<QtyCentering, dimension> const& destination,
                    GridLayout const& layout)
        : field_{field}
        , source_{layout.centering(field)}
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
    /// recursive helper to compute the tensor product of the 1D recentering operation
    template<std::size_t d>
    value_type averageFrom_(std::array<std::uint32_t, dimension> const& node) const
    {
        if constexpr (d == dimension)
            // end of the recursion
            return field_(node);
        else
        {
            if (source_[d] == destination_[d])
                // no recentering needed, return the result of the next recursion
                return averageFrom_<d + 1>(node);

            // the centering differs in the current direction `d`; so return the average
            // of the results of the next recursion, applied to lower and upper nodes.
            auto lower = node;
            auto upper = node;
            if (source_[d] == QtyCentering::primal)
                upper[d] += 1;
            else
                lower[d] -= 1;
            return 0.5 * (averageFrom_<d + 1>(lower) + averageFrom_<d + 1>(upper));
        }
    }

    Field field_;                                     ///< source field
    std::array<QtyCentering, dimension> source_;      ///< source field's centering
    std::array<QtyCentering, dimension> destination_; ///< destination centering
};

} // namespace PHARE::core

#endif
