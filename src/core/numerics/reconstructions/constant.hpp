#ifndef CORE_NUMERICS_RECONSTRUCTION_CONSTANT_HPP
#define CORE_NUMERICS_RECONSTRUCTION_CONSTANT_HPP

#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"
#include <utility>

namespace PHARE::core
{
template<typename GridLayout, typename SlopeLimiter = void>
class ConstantReconstruction
{
public:
    using GridLayout_t = GridLayout;

    template<auto direction, typename Field>
    static auto reconstruct(Field const& F, MeshIndex<Field::dimension> index)
    {
        return std::make_pair(F(GridLayout::template previous<direction>(index)), F(index));
    }

    template<auto direction, typename Field>
    static auto center_reconstruct(Field const& U, MeshIndex<Field::dimension> index,
                                   auto projection)
    {
        auto u_1
            = GridLayout::project(U, GridLayout::template previous<direction>(index), projection);
        auto u = GridLayout::project(U, index, projection);

        return std::make_pair(u_1, u);
    }
};

} // namespace PHARE::core
#endif
