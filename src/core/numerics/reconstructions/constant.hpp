#ifndef CORE_NUMERICS_RECONSTRUCTION_CONSTANT_HPP
#define CORE_NUMERICS_RECONSTRUCTION_CONSTANT_HPP

#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"
#include "core/numerics/reconstructions/stencil.hpp"
#include <utility>

namespace PHARE::core
{
template<typename GridLayout, typename SlopeLimiter = void>
class ConstantReconstruction
{
public:
    static constexpr auto nghosts = 1;

    using GridLayout_t = GridLayout;

    template<auto direction, typename Field>
    static auto reconstruct(Field const& F, MeshIndex<Field::dimension> index)
    {
        return std::make_pair(F(GridLayout::template previous<direction>(index)), F(index));
    }

    template<auto direction, auto Projection, typename Field>
    static auto center_reconstruct(Field const& U, MeshIndex<Field::dimension> index)
    {
        auto u_1 = GridLayout::template project<Projection>(
            U, GridLayout::template previous<direction>(index));
        auto u = GridLayout::template project<Projection>(U, index);

        return std::make_pair(u_1, u);
    }
    // values at the lower and upper faces of cell `index`, from the one stencil centred on it
    template<auto direction, typename Field>
    static auto cell_faces(Field const& F, MeshIndex<Field::dimension> index)
    {
        return faces_of(centred_stencil<GridLayout, nghosts - 1, direction>(F, index));
    }

    // lower and upper face values of the centre cell of a 2 * (nghosts - 1) + 1 stencil
    static auto faces_of(auto const& u) { return std::make_pair(u[0], u[0]); }
};

} // namespace PHARE::core
#endif
