#ifndef VECFIELD_INITIALIZER_HPP
#define VECFIELD_INITIALIZER_HPP

#include "core/def/types.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "initializer/data_provider.hpp"

#include <tuple>
#include <cstdint>


namespace PHARE::core
{
template<typename Field, typename GridLayout, typename InitFunction>
void initialize_field(Field& field, GridLayout const& layout, InitFunction const& init)
{
    auto const indices = layout.ghostStartToEndIndices(field, /*includeEnd=*/true);
    auto const coords  = layout.template indexesToCoordVectors</*WithField=*/true>(
        indices, field, [](auto& gridLayout, auto& field_, auto const&... args) {
            return gridLayout.fieldNodeCoordinates(field_, gridLayout.origin(), args...);
        });

    // keep grid data alive
    auto grid = std::apply([&](auto const&... args) { return init(args...); }, coords);
    assert(field.size() == grid->size());

    if constexpr (Field::is_host_mem) // could this be a memcpy?
        for (std::size_t cell_idx = 0; cell_idx < indices.size(); cell_idx++)
            std::apply([&](auto&... args) { field(args...) = (*grid)[cell_idx]; },
                       indices[cell_idx]);
    else
    {
        if constexpr (CompileOptions::WithUmpire and CompileOptions::WithRAJA)
        {
            PHARE_WITH_RAJA(PHARE::core::raja::copy(field.data(), grid->data(), field.size()));
        }
    }
}
} // namespace PHARE::core

namespace PHARE
{
namespace core
{
    template<std::size_t dimension>
    class VecFieldInitializer
    {
    public:
        VecFieldInitializer() = default;

        VecFieldInitializer(initializer::PHAREDict const& dict)
            : x_{dict["x_component"].template to<initializer::InitFunction<dimension>>()}
            , y_{dict["y_component"].template to<initializer::InitFunction<dimension>>()}
            , z_{dict["z_component"].template to<initializer::InitFunction<dimension>>()}
        {
        }


        template<typename VecField, typename GridLayout>
        void initialize(VecField& v, GridLayout const& layout)
        {
            static_assert(GridLayout::dimension == VecField::dimension,
                          "dimension mismatch between vecfield and gridlayout");

            initialize_field(v.getComponent(Component::X), layout, x_);
            initialize_field(v.getComponent(Component::Y), layout, y_);
            initialize_field(v.getComponent(Component::Z), layout, z_);
        }

    private:
        initializer::InitFunction<dimension> x_;
        initializer::InitFunction<dimension> y_;
        initializer::InitFunction<dimension> z_;
    };

} // namespace core

} // namespace PHARE

#endif // VECFIELD_INITIALIZER_HPP
