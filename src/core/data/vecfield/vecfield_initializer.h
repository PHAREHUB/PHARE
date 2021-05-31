#ifndef VECFIELD_INITIALIZER_H
#define VECFIELD_INITIALIZER_H

#include "core/data/grid/gridlayoutdefs.h"
#include "core/data/vecfield/vecfield_component.h"
#include "initializer/data_provider.h"

#include <array>

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

            initializeComponent_(v.getComponent(Component::X), layout, x_);
            initializeComponent_(v.getComponent(Component::Y), layout, y_);
            initializeComponent_(v.getComponent(Component::Z), layout, z_);
        }

    private:
        template<typename Field, typename GridLayout>
        void initializeComponent_(Field& field, GridLayout const& layout,
                                  initializer::InitFunction<dimension> const& init)
        {
            auto const indices = layout.ghostStartToEndIndices(field, /*includeEnd=*/true);
            auto const coords  = layout.template indexesToCoordVectors</*WithField=*/true>(
                indices, field, [](auto& gridLayout, auto& field_, auto const&... args) {
                    return gridLayout.fieldNodeCoordinates(field_, gridLayout.origin(), args...);
                });

            std::shared_ptr<Span<double>> gridPtr // keep grid data alive
                = std::apply([&](auto&... args) { return init(args...); }, coords);
            Span<double>& grid = *gridPtr;

            for (std::size_t cell_idx = 0; cell_idx < indices.size(); cell_idx++)
                std::apply([&](auto&... args) { field(args...) = grid[cell_idx]; },
                           indices[cell_idx]);
        }



        initializer::InitFunction<dimension> x_;
        initializer::InitFunction<dimension> y_;
        initializer::InitFunction<dimension> z_;
    };

} // namespace core

} // namespace PHARE

#endif // VECFIELD_INITIALIZER_H
