#ifndef VECFIELD_INITIALIZER_H
#define VECFIELD_INITIALIZER_H

#include "data/grid/gridlayoutdefs.h"
#include "data/vecfield/vecfield_component.h"
#include "data_provider.h"

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

        VecFieldInitializer(initializer::PHAREDict<dimension> dict)
            : x_{dict["x_component"].template to<initializer::ScalarFunction<dimension>>()}
            , y_{dict["y_component"].template to<initializer::ScalarFunction<dimension>>()}
            , z_{dict["z_component"].template to<initializer::ScalarFunction<dimension>>()}
        {
        }


        template<typename VecField, typename GridLayout>
        void initialize(VecField& v, GridLayout const& layout)
        {
            static_assert(GridLayout::dimension == VecField::dimension,
                          "dimension mismatch between vecfield and gridlayout");

            auto& vx = v.getComponent(Component::X);
            auto& vy = v.getComponent(Component::Y);
            auto& vz = v.getComponent(Component::Z);

            initializeComponent_(vx, layout, x_);
            initializeComponent_(vy, layout, y_);
            initializeComponent_(vz, layout, z_);
        }

    private:
        template<typename Field, typename GridLayout>
        void initializeComponent_(Field& field, GridLayout const& layout,
                                  initializer::ScalarFunction<dimension> const& init)
        {
            auto psi_X = layout.ghostStartIndex(field, Direction::X);
            auto pei_X = layout.ghostEndIndex(field, Direction::X);

            for (auto ix = psi_X; ix <= pei_X; ++ix)
            {
                if constexpr (dimension > 1)
                {
                    auto psi_Y = layout.ghostStartIndex(field, Direction::Y);
                    auto pei_Y = layout.ghostEndIndex(field, Direction::Y);

                    for (auto iy = psi_Y; iy <= pei_Y; ++iy)
                    {
                        if constexpr (dimension > 2)
                        {
                            auto psi_Z = layout.ghostStartIndex(field, Direction::Z);
                            auto pei_Z = layout.ghostEndIndex(field, Direction::Z);

                            for (auto iz = psi_Z; iz <= pei_Z; ++iz)
                            {
                                auto pos = layout.fieldNodeCoordinates(field, layout.origin(), ix,
                                                                       iy, iz);
                                field(ix, iy, iz) = init(pos[0], pos[1], pos[2]);
                            }
                        }
                        else // 2D
                        {
                            auto pos = layout.fieldNodeCoordinates(field, layout.origin(), ix, iy);
                            field(ix, iy) = init(pos[0], pos[1]);
                        }
                    }
                }
                else // 1D
                {
                    auto pos  = layout.fieldNodeCoordinates(field, layout.origin(), ix);
                    field(ix) = init(pos[0]);
                }
            }
        }




        initializer::ScalarFunction<dimension> x_;
        initializer::ScalarFunction<dimension> y_;
        initializer::ScalarFunction<dimension> z_;
    };

} // namespace core

} // namespace PHARE

#endif // VECFIELD_INITIALIZER_H
