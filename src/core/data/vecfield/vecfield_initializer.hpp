#ifndef VECFIELD_INITIALIZER_HPP
#define VECFIELD_INITIALIZER_HPP

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "initializer/data_provider.hpp"

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
          
            v.getComponent(Component::X).initialize(layout, x_);
            v.getComponent(Component::Y).initialize(layout, y_);
            v.getComponent(Component::Z).initialize(layout, z_);
        }

    private:
        initializer::InitFunction<dimension> x_;
        initializer::InitFunction<dimension> y_;
        initializer::InitFunction<dimension> z_;
    };

} // namespace core

} // namespace PHARE

#endif // VECFIELD_INITIALIZER_HPP
