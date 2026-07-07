#ifndef VECFIELD_INITIALIZER_HPP
#define VECFIELD_INITIALIZER_HPP

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "initializer/data_provider.hpp"
#include "core/data/field/initializers/field_user_initializer.hpp"

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

            FieldUserFunctionInitializer::initialize(v.getComponent(Component::X), layout, x_);
            FieldUserFunctionInitializer::initialize(v.getComponent(Component::Y), layout, y_);
            FieldUserFunctionInitializer::initialize(v.getComponent(Component::Z), layout, z_);
        }

    private:
        initializer::InitFunction<dimension> x_;
        initializer::InitFunction<dimension> y_;
        initializer::InitFunction<dimension> z_;
    };


    // Like VecFieldInitializer but for space+time component functions f(x,t). Used to (re-)stamp a
    // time-dependent background field B0(x,t) (or its time derivative dB0/dt) at a given time.
    template<std::size_t dimension>
    class SpaceTimeVecFieldInitializer
    {
    public:
        SpaceTimeVecFieldInitializer() = default;

        SpaceTimeVecFieldInitializer(initializer::PHAREDict const& dict)
            : x_{dict["x_component"].template to<initializer::SpaceTimeFunction<dimension>>()}
            , y_{dict["y_component"].template to<initializer::SpaceTimeFunction<dimension>>()}
            , z_{dict["z_component"].template to<initializer::SpaceTimeFunction<dimension>>()}
        {
        }

        template<typename VecField, typename GridLayout>
        void initialize(VecField& v, GridLayout const& layout, double time)
        {
            static_assert(GridLayout::dimension == VecField::dimension,
                          "dimension mismatch between vecfield and gridlayout");

            FieldUserFunctionInitializer::initialize(v.getComponent(Component::X), layout, x_, time);
            FieldUserFunctionInitializer::initialize(v.getComponent(Component::Y), layout, y_, time);
            FieldUserFunctionInitializer::initialize(v.getComponent(Component::Z), layout, z_, time);
        }

    private:
        initializer::SpaceTimeFunction<dimension> x_;
        initializer::SpaceTimeFunction<dimension> y_;
        initializer::SpaceTimeFunction<dimension> z_;
    };

} // namespace core

} // namespace PHARE

#endif // VECFIELD_INITIALIZER_HPP
