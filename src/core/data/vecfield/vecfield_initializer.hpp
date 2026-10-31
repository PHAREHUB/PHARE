#ifndef VECFIELD_INITIALIZER_HPP
#define VECFIELD_INITIALIZER_HPP

#include "initializer/data_provider.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/data/field/initializers/field_user_initializer.hpp"

#include <string>
#include <cstddef>
#include <stdexcept>


namespace PHARE
{
namespace core
{
    template<std::size_t dimension>
    class VecFieldInitializer
    {
        static bool validate_dict(initializer::PHAREDict const& dict);

        VecFieldInitializer(initializer::PHAREDict const& dict, bool)
            : x_{dict["x_component"].template to<initializer::InitFunction<dimension>>()}
            , y_{dict["y_component"].template to<initializer::InitFunction<dimension>>()}
            , z_{dict["z_component"].template to<initializer::InitFunction<dimension>>()}
        {
        }

    public:
        VecFieldInitializer() = default;

        VecFieldInitializer(initializer::PHAREDict const& dict)
            : VecFieldInitializer{dict, validate_dict(dict)}
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

} // namespace core

} // namespace PHARE


namespace PHARE::core
{
template<std::size_t dim>
bool VecFieldInitializer<dim>::validate_dict(initializer::PHAREDict const& dict)
{
    for (auto const* key : {"x_component", "y_component", "z_component"})
        if (!dict.contains(key))
            throw std::runtime_error(std::string{"VecFieldInitializer: missing key "} + key);
    return true;
}

} // namespace PHARE::core

#endif // VECFIELD_INITIALIZER_HPP
