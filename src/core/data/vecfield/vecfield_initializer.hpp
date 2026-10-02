#ifndef VECFIELD_INITIALIZER_HPP
#define VECFIELD_INITIALIZER_HPP

#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/numerics/curl/edge_to_face_curl.hpp"
#include "initializer/data_provider.hpp"
#include "core/data/field/initializers/field_user_initializer.hpp"

#include <array>
#include <optional>
#include <string>

namespace PHARE
{
namespace core
{
    /** @brief Initializes a VecField from user functions.
     *
     * B mode (default): each of x/y/z_component is sampled at the component's centring.
     * A mode (dict holds a "vector_potential" subtree): the vector potential components are
     * sampled at Yee edges and the VecField is set to their discrete curl. Any x/y/z_component
     * also present then overrides the curl result for that component (2D out-of-plane B).
     */
    template<std::size_t dimension>
    class VecFieldInitializer
    {
        using InitFn = initializer::InitFunction<dimension>;

    public:
        VecFieldInitializer() = default;

        VecFieldInitializer(initializer::PHAREDict const& dict)
            : direct_{read_optional_(dict, "x_component"), read_optional_(dict, "y_component"),
                      read_optional_(dict, "z_component")}
        {
            if (dict.contains("vector_potential"))
            {
                auto const& vp = dict["vector_potential"];
                vecpot_        = std::array<InitFn, 3>{vp["x_component"].template to<InitFn>(),
                                                       vp["y_component"].template to<InitFn>(),
                                                       vp["z_component"].template to<InitFn>()};
            }
            else
            {
                for (auto const& component : direct_)
                    if (!component)
                        throw std::runtime_error(
                            "VecFieldInitializer: x/y/z_component all required without "
                            "vector_potential");
            }
        }


        template<typename VecField, typename GridLayout>
        void initialize(VecField& v, GridLayout const& layout)
        {
            static_assert(GridLayout::dimension == VecField::dimension,
                          "dimension mismatch between vecfield and gridlayout");

            if (!vecpot_)
            {
                FieldUserFunctionInitializer::initialize(v.getComponent(Component::X), layout,
                                                         direct_[0].value());
                FieldUserFunctionInitializer::initialize(v.getComponent(Component::Y), layout,
                                                         direct_[1].value());
                FieldUserFunctionInitializer::initialize(v.getComponent(Component::Z), layout,
                                                         direct_[2].value());
                return;
            }

            if constexpr (dimension == 1)
                throw std::runtime_error("VecFieldInitializer: vector potential needs dim > 1");
            else
                initialize_from_vector_potential_(v, layout);
        }

    private:
        static std::optional<InitFn> read_optional_(initializer::PHAREDict const& dict,
                                                    std::string const& key)
        {
            if (dict.contains(key))
                return dict[key].template to<InitFn>();
            return std::nullopt;
        }

        template<typename VecField, typename GridLayout>
        void initialize_from_vector_potential_(VecField& v, GridLayout const& layout) const
        {
            using Scalar = typename VecField::field_type::physical_quantity_type;
            using Grid_t = Grid<NdArrayVector<dimension, double>, Scalar>;

            auto constexpr qts = edge_quantities<Scalar>();
            Grid_t Ax{"Ax", layout, qts[0]};
            Grid_t Ay{"Ay", layout, qts[1]};
            Grid_t Az{"Az", layout, qts[2]};

            FieldUserFunctionInitializer::initialize(*&Ax, layout, (*vecpot_)[0]);
            FieldUserFunctionInitializer::initialize(*&Ay, layout, (*vecpot_)[1]);
            FieldUserFunctionInitializer::initialize(*&Az, layout, (*vecpot_)[2]);

            curl_edges_to_faces(layout, std::array{*&Ax, *&Ay, *&Az}, v);

            auto constexpr components = std::array{Component::X, Component::Y, Component::Z};
            for (std::size_t c = 0; c < 3; ++c)
                if (direct_[c])
                    FieldUserFunctionInitializer::initialize(v.getComponent(components[c]),
                                                             layout, *direct_[c]);
        }

        std::array<std::optional<InitFn>, 3> direct_;
        std::optional<std::array<InitFn, 3>> vecpot_;
    };

} // namespace core

} // namespace PHARE

#endif // VECFIELD_INITIALIZER_HPP
