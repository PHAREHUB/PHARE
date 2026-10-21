#ifndef PHARE_CORE_DATA_VECFIELD_VECTOR_POTENTIAL_INITIALIZER_HPP
#define PHARE_CORE_DATA_VECFIELD_VECTOR_POTENTIAL_INITIALIZER_HPP

#include "core/data/grid/grid.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/field/initializers/field_user_initializer.hpp"
#include "core/numerics/curl/edge_to_face_curl.hpp"
#include "initializer/data_provider.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>

namespace PHARE::core
{
template<std::size_t dimension>
class VectorPotentialInitializer
{
    using InitFn = initializer::InitFunction<dimension>;

public:
    VectorPotentialInitializer(initializer::PHAREDict const& dict)
        : potential_{dict["vector_potential"]["x_component"].template to<InitFn>(),
                     dict["vector_potential"]["y_component"].template to<InitFn>(),
                     dict["vector_potential"]["z_component"].template to<InitFn>()}
        , direct_{read_optional_(dict, "x_component"), read_optional_(dict, "y_component"),
                  read_optional_(dict, "z_component")}
    {
        if constexpr (dimension == 1)
            throw std::runtime_error("VectorPotentialInitializer: needs dimension > 1");
    }

    template<typename VecField, typename GridLayout>
    void initialize(VecField& B, GridLayout const& layout) const
    {
        static_assert(GridLayout::dimension == VecField::dimension,
                      "dimension mismatch between vecfield and gridlayout");

        if constexpr (dimension > 1)
        {
            using Scalar = VecField::field_type::physical_quantity_type;
            using Grid_t = Grid<NdArrayVector<dimension, double>, Scalar>;

            auto constexpr qts = edge_quantities<Scalar>();
            std::array<Grid_t, 3> A{Grid_t{"Ax", layout, qts[0]}, Grid_t{"Ay", layout, qts[1]},
                                    Grid_t{"Az", layout, qts[2]}};

            for (std::size_t c = 0; c < 3; ++c)
                FieldUserFunctionInitializer::initialize(*&A[c], layout, potential_[c]);

            curl_edges_to_faces(layout, std::array{*&A[0], *&A[1], *&A[2]}, B);

            for (std::size_t c = 0; c < 3; ++c)
                if (direct_[c])
                    FieldUserFunctionInitializer::initialize(B[c], layout, *direct_[c]);
        }
    }

private:
    static std::optional<InitFn> read_optional_(initializer::PHAREDict const& dict,
                                                std::string const& key)
    {
        if (!dict.contains(key))
            return std::nullopt;
        return dict[key].template to<InitFn>();
    }

    std::array<InitFn, 3> potential_;
    std::array<std::optional<InitFn>, 3> direct_;
};

} // namespace PHARE::core

#endif // PHARE_CORE_DATA_VECFIELD_VECTOR_POTENTIAL_INITIALIZER_HPP
