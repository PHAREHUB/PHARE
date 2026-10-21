#ifndef PHARE_CORE_DATA_VECFIELD_MAGNETIC_FIELD_INITIALIZER_HPP
#define PHARE_CORE_DATA_VECFIELD_MAGNETIC_FIELD_INITIALIZER_HPP

#include "core/data/vecfield/vecfield_initializer.hpp"
#include "core/data/vecfield/vector_potential_initializer.hpp"
#include "initializer/data_provider.hpp"

#include <cstddef>
#include <variant>

namespace PHARE::core
{
template<std::size_t dimension>
class MagneticFieldInitializer
{
    using Components      = VecFieldInitializer<dimension>;
    using VectorPotential = VectorPotentialInitializer<dimension>;

public:
    MagneticFieldInitializer() = default;

    MagneticFieldInitializer(initializer::PHAREDict const& dict)
        : impl_{create(dict)}
    {
    }

    static std::variant<Components, VectorPotential> create(initializer::PHAREDict const& dict)
    {
        if (dict.contains("vector_potential"))
            return VectorPotential{dict};
        return Components{dict};
    }

    template<typename VecField, typename GridLayout>
    void initialize(VecField& B, GridLayout const& layout)
    {
        std::visit([&](auto& init) { init.initialize(B, layout); }, impl_);
    }

private:
    std::variant<Components, VectorPotential> impl_;
};

} // namespace PHARE::core

#endif // PHARE_CORE_DATA_VECFIELD_MAGNETIC_FIELD_INITIALIZER_HPP
