#ifndef PHARE_CORE_INNER_BOUNDARY_GHOST_ELEM_PACK_HPP
#define PHARE_CORE_INNER_BOUNDARY_GHOST_ELEM_PACK_HPP

#include "core/def.hpp"
#include "core/inner_boundary/inner_boundary_defs.hpp"

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace PHARE::core
{

/**
 * @brief Per-patch view over the precomputed ghost-element vectors.
 *
 * Mirrors the role of ParticlesPack for ParticlesData: a thin view object held
 * by the ResourcesUser (here InnerBoundaryMeshData) whose _data pointer is
 * rebound by ResourcesManager::setOnPatch to the array owned by the patch's
 * GhostElemPatchData.
 */
template<std::size_t dim>
struct GhostElemPack
{
    static constexpr std::size_t dimension      = dim;
    static constexpr std::size_t num_elem_types = (std::size_t{1} << dim);

    using ghost_elem_data_type  = GhostElemData<dim>;
    using ghost_elem_array_type = std::array<std::vector<ghost_elem_data_type>, num_elem_types>;

    std::string             _name;
    ghost_elem_array_type*  _data{nullptr};

    GhostElemPack() = default;
    explicit GhostElemPack(std::string name)
        : _name{std::move(name)}
    {
    }

    void setBuffer(GhostElemPack* source)
    {
        (*this) = source ? *source : GhostElemPack{_name};
    }

    auto& name() const { return _name; }

    NO_DISCARD bool isUsable() const { return _data != nullptr; }
    NO_DISCARD bool isSettable() const { return _data == nullptr; }

    NO_DISCARD std::vector<ghost_elem_data_type>& operator[](std::size_t i)
    {
        return data_or_throw_()[i];
    }
    NO_DISCARD std::vector<ghost_elem_data_type> const& operator[](std::size_t i) const
    {
        return data_or_throw_()[i];
    }

    NO_DISCARD auto begin() { return data_or_throw_().begin(); }
    NO_DISCARD auto end()   { return data_or_throw_().end(); }
    NO_DISCARD auto begin() const { return data_or_throw_().begin(); }
    NO_DISCARD auto end()   const { return data_or_throw_().end(); }

private:
    ghost_elem_array_type& data_or_throw_()
    {
        if (!_data)
            throw std::runtime_error("GhostElemPack '" + _name + "' is not bound to a patch");
        return *_data;
    }
    ghost_elem_array_type const& data_or_throw_() const
    {
        if (!_data)
            throw std::runtime_error("GhostElemPack '" + _name + "' is not bound to a patch");
        return *_data;
    }
};

} // namespace PHARE::core

#endif // PHARE_CORE_INNER_BOUNDARY_GHOST_ELEM_PACK_HPP
