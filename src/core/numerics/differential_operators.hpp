#ifndef PHARE_CORE_NUMERICS_DIFFERENTIAL_OPERATORS_HPP
#define PHARE_CORE_NUMERICS_DIFFERENTIAL_OPERATORS_HPP

#include "core/def.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/constants.hpp"
#include "core/utilities/index/index.hpp"
#include "core/utilities/types.hpp"

#include <array>
#include <cstddef>
#include <utility>
#include <stdexcept>
#include <type_traits>

namespace PHARE::core::detail
{
template<std::size_t d, typename Array, std::size_t dim>
double shifted(Array const& array, MeshIndex<dim> const& idx, int const shift)
{
    return [&]<std::size_t... i>(std::index_sequence<i...>) {
        return array((i == d ? idx[i] + shift : idx[i])...);
    }(std::make_index_sequence<dim>{});
}

template<typename GridLayout, typename Dst>
auto destination_centering(GridLayout const& layout, Dst const& dst)
{
    if constexpr (std::is_same_v<Dst, std::array<QtyCentering, GridLayout::dimension>>)
        return dst;
    else
        return layout.centering(dst);
}

template<Component c>
constexpr std::size_t component_index()
{
    static_assert(c == Component::X or c == Component::Y or c == Component::Z);
    return static_cast<std::size_t>(c);
}

template<auto dir, std::size_t dim>
constexpr bool direction_exists = static_cast<std::size_t>(dir) < dim;

struct Stencil
{
    int next;
    int prev;
    double coefficient;
};

constexpr Stencil stencil(QtyCentering const src, QtyCentering const dst)
{
    if (src == dst)
        return {1, -1, 0.5};
    if (src == QtyCentering::primal)
        return {1, 0, 1.};
    return {0, -1, 1.};
}
} // namespace PHARE::core::detail



namespace PHARE::core
{
template<auto dir, std::size_t dim>
class Partial;

template<auto dir, std::size_t dim>
    requires(detail::direction_exists<dir, dim>)
class Partial<dir, dim>
{
    static constexpr auto d = static_cast<std::size_t>(dir);

public:
    template<typename GridLayout, typename Src>
    Partial(GridLayout const& layout, Src const& src, std::array<QtyCentering, dim> const& dst)
        : view_{src.data(), src.shape()}
    {
        if (src.data() == nullptr)
            throw std::runtime_error("partial: the source has no data");

        auto const srcCentering = layout.centering(src);
        for (std::size_t t = 0; t < dim; ++t)
            if (t != d and srcCentering[t] != dst[t])
                throw std::runtime_error(
                    "partial: source and destination centerings differ in a transverse direction");

        auto const stencil = detail::stencil(srcCentering[d], dst[d]);
        next_              = stencil.next;
        prev_              = stencil.prev;
        scale_             = layout.inverseMeshSize(dir) * stencil.coefficient;
    }

    NO_DISCARD double operator()(MeshIndex<dim> const idx) const
    {
        return scale_
               * (detail::shifted<d>(view_, idx, next_) - detail::shifted<d>(view_, idx, prev_));
    }

private:
    NdArrayView<dim, double const> view_;
    int next_;
    int prev_;
    double scale_;
};

template<auto dir, std::size_t dim>
    requires(not detail::direction_exists<dir, dim>)
class Partial<dir, dim>
{
public:
    template<typename GridLayout, typename Src>
    Partial(GridLayout const&, Src const&, std::array<QtyCentering, dim> const&)
    {
    }

    NO_DISCARD double operator()(MeshIndex<dim> const) const { return 0.; }
};


template<Component c, std::size_t dim>
class Curl
{
    static constexpr auto i = detail::component_index<c>();
    static constexpr auto j = (i + 1) % NBR_COMPO;
    static constexpr auto k = (i + 2) % NBR_COMPO;

public:
    template<typename GridLayout, typename VecField>
    Curl(GridLayout const& layout, VecField const& v, std::array<QtyCentering, dim> const& dst)
        : dj_vk_{layout, v[k], dst}
        , dk_vj_{layout, v[j], dst}
    {
    }

    NO_DISCARD double operator()(MeshIndex<dim> const idx) const
    {
        if constexpr (j >= dim && k < dim)
            return -dk_vj_(idx);
        else
            return dj_vk_(idx) - dk_vj_(idx);
    }

private:
    Partial<static_cast<Direction>(j), dim> dj_vk_;
    Partial<static_cast<Direction>(k), dim> dk_vj_;
};


template<std::size_t dim>
class Div
{
public:
    template<typename GridLayout, typename VecField>
    Div(GridLayout const& layout, VecField const& v, std::array<QtyCentering, dim> const& dst)
        : dx_vx_{layout, v[0], dst}
        , dy_vy_{layout, v[1], dst}
        , dz_vz_{layout, v[2], dst}
    {
    }

    NO_DISCARD double operator()(MeshIndex<dim> const idx) const
    {
        return dx_vx_(idx) + dy_vy_(idx) + dz_vz_(idx);
    }

private:
    Partial<Direction::X, dim> dx_vx_;
    Partial<Direction::Y, dim> dy_vy_;
    Partial<Direction::Z, dim> dz_vz_;
};



template<auto dir, typename GridLayout, typename Src, typename Dst>
auto partial(GridLayout const& layout, Src const& src, Dst const& dst)
{
    return Partial<dir, GridLayout::dimension>{layout, src,
                                               detail::destination_centering(layout, dst)};
}

template<Component c, typename GridLayout, typename VecField, typename Dst>
auto curl(GridLayout const& layout, VecField const& v, Dst const& dst)
{
    return Curl<c, GridLayout::dimension>{layout, v, detail::destination_centering(layout, dst)};
}

template<Component c, typename GridLayout, typename Src, typename Dst>
auto grad(GridLayout const& layout, Src const& src, Dst const& dst)
{
    return partial<static_cast<Direction>(detail::component_index<c>())>(layout, src, dst);
}

template<typename GridLayout, typename VecField, typename Dst>
auto div(GridLayout const& layout, VecField const& v, Dst const& dst)
{
    return Div<GridLayout::dimension>{layout, v, detail::destination_centering(layout, dst)};
}

template<typename GridLayout, typename Src>
NO_DISCARD double laplacian(GridLayout const& layout, Src const& f,
                            MeshIndex<GridLayout::dimension> const idx)
{
    static_assert(Src::dimension == GridLayout::dimension,
                  "field and layout dimensions must match");

    return sum(for_N_make_array<GridLayout::dimension>([&](auto ic) {
        constexpr auto d = ic();
        auto const inv   = layout.inverseMeshSize(static_cast<Direction>(d));
        return inv * inv
               * (detail::shifted<d>(f, idx, 1) - 2.0 * detail::shifted<d>(f, idx, 0)
                  + detail::shifted<d>(f, idx, -1));
    }));
}

} // namespace PHARE::core

#endif
