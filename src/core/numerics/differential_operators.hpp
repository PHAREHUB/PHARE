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
#include <stdexcept>

namespace PHARE::core::detail
{
template<std::size_t d, typename Array, std::size_t dim>
double shifted(Array const& array, MeshIndex<dim> const& idx, int const shift)
{
    return for_N_fold<dim>([&](auto... i) { return array((i == d ? idx[i] + shift : idx[i])...); });
}

template<Component c>
constexpr std::size_t vector_component_index()
{
    static_assert(c == Component::X or c == Component::Y or c == Component::Z);
    return static_cast<std::size_t>(c);
}

template<auto dir, std::size_t dim>
constexpr bool direction_exists = static_cast<std::size_t>(dir) < dim;

constexpr QtyCentering flipped(QtyCentering const centering)
{
    return centering == QtyCentering::primal ? QtyCentering::dual : QtyCentering::primal;
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
    Partial(GridLayout const& layout, Src const& src)
        : view_{src.data(), src.shape()}
    {
        if (src.data() == nullptr)
            throw std::runtime_error("partial: the source has no data");

        auto const srcCentering = layout.centering(src);
        shift_                  = srcCentering[d] == QtyCentering::primal ? 0 : -1;
        scale_                  = layout.inverseMeshSize(dir);
        centering_              = srcCentering;
        centering_[d]           = detail::flipped(srcCentering[d]);
    }

    NO_DISCARD std::array<QtyCentering, dim> centering() const { return centering_; }

    NO_DISCARD double operator()(MeshIndex<dim> const idx) const
    {
        return scale_
               * (detail::shifted<d>(view_, idx, shift_ + 1)
                  - detail::shifted<d>(view_, idx, shift_));
    }

private:
    NdArrayView<dim, double const> view_;
    int shift_;
    double scale_;
    std::array<QtyCentering, dim> centering_;
};

template<auto dir, std::size_t dim>
    requires(not detail::direction_exists<dir, dim>)
class Partial<dir, dim>
{
public:
    template<typename GridLayout, typename Src>
    Partial(GridLayout const& layout, Src const& src)
        : centering_{layout.centering(src)}
    {
    }

    NO_DISCARD std::array<QtyCentering, dim> centering() const { return centering_; }

    NO_DISCARD double operator()(MeshIndex<dim> const) const { return 0.; }

private:
    std::array<QtyCentering, dim> centering_;
};


namespace detail
{
    template<Component c, std::size_t dim>
    class CurlComponent
    {
        static constexpr auto i = vector_component_index<c>();
        static constexpr auto j = (i + 1) % NBR_COMPO;
        static constexpr auto k = (i + 2) % NBR_COMPO;

    public:
        template<typename GridLayout, typename VecField>
        CurlComponent(GridLayout const& layout, VecField const& v)
            : dj_vk_{layout, v[k]}
            , dk_vj_{layout, v[j]}
        {
            if (dj_vk_.centering() != dk_vj_.centering())
                throw std::runtime_error("curl: the two derivatives land on different centerings");
        }

        NO_DISCARD std::array<QtyCentering, dim> centering() const { return dj_vk_.centering(); }

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
} // namespace detail


template<std::size_t dim>
struct Curl
{
    detail::CurlComponent<Component::X, dim> x;
    detail::CurlComponent<Component::Y, dim> y;
    detail::CurlComponent<Component::Z, dim> z;
};


template<std::size_t dim>
struct Grad
{
    Partial<Direction::X, dim> x;
    Partial<Direction::Y, dim> y;
    Partial<Direction::Z, dim> z;
};


template<std::size_t dim>
class Div
{
public:
    template<typename GridLayout, typename VecField>
    Div(GridLayout const& layout, VecField const& v)
        : dx_vx_{layout, v[0]}
        , dy_vy_{layout, v[1]}
        , dz_vz_{layout, v[2]}
    {
        if (dx_vx_.centering() != dy_vy_.centering() or dx_vx_.centering() != dz_vz_.centering())
            throw std::runtime_error("div: the three derivatives land on different centerings");
    }

    NO_DISCARD std::array<QtyCentering, dim> centering() const { return dx_vx_.centering(); }

    NO_DISCARD double operator()(MeshIndex<dim> const idx) const
    {
        return dx_vx_(idx) + dy_vy_(idx) + dz_vz_(idx);
    }

private:
    Partial<Direction::X, dim> dx_vx_;
    Partial<Direction::Y, dim> dy_vy_;
    Partial<Direction::Z, dim> dz_vz_;
};



template<auto dir, typename GridLayout, typename Src>
auto partial(GridLayout const& layout, Src const& src)
{
    return Partial<dir, GridLayout::dimension>{layout, src};
}

template<typename GridLayout, typename VecField>
auto curl(GridLayout const& layout, VecField const& v)
{
    return Curl<GridLayout::dimension>{{layout, v}, {layout, v}, {layout, v}};
}

template<typename GridLayout, typename Src>
auto grad(GridLayout const& layout, Src const& src)
{
    return Grad<GridLayout::dimension>{{layout, src}, {layout, src}, {layout, src}};
}

template<typename GridLayout, typename VecField>
auto div(GridLayout const& layout, VecField const& v)
{
    return Div<GridLayout::dimension>{layout, v};
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
