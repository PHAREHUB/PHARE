#ifndef PHARE_AMR_UTILITIES_BOX_BOX_HPP
#define PHARE_AMR_UTILITIES_BOX_BOX_HPP


#include "core/def.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"

#include "amr/amr_constants.hpp"

#include "SAMRAI/hier/Box.h"


namespace PHARE::amr
{
template<typename Type, std::size_t dim>
NO_DISCARD auto samrai_box_from(PHARE::core::Box<Type, dim> const& box, int samrai_blockId = 0)
{
    SAMRAI::tbox::Dimension dimension{dim};
    SAMRAI::hier::BlockId blockId{samrai_blockId};
    return SAMRAI::hier::Box{SAMRAI::hier::Index{dimension, &box.lower[0]},
                             SAMRAI::hier::Index{dimension, &box.upper[0]}, blockId};
}

template<std::size_t dim, typename Type = int>
NO_DISCARD auto phare_box_from(SAMRAI::hier::Box const& box)
{
    std::array<Type, dim> lower = *reinterpret_cast<std::array<int, dim> const*>(&box.lower()[0]);
    std::array<Type, dim> upper = *reinterpret_cast<std::array<int, dim> const*>(&box.upper()[0]);

    return PHARE::core::Box<Type, dim>{core::Point{lower}, core::Point{upper}};
}

template<std::size_t dim>
NO_DISCARD auto as_unsigned_phare_box(SAMRAI::hier::Box const& box)
{
    auto const& amr_box = phare_box_from<dim>(box);
    return PHARE::core::Box<std::uint32_t, dim>{core::Point{amr_box.lower}.as_unsigned(),
                                                core::Point{amr_box.upper}.as_unsigned()};
}

NO_DISCARD inline bool operator==(SAMRAI::hier::Box const& b1, SAMRAI::hier::Box const& b2)
{
    auto dim1 = b1.getDim().getValue();
    auto dim2 = b2.getDim().getValue();

    if (dim1 != dim2)
        return false;

    bool boxesAreEqual = true;
    for (auto i = 0u; i < dim1; ++i)
    {
        boxesAreEqual &= b1.lower(i) == b2.lower(i);
        boxesAreEqual &= b1.upper(i) == b2.upper(i);
    }
    return boxesAreEqual;
}

template<typename Type, std::size_t dim>
struct Box : public PHARE::core::Box<Type, dim>
{
    using Super = PHARE::core::Box<Type, dim>;
    using Super::shape;
    using Super::size;

    using Super::lower;
    using Super::upper;

    Box() = default;

    Box(std::array<Type, dim> _lower, std::array<Type, dim> _upper)
        : Super{core::Point{_lower}, core::Point{_upper}}
    {
    }

    template<typename T, std::size_t s>
    Box(core::Point<T, s> _lower, core::Point<T, s> _upper)
        : Super{_lower, _upper}
    {
    }

    Box(SAMRAI::hier::Box const& box)
        : Super{phare_box_from<dim, Type>(box)}
    {
    }

    NO_DISCARD operator SAMRAI::hier::Box() const { return samrai_box_from(*this); }

    NO_DISCARD bool operator==(SAMRAI::hier::Box const& that) const
    {
        bool eq = 1;

        for (std::size_t i = 0u; i < dim; ++i)
            eq &= (this->lower[i] == that.lower()[i]) and (this->upper[i] == that.upper()[i]);

        return eq;
    }
};


template<typename T, std::size_t dim>
inline bool isInBox(SAMRAI::hier::Box const& box, std::array<T, dim> const& iCell)
{
    auto const& lower = box.lower();
    auto const& upper = box.upper();
    for (std::size_t i = 0; i < dim; ++i)
        if (iCell[i] < lower(i) || iCell[i] > upper(i))
            return false;
    return true;
}


template<typename Particle>
inline bool isInBox(SAMRAI::hier::Box const& box, Particle const& particle)
{
    return isInBox(box, particle.iCell());
}


template<std::size_t dim>
auto as_point(SAMRAI::hier::IntVector const& vec)
{
    return core::Point{core::for_N_make_array<dim>([&](auto i) { return vec[i]; })};
}


template<std::size_t dim>
auto as_point(SAMRAI::hier::Transformation const& tform)
{
    return as_point<dim>(tform.getOffset());
}



template<typename Type, std::size_t dim>
NO_DISCARD core::Box<Type, dim> shift(core::Box<Type, dim> const& box,
                                      SAMRAI::hier::Transformation const& tform)
{
    return core::shift(box, as_point<dim>(tform));
}


template<typename Box_t>
Box_t refine_box(Box_t const& box)
{
    return Box_t{((box.lower) * refinementRatio) + 1, ((box.upper) * refinementRatio) + 1};
}


template<typename Box_t>
Box_t coarsen_box(Box_t const& box)
{
    auto static constexpr dim = Box_t::dimension;

    return Box_t{core::Point{core::for_N_make_array<dim>(
                     [&](auto i) { return box.lower[i] / refinementRatio; })},
                 core::Point{core::for_N_make_array<dim>(
                     [&](auto i) { return box.upper[i] / refinementRatio; })}};
}


} // namespace PHARE::amr

#endif
