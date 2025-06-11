#ifndef PHARE_AMR_UTILITIES_BOX_BOX_HPP
#define PHARE_AMR_UTILITIES_BOX_BOX_HPP


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep


#include "SAMRAI/hier/Box.h"
#include "core/utilities/box/box.hpp"
#include "core/def.hpp"


namespace PHARE::amr
{
template<typename Type, std::size_t dim>
NO_DISCARD auto samrai_box_from(PHARE::core::Box<Type, dim> const& box, int samrai_blockId = 0)
{
    SAMRAI::tbox::Dimension dimension{dim};
    SAMRAI::hier::BlockId blockId{samrai_blockId};
    return SAMRAI::hier::Box{SAMRAI::hier::Index{dimension, (*box.lower).data()},
                             SAMRAI::hier::Index{dimension, (*box.upper).data()}, blockId};
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


template<typename Particle>
NO_DISCARD inline bool isInBox(SAMRAI::hier::Box const& box, Particle const& particle)
{
    constexpr auto dim = Particle::dimension;
    auto const& iCell  = particle.iCell;
    auto const& lower  = box.lower();
    auto const& upper  = box.upper();
    for (std::size_t i = 0; i < dim; ++i)
        if (iCell[i] < lower(i) || iCell[i] > upper(i))
            return false;
    return true;
}


template<std::size_t dim>
auto as_point(SAMRAI::hier::IntVector const& vec)
{
    return core::Point{
        core::for_N<dim, core::for_N_R_mode::make_array>([&](auto i) { return vec[i]; })};
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

} // namespace PHARE::amr

#endif
