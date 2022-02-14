
#ifndef PHARE_AMR_UTILITIES_BOX_BOX_HPP
#define PHARE_AMR_UTILITIES_BOX_BOX_HPP


#include "SAMRAI/hier/Box.h"
#include "core/utilities/box/box.hpp"


namespace PHARE::amr
{
template<typename Type, std::size_t dim>
auto samrai_box_from(PHARE::core::Box<Type, dim> const& box, int samrai_blockId = 0)
{
    SAMRAI::tbox::Dimension dimension{dim};
    SAMRAI::hier::BlockId blockId{samrai_blockId};
    return SAMRAI::hier::Box{SAMRAI::hier::Index{dimension, &box.lower[0]},
                             SAMRAI::hier::Index{dimension, &box.upper[0]}, blockId};
}

template<typename Type, std::size_t dim>
auto phare_box_from(SAMRAI::hier::Box const& box)
{
    std::array<Type, dim> lower = *reinterpret_cast<std::array<int, dim> const*>(&box.lower()[0]);
    std::array<Type, dim> upper = *reinterpret_cast<std::array<int, dim> const*>(&box.upper()[0]);

    return PHARE::core::Box<Type, dim>{core::Point{lower}, core::Point{upper}};
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
        : Super{phare_box_from<Type, dim>(box)}
    {
    }

    operator SAMRAI::hier::Box() const { return samrai_box_from(*this); }

    bool operator==(SAMRAI::hier::Box const& that) const
    {
        bool eq = 1;

        for (std::size_t i = 0u; i < dim; ++i)
            eq &= (this->lower[i] == that.lower()[i]) and (this->upper[i] == that.upper()[i]);

        return eq;
    }
};




template<std::size_t dim>
inline bool isInBox(SAMRAI::hier::Box const& box, std::array<int, dim> const& iCell)
{
    auto const& lower = box.lower();
    auto const& upper = box.upper();

    if (iCell[0] >= lower(0) && iCell[0] <= upper(0))
    {
        if constexpr (dim > 1)
        {
            if (iCell[1] >= lower(1) && iCell[1] <= upper(1))
            {
                if constexpr (dim > 2)
                {
                    if (iCell[2] >= lower(2) && iCell[2] <= upper(2))
                    {
                        return true;
                    }
                }
                else
                {
                    return true;
                }
            }
        }
        else
        {
            return true;
        }
    }
    return false;
}


template<typename Particle>
inline bool isInBox(SAMRAI::hier::Box const& box, Particle const& particle)
{
    return isInBox(box, particle.iCell());
}


} // namespace PHARE::amr

#endif
