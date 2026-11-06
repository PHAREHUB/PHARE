#ifndef PHARE_AMR_DATA_FIELD_CENTERED_FACTORY_HPP
#define PHARE_AMR_DATA_FIELD_CENTERED_FACTORY_HPP

#include "phare_mpi.hpp" // IWYU pragma: keep

#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/utilities/meta/enum_meta.hpp"

#include <SAMRAI/hier/PatchDataFactory.h>

#include <array>
#include <cstdint>
#include <memory>
#include <type_traits>
#include <variant>

/*
 * SAMRAI groups the items of a refine or coarsen algorithm by the type of their patch data factory
 * and their ghost width (RefineClasses::patchDataMatch, CoarsenClasses::patchDataMatch), and
 * computes the overlaps of a group once, from its first item. Field and tensor field factories are
 * therefore given a distinct type per centering, so that items of different centerings can share
 * an algorithm, hence a schedule.
 *
 * The split is by factory type rather than by fill pattern type because, when a fill needs a coarse
 * interpolation level, SAMRAI replaces every fill pattern with a BoxGeometryVariableFillPattern and
 * regroups the items (RefineSchedule.cpp, coarse interpolation schedule): classes split by fill
 * pattern type merge there, classes split by factory type do not.
 */
namespace PHARE::amr
{
// one bit per primal direction
template<std::size_t dim>
constexpr std::uint32_t centering_key(std::array<core::QtyCentering, dim> const& centering)
{
    std::uint32_t key = 0;
    for (std::size_t d = 0; d < dim; ++d)
        if (centering[d] == core::QtyCentering::primal)
            key |= 1u << d;
    return key;
}

// component keys packed dim bits apart
template<std::size_t dim, std::size_t n_components>
constexpr std::uint32_t
centering_key(std::array<std::array<core::QtyCentering, dim>, n_components> const& centerings)
{
    std::uint32_t key = 0;
    for (std::size_t c = 0; c < n_components; ++c)
        key |= centering_key(centerings[c]) << (c * dim);
    return key;
}

template<typename BaseFactory, std::uint32_t key>
class CenteredFactory final : public BaseFactory
{
public:
    using BaseFactory::BaseFactory;

    // the base clone returns a BaseFactory, which would drop the centering from the type
    std::shared_ptr<SAMRAI::hier::PatchDataFactory>
    cloneFactory(SAMRAI::hier::IntVector const& /*ghost*/) final
    {
        return std::make_shared<CenteredFactory>(this->fineBoundaryRepresentsVariable(),
                                                 this->dataLivesOnPatchBorder(), this->name_,
                                                 this->quantity_);
    }
};

/**
 * @brief makes the BaseFactory subtype of the centering of qty. The type depends on the centering
 * only: quantities of one centering share it, hence their equivalence classes.
 */
template<typename BaseFactory, typename GridLayoutT, typename Quantity>
std::shared_ptr<SAMRAI::hier::PatchDataFactory> make_centered_factory(Quantity const qty,
                                                                      auto&&... args)
{
    return std::visit(
        [&]<Quantity q>(std::integral_constant<Quantity, q>)
            -> std::shared_ptr<SAMRAI::hier::PatchDataFactory> {
            return std::make_shared<
                CenteredFactory<BaseFactory, centering_key(GridLayoutT::centering(q))>>(args...);
        },
        core::asEnumConstant(qty));
}

} // namespace PHARE::amr

#endif
