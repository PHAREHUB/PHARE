#ifndef PHARE_TESTS_AMR_TEST_LEVEL_FIXTURES_HPP
#define PHARE_TESTS_AMR_TEST_LEVEL_FIXTURES_HPP

#include "core/utilities/box/box.hpp"
#include "amr/utilities/box/amr_box.hpp"

#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/hier/BoxLevel.h>
#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include <array>
#include <memory>
#include <vector>

namespace PHARE::test
{
// a hierarchy over `domain`, of cell width 1, whose level 0 is made of exactly `boxes`, all on
// this rank: same-level layouts the load balancer would rarely produce on demand
template<std::size_t dim>
auto hierarchyWithLevel0(core::Box<int, dim> const& domain,
                         std::vector<core::Box<int, dim>> const& boxes)
{
    SAMRAI::tbox::Dimension const samraiDim{dim};
    std::array<double, dim> lower{}, upper{};
    for (std::size_t i = 0; i < dim; ++i)
        upper[i] = domain.upper[i] - domain.lower[i] + 1;
    SAMRAI::hier::BoxContainer domainBoxes{amr::samrai_box_from(domain)};
    auto geometry = std::make_shared<SAMRAI::geom::CartesianGridGeometry>(
        "geometry", lower.data(), upper.data(), domainBoxes);
    auto hierarchy = std::make_shared<SAMRAI::hier::PatchHierarchy>("hierarchy", geometry);

    SAMRAI::hier::BoxLevel level0{SAMRAI::hier::IntVector::getOne(samraiDim), geometry};
    auto const rank = SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld().getRank();
    for (std::size_t i = 0; i < boxes.size(); ++i)
        level0.addBox(SAMRAI::hier::Box{amr::samrai_box_from(boxes[i]),
                                        SAMRAI::hier::LocalId{static_cast<int>(i)}, rank});
    hierarchy->makeNewPatchLevel(0, level0);
    return hierarchy;
}

} // namespace PHARE::test

#endif /* PHARE_TESTS_AMR_TEST_LEVEL_FIXTURES_HPP */
