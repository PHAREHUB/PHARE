#ifndef PHARE_TEST_CORE_DATA_TEST_VECFIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_VECFIELD_FIXTURES_HPP

#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"

namespace PHARE::core
{
template<std::size_t dim>
class UsableVecField : public UsableTensorField<dim, /*rank=*/1>
{
public:
    auto static constexpr dimension = dim;
    using Super                     = UsableTensorField<dim, /*rank=*/1>;

    template<typename GridLayout>
    UsableVecField(std::string const& name, GridLayout const& layout, HybridQuantity::Vector qty)
        : Super{name, layout, qty}
    {
    }
};

} // namespace PHARE::core

#endif /*PHARE_TEST_CORE_DATA_TEST_VECFIELD_FIXTURES_HPP*/
