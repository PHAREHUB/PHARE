#ifndef PHARE_TEST_CORE_DATA_TEST_VECFIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_VECFIELD_FIXTURES_HPP

#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"

namespace PHARE::core
{

template<std::size_t dim>
using VecField_t = VecField<Field_t<dim>, HybridQuantity>;

template<std::size_t dim>
using UsableVecField = UsableTensorField<dim, /*rank=*/1>;


} // namespace PHARE::core

#endif /*PHARE_TEST_CORE_DATA_TEST_VECFIELD_FIXTURES_HPP*/
