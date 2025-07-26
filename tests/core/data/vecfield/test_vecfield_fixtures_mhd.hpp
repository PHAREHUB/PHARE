#ifndef PHARE_TEST_CORE_DATA_TEST_VECFIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_VECFIELD_FIXTURES_HPP

#include "tests/core/data/field/test_field_fixtures_mhd.hpp"
#include "tests/core/data/tensorfield/test_tensorfield_fixtures_mhd.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/data/vecfield/vecfield.hpp"

namespace PHARE::core
{

template<std::size_t dim>
using VecFieldMHD = VecField<FieldMHD<dim>, MHDQuantity>;

template<std::size_t dim>
using UsableVecFieldMHD = UsableTensorFieldMHD<dim, /*rank=*/1>;


} // namespace PHARE::core

#endif /*PHARE_TEST_CORE_DATA_TEST_VECFIELD_FIXTURES_HPP*/
