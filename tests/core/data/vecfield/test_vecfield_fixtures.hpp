#ifndef PHARE_TEST_CORE_DATA_TEST_VECFIELD_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_VECFIELD_FIXTURES_HPP

#include "core/def/phare_config.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include "tests/core/data/field/test_field_fixtures.hpp"
#include "tests/core/data/tensorfield/test_tensorfield_fixtures.hpp"

namespace PHARE::core
{


template<auto opts>
using UsableVecField = UsableTensorField<opts.with_rank(1)>;


} // namespace PHARE::core

#endif /*PHARE_TEST_CORE_DATA_TEST_VECFIELD_FIXTURES_HPP*/
