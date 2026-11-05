#ifndef PHARE_CORE_DATA_FIELD_FIELD_SCRATCH_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_SCRATCH_HPP

#include "core/data/vector.hpp"

namespace PHARE::core
{
// The one scratch buffer for temporary field data, rather than static vectors per class or
// function. Defined in field_scratch.cpp so there is one instance per linked phare_core, not one
// per template instantiation: one per process if phare_core is shared, otherwise one per module
// linking it, of which only the running simulator's is used. One user at a time: not for nested
// use.
MinimizingVector<double>& field_scratch();

void release_field_scratch() noexcept;

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_FIELD_FIELD_SCRATCH_HPP */
