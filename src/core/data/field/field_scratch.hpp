#ifndef PHARE_CORE_DATA_FIELD_FIELD_SCRATCH_HPP
#define PHARE_CORE_DATA_FIELD_FIELD_SCRATCH_HPP

#include "core/data/vector.hpp"

#include <span>
#include <stdexcept>

namespace PHARE::core
{
// The one scratch buffer for temporary field data, rather than static vectors per class or
// function. Defined in field_scratch.cpp so there is one instance per linked phare_core, not one
// per template instantiation: one per process if phare_core is shared, otherwise one per module
// linking it, of which only the running simulator's is used. Access it through ScratchWrite or
// ScratchRead below, which count the live holders so nested use cannot silently overwrite it.
MinimizingVector<double>& field_scratch();

// number of live ScratchWrite/ScratchRead handles
int& field_scratch_holders();

void release_field_scratch() noexcept;


// exclusive: only when nobody else holds the buffer
class ScratchWrite
{
public:
    ScratchWrite()
    {
        if (field_scratch_holders() != 0)
            throw std::runtime_error("field_scratch: write while already held");
        ++field_scratch_holders();
    }
    ~ScratchWrite() { --field_scratch_holders(); }

    ScratchWrite(ScratchWrite const&)            = delete;
    ScratchWrite& operator=(ScratchWrite const&) = delete;

    MinimizingVector<double>& get() { return field_scratch(); }
};


// shared, read-only view: nested code may read the holder's data
class ScratchRead
{
public:
    ScratchRead() { ++field_scratch_holders(); }
    ~ScratchRead() { --field_scratch_holders(); }

    ScratchRead(ScratchRead const&)            = delete;
    ScratchRead& operator=(ScratchRead const&) = delete;

    std::span<double const> get() const { return field_scratch()(); }
};

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_FIELD_FIELD_SCRATCH_HPP */
