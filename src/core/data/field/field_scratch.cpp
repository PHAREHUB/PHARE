#include "field_scratch.hpp"

namespace PHARE::core
{
MinimizingVector<double>& field_scratch()
{
    static MinimizingVector<double> buffer;
    return buffer;
}

void release_field_scratch() noexcept
{
    field_scratch().destroy();
}

} // namespace PHARE::core
