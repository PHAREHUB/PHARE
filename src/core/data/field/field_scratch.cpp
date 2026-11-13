#include "field_scratch.hpp"

namespace PHARE::core
{
MinimizingVector<double>& field_scratch()
{
    static MinimizingVector<double> buffer;
    return buffer;
}

int& field_scratch_holders()
{
    static int holders = 0;
    return holders;
}

void release_field_scratch() noexcept
{
    field_scratch().destroy();
}

} // namespace PHARE::core
