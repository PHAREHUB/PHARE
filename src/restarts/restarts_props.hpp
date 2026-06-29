#ifndef RESTART_DAO_HPP
#define RESTART_DAO_HPP

#include <cstddef>
#include <vector>

#include "dict.hpp"

namespace PHARE::restarts
{
struct RestartsProperties
{
    using FileAttributes = cppdict::Dict<std::string>;

    std::vector<double> writeTimestamps, elapsedTimestamps;

    // write a restart every writeNiterPeriod coarse iterations (0 = disabled, use timestamps).
    // Iteration cadence is the only timestamp-free option valid under adaptive dt.
    std::size_t writeNiterPeriod = 0;

    FileAttributes fileAttributes{};
};

} // namespace PHARE::restarts

#endif // RESTART_DAO_HPP
