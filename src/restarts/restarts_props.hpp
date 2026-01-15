#ifndef RESTART_DAO_HPP
#define RESTART_DAO_HPP

#include <vector>

#include "dict.hpp"

namespace PHARE::restarts
{
struct RestartsProperties
{
    using FileAttributes = cppdict::Dict<std::string>;

    std::vector<double> writeTimestamps, elapsedTimestamps;

    FileAttributes fileAttributes{};
};

} // namespace PHARE::restarts

#endif // RESTART_DAO_HPP
