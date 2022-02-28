#ifndef RESTART_DAO_HPP
#define RESTART_DAO_HPP

#include <string>
#include <vector>
#include <cstddef>

#include "cppdict/include/dict.hpp"

namespace PHARE::restarts
{
struct RestartsProperties
{
    std::vector<double> writeTimestamps;
};

} // namespace PHARE::restarts

#endif // RESTART_DAO_HPP
