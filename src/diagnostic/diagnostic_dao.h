#ifndef DIAGNOSTIC_DAO_H
#define DIAGNOSTIC_DAO_H

#include <cstddef>
#include <string>
#include <vector>


namespace PHARE::diagnostic
{
struct DiagnosticDAO // DAO = DataAccessObject
{
    std::size_t compute_every = 1, write_every = 1;
    std::size_t start_iteration = 0, last_iteration = 100; /* likely to be time rather than index*/
                                                           /* do we allow ranges?*/
    size_t iterator;
    std::vector<size_t> activeTimesteps;
    std::string type, subtype;
};

} // namespace PHARE::diagnostic

#endif // DIAGNOSTIC_DAO_H
