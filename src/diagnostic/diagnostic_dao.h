#ifndef DIAGNOSTIC_DAO_H
#define DIAGNOSTIC_DAO_H

#include <cstddef>
#include <string>


namespace PHARE
{
namespace diagnostic
{
    struct DiagnosticDAO // DAO = DataAccessObject
    {
        std::size_t compute_every = 1, write_every = 1;
        std::size_t start_iteration = 0,
                    end_iteration   = 100; /* likely to be time rather than index*/
                                           /* do we allow ranges?*/
        std::string name, species, type;
    };

} // namespace diagnostic
} // namespace PHARE

#endif // DIAGNOSTIC_DAO_H
