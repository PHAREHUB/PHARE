#ifndef DIAGNOSTIC_DAO_H
#define DIAGNOSTIC_DAO_H

#include <cstddef>
#include <string>
#include <vector>


namespace PHARE::diagnostic
{
struct DiagnosticDAO // DAO = DataAccessObject
{
    size_t lastWrite = 0, lastCompute = 0;
    std::vector<double> writeTimestamps, computeTimestamps;
    std::string category, type;
};

} // namespace PHARE::diagnostic

#endif // DIAGNOSTIC_DAO_H
