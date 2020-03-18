#ifndef DIAGNOSTIC_DAO_H
#define DIAGNOSTIC_DAO_H

#include <cstddef>
#include <string>
#include <vector>


namespace PHARE::diagnostic
{
struct DiagnosticDAO // DAO = DataAccessObject
{
    double start_iteration = 0, last_iteration = 100;

    size_t lastWrite = 0, lastCompute = 0;
    std::vector<double> writeTimestamps, computeTimestamps;
    std::string category, type;
};

} // namespace PHARE::diagnostic

#endif // DIAGNOSTIC_DAO_H
