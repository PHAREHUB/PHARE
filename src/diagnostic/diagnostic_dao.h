#ifndef DIAGNOSTIC_DAO_H
#define DIAGNOSTIC_DAO_H

#include <cstddef>
#include <string>
#include <vector>


namespace PHARE::diagnostic
{
struct DiagnosticProperties
{
    size_t lastWrite = 0, lastCompute = 0;
    std::vector<double> writeTimestamps, computeTimestamps;
    std::string type, quantity;

    bool needsWrite(double timeStamp, double timeStep) const
    {
        return writeTimestamps[lastWrite] + timeStep > timeStamp;
    }

    bool needsCompute(double timeStamp, double timeStep) const
    {
        return computeTimestamps[lastCompute] + timeStep > timeStamp;
    }
};

} // namespace PHARE::diagnostic

#endif // DIAGNOSTIC_DAO_H
