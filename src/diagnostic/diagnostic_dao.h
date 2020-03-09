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
    double start_iteration = 0, last_iteration = 100;

    size_t lastWrite = 0, lastCompute = 0;
    std::vector<double> writeTimestamps, computeTimestamps;
    std::string category, type;

    void compute_active_timestamps(double time_step)
    {
        auto if_active_add = [&](auto& action, auto& every, auto& timestamps, double current) {
            action++;
            if (action % every == 0)
            {
                timestamps.push_back(current);
                action = 0;
            }
        };

        size_t write = 0, compute = 0;

        for (double curr = start_iteration; curr < last_iteration; curr += time_step)
        {
            if_active_add(compute, compute_every, computeTimestamps, curr);
            if_active_add(write, write_every, writeTimestamps, curr);
        }
    }
};

} // namespace PHARE::diagnostic

#endif // DIAGNOSTIC_DAO_H
