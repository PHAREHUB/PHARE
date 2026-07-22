#ifndef PHARE_CORE_UTILITIES_CADENCE_HPP
#define PHARE_CORE_UTILITIES_CADENCE_HPP

#include <cstddef>
#include <vector>

#include "core/utilities/mpi_utils.hpp"

namespace PHARE::core
{
// Shared by DiagnosticsManager and RestartsManager: both schedule dumps/checkpoints against an
// array of absolute target times, using the same "reached + catch-up" semantics.

// A scheduled time is "reached" if it lies before the next step boundary: scheduledTime is
// within (or behind) the current step [timeStamp, timeStamp+timeStep). Compared as
// (scheduledTime - timeStamp) < timeStep with a float cast to truncate trailing fp imprecision:
// a time exactly one step ahead (scheduledTime == timeStamp+timeStep) is NOT consumed now -- it
// belongs to the next step. There is no abs(): a time already behind timeStamp yields a large
// negative difference, so it stays "reached" and the catch-up loop keeps advancing instead of
// freezing.
inline bool cadence_reached(double scheduledTime, double timeStamp, double timeStep)
{
    return static_cast<float>(scheduledTime - timeStamp) < static_cast<float>(timeStep);
}

// Advance idx past every scheduled time the current step has reached; return true if any.
// The while-loop (vs a single ++) keeps the cadence from freezing when one step overshoots
// several scheduled times (dt > period, or adaptive dt growing past the period): a single ++
// would let idx fall more than a step behind currentTime, after which cadence_reached() never
// triggers again and all remaining scheduled times are silently lost.
inline bool cadence_catch_up(std::vector<double> const& times, std::size_t& idx,
                             double timeStamp, double timeStep)
{
    bool acted = false;
    while (idx < times.size() and cadence_reached(times[idx], timeStamp, timeStep))
    {
        acted = true;
        ++idx;
    }
    return acted;
}

inline bool cadence_elapsed(double nextTime)
{
    return core::mpi::unix_timestamp_now() > nextTime;
}

} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_CADENCE_HPP
