#ifndef PHARE_RESTART_MANAGER_HPP_
#define PHARE_RESTART_MANAGER_HPP_


#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/utilities/mpi_utils.hpp"

#include "initializer/data_provider.hpp"

#include "restarts_props.hpp"


#include <cmath>
#include <memory>
#include <utility>



namespace PHARE::restarts
{
class IRestartsManager
{
public:
    virtual bool dump(double timeStamp, double timeStep) = 0;
    inline virtual ~IRestartsManager();
};
IRestartsManager::~IRestartsManager() {}



template<typename Writer>
class RestartsManager : public IRestartsManager
{
public:
    bool dump(double timeStamp, double timeStep) override;



    RestartsManager(std::unique_ptr<Writer>&& writer_ptr)
        : writer_{std::move(writer_ptr)}
    {
        if (!writer_)
            throw std::runtime_error("Error: RestartsManager received null Writer");
    }


    template<typename Hierarchy, typename ResourceManager_t>
    NO_DISCARD static std::unique_ptr<RestartsManager>
    make_unique(Hierarchy& hier, ResourceManager_t& resman, initializer::PHAREDict const& dict)
    {
        auto rMan = std::make_unique<RestartsManager>(Writer::make_unique(hier, resman, dict));
        auto restarts_are_written = core::any(
            core::generate([&](auto const& v) { return dict.contains(v); },
                           std::vector<std::string>{"write_timestamps", "elapsed_timestamps"}));
        // niter_period leaves no timestamp arrays; activate on a non-zero iteration cadence too.
        if (dict.contains("write_niter_period")
            and dict["write_niter_period"].template to<std::size_t>() > 0)
            restarts_are_written = true;
        if (restarts_are_written) // else is only loading not saving restarts
            rMan->addRestartDict(dict);
        return rMan;
    }



    RestartsManager& addRestartDict(initializer::PHAREDict const& dict);
    RestartsManager& addRestartDict(initializer::PHAREDict&& dict) { return addRestartDict(dict); }


    NO_DISCARD Writer& writer() { return *writer_.get(); }


    RestartsManager(RestartsManager const&)            = delete;
    RestartsManager(RestartsManager&&)                 = delete;
    RestartsManager& operator=(RestartsManager const&) = delete;
    RestartsManager& operator=(RestartsManager&&)      = delete;

private:
    // A scheduled time is "reached" if it lies before the next step boundary: scheduledTime is
    // within (or behind) the current step [timeStamp, timeStamp+timeStep). Compared as
    // (scheduledTime - timeStamp) < timeStep with a float cast to truncate trailing fp imprecision
    // (same robustness as the historical |nextTime-timeStamp| < timeStep): a time exactly one step
    // ahead (scheduledTime == timeStamp+timeStep) is NOT consumed now -- it belongs to the next
    // step. There is no abs(): a time already behind timeStamp yields a large negative difference,
    // so it stays "reached" and the catch-up loop keeps advancing instead of freezing.
    bool reached_(double const scheduledTime, double const timeStamp, double const timeStep) const
    {
        return static_cast<float>(scheduledTime - timeStamp) < static_cast<float>(timeStep);
    }

    // Advance idx past every scheduled time the current step has reached; return true if any.
    // The while-loop (vs a single ++) keeps the cadence from freezing when one step overshoots
    // several scheduled times (dt > period, or adaptive dt growing past the period): a single ++
    // would let nextTime fall more than a step behind currentTime, after which
    // |nextTime-currentTime| never drops below dt again and all remaining checkpoints are silently
    // lost.
    bool catchUp_(std::vector<double> const& times, std::size_t& idx, double const timeStamp,
                  double const timeStep) const
    {
        bool acted = false;
        while (idx < times.size() and reached_(times[idx], timeStamp, timeStep))
        {
            acted = true;
            ++idx;
        }
        return acted;
    }

    bool needsElapsedAction_(double const nextTime) const
    {
        return core::mpi::unix_timestamp_now() > nextTime;
    }


    bool needsWrite_(RestartsProperties const& rest, double const timeStamp, double const timeStep)
    {
        auto const simUnit = catchUp_(rest.writeTimestamps, nextWriteSimUnit_, timeStamp, timeStep);

        auto const elapsed = nextWriteElapsed_ < rest.elapsedTimestamps.size()
                             and needsElapsedAction_(rest.elapsedTimestamps[nextWriteElapsed_]);

        if (elapsed)
            ++nextWriteElapsed_;

        return simUnit || elapsed;
    }


    std::unique_ptr<RestartsProperties> restarts_properties_;
    std::unique_ptr<Writer> writer_;
    std::size_t nextWriteSimUnit_ = 0;
    std::size_t nextWriteElapsed_ = 0;
    std::size_t iteration_        = 0; ///< coarse-step counter for writeNiterPeriod cadence

    std::time_t const start_time_{core::mpi::unix_timestamp_now()};
};



template<typename Writer>
RestartsManager<Writer>&
RestartsManager<Writer>::addRestartDict(initializer::PHAREDict const& params)
{
    if (restarts_properties_)
        throw std::runtime_error("Restarts invalid, properties already set");

    restarts_properties_ = std::make_unique<RestartsProperties>();

    if (params.contains("write_timestamps"))
        restarts_properties_->writeTimestamps
            = params["write_timestamps"].template to<std::vector<double>>();

    if (params.contains("write_niter_period"))
        restarts_properties_->writeNiterPeriod
            = params["write_niter_period"].template to<std::size_t>();

    if (params.contains("elapsed_timestamps"))
    {
        restarts_properties_->elapsedTimestamps
            = params["elapsed_timestamps"].template to<std::vector<double>>();
        for (auto& time : restarts_properties_->elapsedTimestamps)
            time += start_time_; // expected for comparison later
    }

    assert(params.contains("serialized_simulation"));

    restarts_properties_->fileAttributes["serialized_simulation"]
        = params["serialized_simulation"].template to<std::string>();

    return *this;
}




template<typename Writer>
bool RestartsManager<Writer>::dump(double timeStamp, double timeStep)
{
    if (!restarts_properties_)
        return false; // not active

    // iteration-based cadence: write every writeNiterPeriod coarse steps (the only timestamp-free
    // option valid under adaptive dt). needsWrite_ is called unconditionally so its timestamp
    // index still advances.
    bool const niterNow = restarts_properties_->writeNiterPeriod > 0
                          and (iteration_ % restarts_properties_->writeNiterPeriod == 0);
    bool const writeNow = needsWrite_(*restarts_properties_, timeStamp, timeStep);
    ++iteration_;

    if (!(writeNow || niterNow))
        return false; // not needed now

    PHARE_LOG_SCOPE(3, "RestartsManager::dump");
    writer_->dump(*restarts_properties_, timeStamp);
    return true;
}


} // namespace PHARE::restarts


#endif /* PHARE_RESTART_MANAGER_HPP_ */
