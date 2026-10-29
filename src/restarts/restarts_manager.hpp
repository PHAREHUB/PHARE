#ifndef PHARE_RESTART_MANAGER_HPP_
#define PHARE_RESTART_MANAGER_HPP_


#include "core/def.hpp"
#include "core/logger.hpp"

#include "initializer/data_provider.hpp"

#include "mpi/mpi_utils.hpp"

#include "restarts/restarts_props.hpp"

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
        auto restarts_are_written = core::any(core::generate_from(
            [&](auto const& v) { return dict.contains(v); },
            std::vector<std::string>{"write_timestamps", "elapsed_timestamps"}));
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
    // Returns true if nextTime < timeStamp + timeStep
    NO_DISCARD bool needsCadenceAction_(double const nextTime, double const timeStamp,
                                        double const timeStep) const
    {
        // casting to float to truncate double to avoid trailing imprecision
        return static_cast<float>(nextTime - timeStamp) < static_cast<float>(timeStep);
    }

    // Advances id such as times[id] is past timeStamp + timeStep. If there is no such time, id =
    // times.size() on return; returns true if any encountered id needs action.
    NO_DISCARD bool needsCadenceActionAndAdvance_(std::vector<double> const& times, std::size_t& id,
                                                  double const timeStamp,
                                                  double const timeStep) const
    {
        bool acted = false;
        while (id < times.size() and needsCadenceAction_(times[id], timeStamp, timeStep))
        {
            acted = true;
            ++id;
        }
        return acted;
    }

    NO_DISCARD bool needsElapsedAction_(double const nextTime) const
    {
        return mpi::unix_timestamp_now() > nextTime;
    }


    NO_DISCARD bool needsWrite_(RestartsProperties const& rest, double const timeStamp,
                                double const timeStep)
    {
        auto const simUnit = needsCadenceActionAndAdvance_(rest.writeTimestamps, nextWriteSimUnit_,
                                                           timeStamp, timeStep);

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

    std::time_t const start_time_{mpi::unix_timestamp_now()};
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

    assert(params.contains("resources_hash"));

    restarts_properties_->fileAttributes["resources_hash"]
        = params["resources_hash"].template to<std::string>();

    return *this;
}




template<typename Writer>
bool RestartsManager<Writer>::dump(double timeStamp, double timeStep)
{
    if (!restarts_properties_)
        return false; // not active

    if (!needsWrite_(*restarts_properties_, timeStamp, timeStep))
        return false; // not needed now

    PHARE_LOG_SCOPE(3, "RestartsManager::dump");
    writer_->dump(*restarts_properties_, timeStamp);
    return true;
}


} // namespace PHARE::restarts


#endif /* PHARE_RESTART_MANAGER_HPP_ */
