
#ifndef PHARE_RESTART_MANAGER_HPP_
#define PHARE_RESTART_MANAGER_HPP_

#include "core/logger.hpp"
#include "core/data/particles/particle_array.hpp"

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
    virtual void dump(double timeStamp, double timeStep) = 0;
    inline virtual ~IRestartsManager();
};
IRestartsManager::~IRestartsManager() {}



template<typename Writer>
class RestartsManager : public IRestartsManager
{
public:
    void dump(double timeStamp, double timeStep) override;



    RestartsManager(std::unique_ptr<Writer>&& writer_ptr)
        : writer_{std::move(writer_ptr)}
    {
        if (!writer_)
            throw std::runtime_error("Error: RestartsManager received null Writer");
    }


    template<typename Hierarchy, typename Model>
    static std::unique_ptr<RestartsManager> make_unique(Hierarchy& hier, Model& model,
                                                        initializer::PHAREDict const& dict)
    {
        auto rMan = std::make_unique<RestartsManager>(Writer::make_unique(hier, model, dict));
        if (dict.contains("write_timestamps")) // else is only loading not saving restarts
            rMan->addRestartDict(dict);
        return rMan;
    }



    RestartsManager& addRestartDict(initializer::PHAREDict const& dict);
    RestartsManager& addRestartDict(initializer::PHAREDict&& dict) { return addRestartDict(dict); }


    Writer& writer() { return *writer_.get(); }


    RestartsManager(RestartsManager const&) = delete;
    RestartsManager(RestartsManager&&)      = delete;
    RestartsManager& operator=(RestartsManager const&) = delete;
    RestartsManager& operator=(RestartsManager&&) = delete;

private:
    bool needsAction_(double nextTime, double timeStamp, double timeStep)
    {
        // casting to float to truncate double to avoid trailing imprecision
        return static_cast<float>(std::abs(nextTime - timeStamp)) < static_cast<float>(timeStep);
    }


    bool needsWrite_(RestartsProperties const& rest, double const timeStamp, double const timeStep)
    {
        auto const& nextWrite = nextWrite_;

        return nextWrite < rest.writeTimestamps.size()
               and needsAction_(rest.writeTimestamps[nextWrite], timeStamp, timeStep);
    }


    std::unique_ptr<RestartsProperties> restarts_properties_;
    std::unique_ptr<Writer> writer_;
    std::size_t nextWrite_ = 0;
};



template<typename Writer>
RestartsManager<Writer>&
RestartsManager<Writer>::addRestartDict(initializer::PHAREDict const& params)
{
    if (restarts_properties_)
        throw std::runtime_error("Restarts invalid, properties already set");

    restarts_properties_ = std::make_unique<RestartsProperties>();
    restarts_properties_->writeTimestamps
        = params["write_timestamps"].template to<std::vector<double>>();

    return *this;
}




template<typename Writer>
void RestartsManager<Writer>::dump(double timeStamp, double timeStep)
{
    if (!restarts_properties_)
        return; // not active

    if (needsWrite_(*restarts_properties_, timeStamp, timeStep))
    {
        writer_->dump(timeStamp);
        ++nextWrite_;
    }
}


} // namespace PHARE::restarts


#endif /* PHARE_RESTART_MANAGER_HPP_ */
