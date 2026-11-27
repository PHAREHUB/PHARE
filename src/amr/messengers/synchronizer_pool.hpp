#ifndef PHARE_SYNCHRONIZER_POOL_HPP
#define PHARE_SYNCHRONIZER_POOL_HPP

#include "amr/messengers/synchronizer.hpp"

#include <SAMRAI/hier/CoarsenOperator.h>
#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/hier/PatchLevel.h>

#include <map>
#include <memory>
#include <stdexcept>
#include <string>



namespace PHARE::amr
{


template<typename ResourcesManager>
class SynchronizerPool
{
public:
    template<typename Descriptor>
    void add(Descriptor const& descriptor,
             std::shared_ptr<SAMRAI::hier::CoarsenOperator> const& coarsenOp, std::string key)
    {
        // check if key is in synchronizers_
        auto const [it, success] = synchronizers_.insert(
            {key, Synchronizer<ResourcesManager>(descriptor, rm_, coarsenOp)});
        if (!success)
            throw std::runtime_error(key + " is already registered");
    }


    void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                       std::shared_ptr<SAMRAI::hier::PatchLevel> const& level)
    {
        for (auto& [_, synchronizer] : synchronizers_)
            synchronizer.registerLevel(hierarchy, level);
    }



    void removeLevel(int const levelNumber)
    {
        for (auto& [_, synchronizer] : synchronizers_)
            synchronizer.removeLevel(levelNumber);
    }



    void sync(int const levelNumber)
    {
        for (auto& [key, synchronizer] : synchronizers_)
        {
            synchronizer.sync(levelNumber);
        }
    }

    SynchronizerPool(std::shared_ptr<ResourcesManager> const& rm)
        : rm_{rm}
    {
    }

private:
    std::shared_ptr<ResourcesManager> rm_;
    std::map<std::string, Synchronizer<ResourcesManager>> synchronizers_;
};


} // namespace PHARE::amr
#endif
