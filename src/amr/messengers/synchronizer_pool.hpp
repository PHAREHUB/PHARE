#ifndef PHARE_SYNCHRONIZER_POOL_HPP
#define PHARE_SYNCHRONIZER_POOL_HPP

#include "communicator.hpp"



namespace PHARE::amr
{


template<std::size_t dimension>
class SynchronizerPool
{
public:
    template<typename Descriptor, typename ResourcesManager>
    void add(Descriptor const& descriptor,
             std::shared_ptr<SAMRAI::hier::CoarsenOperator> const& coarsenOp, std::string key,
             std::shared_ptr<ResourcesManager> const& rm)
    {
        // check if key is in synchronizers_
        auto const [it, success] = synchronizers_.insert(
            {key, makeSynchronizer<ResourcesManager, dimension>(descriptor, rm, coarsenOp)});
        if (!success)
            throw std::runtime_error(key + " is already registered");
    }




    template<typename ResourcesManager>
    void add(VecFieldDescriptor const& descriptor, std::shared_ptr<ResourcesManager> const& rm,
             std::shared_ptr<SAMRAI::hier::CoarsenOperator> const& coarsenOp, std::string key)
    {
        auto const [it, success] = synchronizers_.insert(
            {key, makeSynchronizer<ResourcesManager, dimension>(descriptor, rm, coarsenOp)});

        if (!success)
            throw std::runtime_error(key + " is already registered");
    }



    void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                       std::shared_ptr<SAMRAI::hier::PatchLevel> const& level)
    {
        auto levelNumber        = level->getLevelNumber();
        auto const& coarseLevel = hierarchy->getPatchLevel(levelNumber - 1);

        for (auto& [_, synchronizer] : synchronizers_)
            for (auto& algo : synchronizer.algos)
                synchronizer.add(algo, algo->createSchedule(coarseLevel, level), levelNumber);
    }



    void sync(int const levelNumber)
    {
        for (auto& [key, synchronizer] : synchronizers_)
        {
            if (synchronizer.algos.size() == 0)
                throw std::runtime_error("Algorithms are not configured");

            for (auto const& algo : synchronizer.algos)
                synchronizer.findSchedule(algo, levelNumber)->coarsenData();
        }
    }

private:
    std::map<std::string, Communicator<Synchronizer, dimension>> synchronizers_;
};


} // namespace PHARE::amr
#endif
