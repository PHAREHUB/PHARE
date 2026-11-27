#ifndef PHARE_SYNCHRONIZER_HPP
#define PHARE_SYNCHRONIZER_HPP

#include "amr/messengers/communicator.hpp"

#include <SAMRAI/hier/CoarsenOperator.h>
#include <SAMRAI/hier/PatchHierarchy.h>
#include <SAMRAI/hier/PatchLevel.h>

#include <memory>
#include <stdexcept>
#include <string>

namespace PHARE::amr
{
template<typename ResourcesManager>
class Synchronizer : private Communicator<SynchronizerTypes, ResourcesManager::dimension>
{
public:
    Synchronizer(std::string const& name, std::shared_ptr<ResourcesManager> const& rm,
                 std::shared_ptr<SAMRAI::hier::CoarsenOperator> coarsenOp)
    {
        auto id = rm->getID(name);
        if (id)
        {
            this->add_algorithm()->registerCoarsen(*id, *id, coarsenOp);
        }
    }



    void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                       std::shared_ptr<SAMRAI::hier::PatchLevel> const& level)
    {
        auto levelNumber        = level->getLevelNumber();
        auto const& coarseLevel = hierarchy->getPatchLevel(levelNumber - 1);

        for (auto& algo : this->algos)
            this->add(algo, algo->createSchedule(coarseLevel, level), levelNumber);
    }


    void removeLevel(int const levelNumber)
    {
        Communicator<SynchronizerTypes, ResourcesManager::dimension>::removeLevel(levelNumber);
    }


    void sync(int const levelNumber)
    {
        if (this->algos.size() == 0)
            throw std::runtime_error("Algorithms are not configured");

        for (auto const& algo : this->algos)
            this->findSchedule(algo, levelNumber)->coarsenData();
    }
};
} // namespace PHARE::amr

#endif
