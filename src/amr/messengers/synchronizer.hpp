#ifndef PHARE_SYNCHRONIZER_HPP
#define PHARE_SYNCHRONIZER_HPP

#include "communicator.hpp"
#include "core/data/vecfield/vecfield.hpp"

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
