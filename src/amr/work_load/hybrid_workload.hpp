
#include <SAMRAI/hier/PatchLevel.h>

#include "core/data/ions/ions.h"


class HybridWorkLoadEstimator : public WorkLoadEstimator
{
    public :
        virtual void estimate(SAMRAI::hier::PatchLevel levels, double* wl, PHARE::core::Ions const& ions)
        {
            for (auto& p : levels)
            {
                auto pd = p.getPatchData(this->getID());
                (void)pd;



                // TODO



            }
        };
        virtual void set_strategy(std::string stratName)
        {
            strat_ = HybridWorkLoadStrategyFactory::create(stratName);
        };
    private:
        std::unique_ptr<HybridWorkLoadEstimatorStrategy> strat_;
};
