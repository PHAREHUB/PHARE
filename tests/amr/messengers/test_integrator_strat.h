#ifndef PHARE_TEST_INTEGRATOR_STRAT_H
#define PHARE_TEST_INTEGRATOR_STRAT_H




#include <SAMRAI/algs/TimeRefinementLevelStrategy.h>




class TestIntegratorStrat : public SAMRAI::algs::TimeRefinementLevelStrategy
{
public:
    // -----------------------------------------------------------------------------------------------
    //
    //                          SAMRAI TimeRefinementLevelStrategy interface
    //
    // -----------------------------------------------------------------------------------------------



    virtual void initializeLevelIntegrator(
        const std::shared_ptr<SAMRAI::mesh::GriddingAlgorithmStrategy>& griddingAlg) override
    {
    }

    virtual double getLevelDt(const std::shared_ptr<SAMRAI::hier::PatchLevel>& level,
                              const double dtTime, const bool initialTime) override
    {
        return dtTime;
    }


    virtual double getMaxFinerLevelDt(const int finerLevelNumber, const double coarseDt,
                                      const SAMRAI::hier::IntVector& ratio) override
    {
        return std::pow(coarseDt, ratio.max());
    }




    virtual double advanceLevel(const std::shared_ptr<SAMRAI::hier::PatchLevel>& level,
                                const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                                const double currentTime, const double newTime,
                                const bool firstStep, const bool lastStep,
                                const bool regridAdvance = false) override
    {
        return newTime;
    }




    virtual void standardLevelSynchronization(
        const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy, const int coarsestLevel,
        const int finestLevel, const double syncTime, const std::vector<double>& oldTimes) override
    {
        // TODO use messengers to sync with coarser
    }

    virtual void
    synchronizeNewLevels(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                         const int coarsestLevel, const int finestLevel, const double syncTime,
                         const bool initialTime) override
    {
    }


    virtual void resetTimeDependentData(const std::shared_ptr<SAMRAI::hier::PatchLevel>& level,
                                        const double newTime, const bool canBeRefined) override
    {
    }

    virtual void
    resetDataToPreadvanceState(const std::shared_ptr<SAMRAI::hier::PatchLevel>& level) override
    {
    }

    virtual bool usingRefinedTimestepping() const override
    {
        // refined time stepping not allowed
        // so that regridding occurs right away at the first advance
        return false;
    }
};




#endif
