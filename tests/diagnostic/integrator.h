
#ifndef PHARE_TEST_DIAGNOSTIC_INTEGRATOR
#define PHARE_TEST_DIAGNOSTIC_INTEGRATOR

class IntegratorStrat : public SAMRAI::algs::TimeRefinementLevelStrategy
{
public:
    // -----------------------------------------------------------------------------------------------
    //
    //                          SAMRAI TimeRefinementLevelStrategy interface
    //
    // -----------------------------------------------------------------------------------------------

    virtual void initializeLevelIntegrator([
        [maybe_unused]] std::shared_ptr<SAMRAI::mesh::GriddingAlgorithmStrategy> const& griddingAlg)
        override
    {
    }

    virtual double
    getLevelDt([[maybe_unused]] const std::shared_ptr<SAMRAI::hier::PatchLevel>& level,
               double const dtTime, [[maybe_unused]] bool const initialTime) override
    {
        return dtTime;
    }

    virtual double getMaxFinerLevelDt([[maybe_unused]] int const finerLevelNumber,
                                      double const coarseDt,
                                      const SAMRAI::hier::IntVector& ratio) override
    {
        return std::pow(coarseDt, ratio.max());
    }

    virtual double
    advanceLevel([[maybe_unused]] std::shared_ptr<SAMRAI::hier::PatchLevel> const& level,
                 [[maybe_unused]] std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                 [[maybe_unused]] double const currentTime, double const newTime,
                 [[maybe_unused]] bool const firstStep, [[maybe_unused]] bool const lastStep,
                 [[maybe_unused]] bool const regridAdvance = false) override
    {
        return newTime;
    }

    virtual void standardLevelSynchronization(
        [[maybe_unused]] std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
        [[maybe_unused]] int const coarsestLevel, [[maybe_unused]] int const finestLevel,
        [[maybe_unused]] double const syncTime,
        [[maybe_unused]] std::vector<double> const& oldTimes) override
    {
        // TODO use messengers to sync with coarser
    }

    virtual void synchronizeNewLevels(
        [[maybe_unused]] std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
        [[maybe_unused]] int const coarsestLevel, [[maybe_unused]] int const finestLevel,
        [[maybe_unused]] double const syncTime, [[maybe_unused]] bool const initialTime) override
    {
    }


    virtual void
    resetTimeDependentData([[maybe_unused]] std::shared_ptr<SAMRAI::hier::PatchLevel> const& level,
                           [[maybe_unused]] double const newTime,
                           [[maybe_unused]] bool const canBeRefined) override
    {
    }

    virtual void resetDataToPreadvanceState([
        [maybe_unused]] std::shared_ptr<SAMRAI::hier::PatchLevel> const& level) override
    {
    }

    virtual bool usingRefinedTimestepping() const override
    {
        // refined time stepping not allowed
        // so that regridding occurs right away at the first advance
        return false;
    }
};

#endif /*PHARE_TEST_DIAGNOSTIC_INTEGRATOR*/
