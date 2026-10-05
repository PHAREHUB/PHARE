#ifndef PHARE_TEST_INTEGRATOR_STRAT_HPP
#define PHARE_TEST_INTEGRATOR_STRAT_HPP




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
        std::shared_ptr<SAMRAI::mesh::GriddingAlgorithmStrategy> const& /*griddingAlg*/) override
    {
    }

    virtual double getLevelDt(std::shared_ptr<SAMRAI::hier::PatchLevel> const& /*level*/,
                              double const dtTime, bool const /*initialTime*/) override
    {
        return dtTime;
    }


    virtual double getMaxFinerLevelDt(int const /*finerLevelNumber*/, double const coarseDt,
                                      SAMRAI::hier::IntVector const& ratio) override
    {
        return std::pow(coarseDt, ratio.max());
    }




    virtual double advanceLevel(std::shared_ptr<SAMRAI::hier::PatchLevel> const& /*level*/,
                                std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& /*hierarchy*/,
                                double const /*currentTime*/, double const newTime,
                                bool const /*firstStep*/, bool const /*lastStep*/,
                                [[maybe_unused]] bool const regridAdvance = false) override
    {
        return newTime;
    }




    virtual void
    standardLevelSynchronization(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& /*hierarchy*/,
                                 int const /*coarsestLevel*/, int const /*finestLevel*/,
                                 double const /*syncTime*/,
                                 std::vector<double> const& /*oldTimes*/) override
    {
        // TODO use messengers to sync with coarser
    }

    virtual void
    synchronizeNewLevels(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& /*hierarchy*/,
                         int const /*coarsestLevel*/, int const /*finestLevel*/,
                         double const /*syncTime*/, bool const /*initialTime*/) override
    {
    }


    virtual void resetTimeDependentData(std::shared_ptr<SAMRAI::hier::PatchLevel> const& /*level*/,
                                        double const /*newTime*/,
                                        bool const /*canBeRefined*/) override
    {
    }

    virtual void
    resetDataToPreadvanceState(std::shared_ptr<SAMRAI::hier::PatchLevel> const& /*level*/) override
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
