
#ifndef PHARE_WORKLOAD_HPP
#define PHARE_WORKLOAD_HPP


#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/pdat/CellVariable.h>
#include <SAMRAI/hier/VariableContext.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/pdat/CellData.h>



namespace PHARE::amr
{
template<typename PHARE_T>
class IWorkLoadEstimator
{
    static constexpr std::string workLoadName{"workLoad_" + PHARE_T::dimension + "_"
                                              + PHARE_T::interp_order + "_"
                                              + PHARE_T::nbRefinedPart};

public:
    IWorkLoadEstimator()
        : variableDatabase_{SAMRAI::hier::VariableDatabase::getDatabase()}
        , dimension_{SAMRAI::tbox::Dimension{PHARE_T::dim}}
        , workLoadVariable_{std::make_shared<SAMRAI::pdat::CellVariable<double>>(
              dimension_, workLoadName)} // TODO make a name with dim, iterporder, nbr_refine
        , context_{variableDatabase_->getContext("default")}
        , id_{variableDatabase_->registerVariableAndContext(
              workLoadVariable_, context_, SAMRAI::hier::IntVector::getZero(dimension_))}
    {
        // TODO have the correct initialization in the ctor above
    }

    virtual void estimate(SAMRAI::hier::PatchLevel, double*,
                          PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types> const&)
        = 0;
    virtual void set_strategy(std::string) = 0;
    int getID() { return id_; };
    virtual std::string name() const = 0;

private:
    SAMRAI::hier::VariableDatabase* variableDatabase_;
    SAMRAI::tbox::Dimension dimension_;
    // std::shared_ptr<SAMRAI::pdat::CellData<double>> workLoad_; now this var is locally defined in
    // advanceLevel
    std::shared_ptr<SAMRAI::pdat::CellVariable<double>> workLoadVariable_;
    // std::string contextName_{"default"}; directly called in the init list of the ctor
    std::shared_ptr<SAMRAI::hier::VariableContext> context_;
    // SAMRAI::hier::IntVector gw0 = SAMRAI::hier::IntVector::getZero(dimension_); also directly
    // called i, the init list of the ctor

protected:
    int const id_;
    // int const id_{12 };   // TODO to be initialized in the ctor
};

} // namespace PHARE::amr

#endif
