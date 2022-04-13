
#ifndef PHARE_WORKLOAD_HPP
#define PHARE_WORKLOAD_HPP


#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/pdat/CellVariable.h>
#include <SAMRAI/hier/VariableContext.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/pdat/CellData.h>



namespace PHARE::amr
{
template<std::size_t dim>
class IWorkLoadEstimator
{
public:
    IWorkLoadEstimator()
        : variableDatabase_{SAMRAI::hier::VariableDatabase::getDatabase()}
        , dimension_{SAMRAI::tbox::Dimension{dim}}
        , workLoadVariable_{std::make_shared<SAMRAI::pdat::CellVariable<double>>(dimension_,
                                                                                 "workLoad")}
        , context_{variableDatabase_->getContext("default")}
    // , id_{variableDatabase_->registerVariableAndContext(workLoadVariable_, context_, gw0)}
    {
        // TODO have the correct initialization in the ctor above
    }

    virtual void estimate(SAMRAI::hier::PatchLevel, double*,
                          PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types> const&)
        = 0;
    virtual void set_strategy(std::string) = 0;
    int getID() { return id_; };
    virtual std::string name() const = 0;

protected:
    int const id_{12}; // TODO to be initialized in the ctor

private:
    SAMRAI::hier::VariableDatabase* variableDatabase_;
    SAMRAI::tbox::Dimension dimension_;
    // std::shared_ptr<SAMRAI::pdat::CellData<double>> workLoad_;
    std::shared_ptr<SAMRAI::pdat::CellVariable<double>> workLoadVariable_;
    // std::string contextName_{"default"};
    std::shared_ptr<SAMRAI::hier::VariableContext> context_;
    SAMRAI::hier::IntVector gw0 = SAMRAI::hier::IntVector::getZero(dimension_);
};

} // namespace PHARE::amr

#endif
