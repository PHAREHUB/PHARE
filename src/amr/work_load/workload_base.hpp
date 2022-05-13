
#ifndef PHARE_WORKLOAD_BASE_HPP
#define PHARE_WORKLOAD_BASE_HPP


#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/pdat/CellVariable.h>
#include <SAMRAI/hier/VariableContext.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/pdat/CellData.h>
#include <string>

#include "amr/types/amr_types.hpp"
#include "amr/solvers/solver.hpp"


namespace PHARE::amr
{
template<typename PHARE_T>
class WorkLoadEstimatorBase
{
public:
    WorkLoadEstimatorBase(std::string physicalTypeOfWorkLoad)
        : workLoadName_{physicalTypeOfWorkLoad + "_" + std::to_string(PHARE_T::dimension) + "_"
                        + std::to_string(PHARE_T::interp_order) + "_"
                        + std::to_string(PHARE_T::nbRefinedPart)}
        , dimension_{SAMRAI::tbox::Dimension{PHARE_T::dimension}}
        , workLoadVariable_{std::make_shared<SAMRAI::pdat::CellVariable<double>>(dimension_,
                                                                                 workLoadName_)}
        , variableDatabase_{SAMRAI::hier::VariableDatabase::getDatabase()}
        , context_{variableDatabase_->getContext("default")}
        , id_{variableDatabase_->registerVariableAndContext(
              workLoadVariable_, context_, SAMRAI::hier::IntVector::getZero(dimension_))}
    {
    }

    ~WorkLoadEstimatorBase()
    {
        assert(variableDatabase_->checkVariableExists(workLoadName_));
        std::cout << variableDatabase_->checkVariableExists(workLoadName_) << std::endl;
        std::cout << "zob" << std::endl;
        variableDatabase_->removeVariable(workLoadName_);
    };

    virtual void estimate(SAMRAI::hier::PatchLevel, double*,
                          PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types> const&)
        = 0;
    virtual void set_strategy(std::string) = 0;
    int getID() { return id_; };
    virtual std::string name() const = 0;

protected:
    std::string workLoadName_;
    SAMRAI::tbox::Dimension dimension_;
    std::shared_ptr<SAMRAI::pdat::CellVariable<double>> workLoadVariable_;
    SAMRAI::hier::VariableDatabase* variableDatabase_;
    std::shared_ptr<SAMRAI::hier::VariableContext> context_;
    int const id_;
};

} // namespace PHARE::amr

#endif
