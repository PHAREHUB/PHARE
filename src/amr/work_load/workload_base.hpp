
#ifndef PHARE_WORKLOAD_BASE_HPP
#define PHARE_WORKLOAD_BASE_HPP


#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/pdat/CellVariable.h>
#include <SAMRAI/hier/VariableContext.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/pdat/CellData.h>
#include <bits/types/time_t.h>
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

    ~WorkLoadEstimatorBase(){}; // TODO need to ke kept ? we have the default dtor in IWLE

    int getID() { return id_; };

    void allocate_(SAMRAI::hier::Patch&, double const);

protected:
    std::string workLoadName_;
    SAMRAI::tbox::Dimension dimension_;
    std::shared_ptr<SAMRAI::pdat::CellVariable<double>> workLoadVariable_;
    SAMRAI::hier::VariableDatabase* variableDatabase_;
    std::shared_ptr<SAMRAI::hier::VariableContext> context_;
    int const id_;
};


template<typename PHARE_T>
void WorkLoadEstimatorBase<PHARE_T>::allocate_(SAMRAI::hier::Patch& patch,
                                               double const allocateTime)
{
    patch.allocatePatchData(id_, allocateTime);
};


} // namespace PHARE::amr

#endif
