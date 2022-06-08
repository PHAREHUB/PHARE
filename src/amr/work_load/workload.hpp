
#ifndef PHARE_WORKLOAD_HPP
#define PHARE_WORKLOAD_HPP


#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/pdat/CellVariable.h>
#include <SAMRAI/hier/VariableContext.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/pdat/CellData.h>

#include "amr/types/amr_types.hpp"
#include "amr/solvers/solver.hpp"



namespace PHARE::amr
{
class IWorkLoadEstimator
{
public:
    virtual void estimate(SAMRAI::hier::PatchLevel,
                          PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types> const&)
        = 0;
    virtual void set_strategy(std::string)    = 0;
    virtual std::string name() const          = 0;
    virtual void terminateWorkLoadEstimator() = 0;
    virtual ~IWorkLoadEstimator()             = default;
};

} // namespace PHARE::amr

#endif
