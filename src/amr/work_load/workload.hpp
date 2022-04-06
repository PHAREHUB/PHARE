
#ifndef PHARE_WORKLOAD_HPP
#define PHARE_WORKLOAD_HPP


#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/pdat/CellVariable.h>
#include <SAMRAI/hier/VariableContext.h>
#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/pdat/CellData.h>

#include "amr/types/amr_types.hpp"


namespace PHARE::amr
{
class IWorkLoadEstimator
{
protected:
    using amr_t = PHARE::amr::SAMRAI_Types;

public:
    IWorkLoadEstimator() { }

    virtual void estimate(SAMRAI::hier::PatchLevel lev, double*, PHARE::solver::IPhysicalModel<amr_t> const& model) = 0;
    virtual void set_strategy(std::string) = 0;
    int getID() { return id_; };

private:
    SAMRAI::tbox::Dimension dim_{2};
    std::shared_ptr<SAMRAI::pdat::CellData<double>> workLoad_;
    std::shared_ptr<SAMRAI::pdat::CellVariable<double>> workLoadVariable_;
    std::string contextName_{"default"};
    std::shared_ptr<SAMRAI::hier::VariableContext> context_;
    SAMRAI::hier::IntVector gw0 = SAMRAI::hier::IntVector::getZero(dim_);
    int const id_{12};
};

} // namespace PHARE::amr

#endif
