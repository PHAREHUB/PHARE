#ifndef PHARE_LOAD_BALANCER_MANAGER_HPP
#define PHARE_LOAD_BALANCER_MANAGER_HPP

#include <memory>
#include <vector>
#include <algorithm>

#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/pdat/CellVariable.h>
// #include "phare_core.hpp"
#include "initializer/data_provider.hpp"
#include "load_balancer_estimator.hpp"
// #include "amr/resources_manager/amr_utils.hpp"
// #include "amr/solvers/solver.hpp"



namespace PHARE::amr
{
template<std::size_t dim>
class LoadBalancerManager
{
public:
    // LoadBalancerManager(int const maxLevelNumber)
    LoadBalancerManager(PHARE::initializer::PHAREDict const& dict)
        : dim_{SAMRAI::tbox::Dimension{dim}}
        , loadBalancerVar_{std::make_shared<SAMRAI::pdat::CellVariable<double>>(
              dim_, "LoadBalancerVariable")}
        , variableDatabase_{SAMRAI::hier::VariableDatabase::getDatabase()}
        , context_{variableDatabase_->getContext("default")}
        , id_{variableDatabase_->registerVariableAndContext(loadBalancerVar_, context_,
                                                            SAMRAI::hier::IntVector::getZero(dim_))}
        , maxLevelNumber_{dict["simulation"]["AMR"]["max_nbr_levels"].template to<int>()}
        // , loadBalancerEstimators_(maxLevelNumber, nullptr){};
        , loadBalancerEstimators_(maxLevelNumber_, nullptr){};

    ~LoadBalancerManager() { variableDatabase_->removeVariable("LoadBalancerVariable"); };

    int getId() const;

    void addLoadBalancerEstimator(int const iLevel_min, int const iLevel_max,
                                  std::shared_ptr<amr::LoadBalancerEstimator> lbe);

    void allocate(SAMRAI::hier::Patch& patch, double const allocateTime);

    void estimate(SAMRAI::hier::PatchLevel& level,
                  PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types>& model);


private:
    SAMRAI::tbox::Dimension dim_;
    std::shared_ptr<SAMRAI::pdat::CellVariable<double>> loadBalancerVar_;
    SAMRAI::hier::VariableDatabase* variableDatabase_;
    std::shared_ptr<SAMRAI::hier::VariableContext> context_;
    int const id_;
    int const maxLevelNumber_;
    std::vector<std::shared_ptr<amr::LoadBalancerEstimator>> loadBalancerEstimators_;
};




template<std::size_t dim>
inline int LoadBalancerManager<dim>::getId() const
{
    return id_;
}




template<std::size_t dim>
inline void
LoadBalancerManager<dim>::addLoadBalancerEstimator(int const iLevel_min, int const iLevel_max,
                                                   std::shared_ptr<amr::LoadBalancerEstimator> lbe)
{
    for (auto ilevel = iLevel_min; ilevel <= iLevel_max; ilevel++)
    {
        loadBalancerEstimators_[ilevel] = lbe;
    }
}




template<std::size_t dim>
inline void LoadBalancerManager<dim>::allocate(SAMRAI::hier::Patch& patch,
                                               double const allocateTime)
{
    patch.allocatePatchData(id_, allocateTime);
}



template<std::size_t dim>
inline void
LoadBalancerManager<dim>::estimate(SAMRAI::hier::PatchLevel& level,
                                   PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types>& model)
{
    auto iLevel = level.getLevelNumber();
    auto lbe    = loadBalancerEstimators_[iLevel];

    lbe->estimate(level, model);
}

} // namespace PHARE::amr

#endif
