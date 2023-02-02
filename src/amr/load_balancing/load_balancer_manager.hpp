#ifndef PHARE_LOAD_BALANCER_MANAGER_HPP
#define PHARE_LOAD_BALANCER_MANAGER_HPP

#include <memory>
#include <vector>

#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/pdat/CellVariable.h>
#include "phare_core.hpp"
#include "load_balancer_estimator.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/solvers/solver.hpp"


namespace PHARE::amr
{
template<std::size_t dim>
class LoadBalancerManager
{
public:
    LoadBalancerManager() // TODO doit prendre en arg le max level number
        : dim_{SAMRAI::tbox::Dimension{dim}}
        , loadBalancerVar_{std::make_shared<SAMRAI::pdat::CellVariable<double>>(
              dim_, "LoadBalancerVariable")}
        , variableDatabase_{SAMRAI::hier::VariableDatabase::getDatabase()}
        , context_{variableDatabase_->getContext("default")}
        , id_{variableDatabase_->registerVariableAndContext(
              loadBalancerVar_, context_, SAMRAI::hier::IntVector::getZero(dim_))} {};

    ~LoadBalancerManager() { variableDatabase_->removeVariable("LoadBalancerVariable"); };

    void addLoadBalancerEstimator(
        std::unique_ptr<amr::LoadBalancerEstimator<dim>>
            lbe); // TODO doit aussi prendre le level sur lequel on l'enregistre

    int numOfEstimators() const;

    // std::shared_ptr<amr::LoadBalancerEstimator> getLoadBalancerEstimator(int estimator_index);

    void allocate(SAMRAI::hier::Patch& patch, double const allocateTime);

    void estimate(SAMRAI::hier::PatchLevel& level,
                  PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types> const& model);

private:
    SAMRAI::tbox::Dimension dim_;
    std::shared_ptr<SAMRAI::pdat::CellVariable<double>> loadBalancerVar_;
    SAMRAI::hier::VariableDatabase* variableDatabase_;
    std::shared_ptr<SAMRAI::hier::VariableContext> context_;
    int const id_;
    std::vector<std::shared_ptr<amr::LoadBalancerEstimator<dim>>>
        loadBalancerEstimators_; // TODO ce vector doit avoir maxLevelNumber cases, qui contiendront
                                 // toutes des LBEHybrid
};



template<std::size_t dim>
inline void LoadBalancerManager<dim>::addLoadBalancerEstimator(
    std::unique_ptr<amr::LoadBalancerEstimator<dim>> lbe)
{
    if (core::notIn(lbe, loadBalancerEstimators_))
    {
        loadBalancerEstimators_.push_back(std::move(lbe));
    }
    else
    {
        throw std::runtime_error("load balancer estimator " + lbe->name() + " already registered");
    }
}



template<std::size_t dim>
inline int LoadBalancerManager<dim>::numOfEstimators() const
{
    return loadBalancerEstimators_.size() - 1;
}



// inline std::shared_ptr<amr::LoadBalancerEstimator>
// LoadBalancerManager::getLoadBalancerEstimator(int estimator_index)
// {
//     if (loadBalancerEstimators_[estimator_index] == nullptr)
//     {
//         throw std::runtime_error("no load balancer assigned to this index");
//     }
//
//    return loadBalancerEstimators_[estimator_index];
//}



template<std::size_t dim>
inline void LoadBalancerManager<dim>::allocate(SAMRAI::hier::Patch& patch,
                                               double const allocateTime)
{
    patch.allocatePatchData(id_, allocateTime);
}



template<std::size_t dim>
inline void LoadBalancerManager<dim>::estimate(
    SAMRAI::hier::PatchLevel& level,
    PHARE::solver::IPhysicalModel<PHARE::amr::SAMRAI_Types> const& model)
{
}

} // namespace PHARE::amr

#endif
