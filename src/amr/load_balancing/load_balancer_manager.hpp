#ifndef PHARE_LOAD_BALANCER_MANAGER_HPP
#define PHARE_LOAD_BALANCER_MANAGER_HPP

#include <memory>
#include <vector>

#include <SAMRAI/hier/PatchLevel.h>
#include "phare_core.hpp"
#include "load_balancer_estimator.hpp"



namespace PHARE::amr
{
template<std::size_t dim>
class LoadBalancerManager
{
public:
    LoadBalancerManager()
        : variableDatabase_{SAMRAI::hier::VariableDatabase::getDatabase()}
        , id_{12000000} {};
    void addLoadBalancerEstimator(std::unique_ptr<amr::LoadBalancerEstimator<dim>> lbe);
    int numOfEstimators() const;
    // std::shared_ptr<amr::LoadBalancerEstimator> getLoadBalancerEstimator(int estimator_index);
    void allocate(SAMRAI::hier::Patch& patch, double const allocateTime);
    // VariableDatabase.h : virtual void removeVariable(const std::string& variable_name); TODO
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

private:
    SAMRAI::hier::VariableDatabase* variableDatabase_;
    int const id_;
    std::vector<std::shared_ptr<amr::LoadBalancerEstimator<dim>>> loadBalancerEstimators_;
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



} // namespace PHARE::amr

#endif
