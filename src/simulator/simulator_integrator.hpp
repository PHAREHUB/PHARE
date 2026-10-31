#ifndef PHARE_SIMULATOR_SIMULATOR_INTEGRATOR_HPP
#define PHARE_SIMULATOR_SIMULATOR_INTEGRATOR_HPP

#include "simulator/simulator_decl.hpp"

#include "amr/tagging/tagger_factory.hpp"
#include "amr/load_balancing/load_balancer_details.hpp"
#include "amr/load_balancing/load_balancer_manager.hpp"
#include "amr/load_balancing/load_balancer_estimator_hybrid.hpp"
#include "amr/load_balancing/load_balancer_estimator_mhd.hpp"

#include <memory>
#include <string>

namespace PHARE
{

// constructing the MultiPhysicsIntegrator emits its vtable, and registering a model creates the
// model's level initializer, both live in this part

template<auto opts>
void Simulator<opts>::multiphysics_init(initializer::PHAREDict const& dict)
{
    multiphysInteg_ = std::make_shared<MultiPhysicsIntegrator>(dict["simulation"], functors_);
}


template<auto opts>
void Simulator<opts>::hybrid_register_model()
{
    if constexpr (has_hybrid_v<opts>)
        multiphysInteg_->registerModel(maxMHDLevel_, maxLevelNumber_ - 1, hyb_.model_);
}


template<auto opts>
void Simulator<opts>::mhd_register_model()
{
    if constexpr (has_mhd_v<opts>)
        multiphysInteg_->registerModel(0, maxMHDLevel_ - 1, mhd_.model_);
}


template<auto opts>
void Simulator<opts>::hybrid_tagger_init(initializer::PHAREDict const& dict)
{
    if constexpr (has_hybrid_v<opts>)
    {
        using HybridModel = PHARETypes::Hybrid::Model_t;

        // hard coded for now, should get some params later from the dict
        if (dict["simulation"]["AMR"]["refinement"].contains("tagging"))
        {
            if (dict["simulation"]["AMR"]["refinement"]["tagging"]["method"]
                    .template to<std::string>()
                != "none")
            {
                auto hybridTagger_ = amr::TaggerFactory<HybridModel>::make(
                    dict["simulation"]["AMR"]["refinement"]["tagging"]);
                multiphysInteg_->registerTagger(maxMHDLevel_, maxLevelNumber_ - 1,
                                                std::move(hybridTagger_));
            }
        }
    }
}


template<auto opts>
void Simulator<opts>::mhd_tagger_init(initializer::PHAREDict const& dict)
{
    if constexpr (has_mhd_v<opts>)
    {
        using MHDModel = PHARETypes::MHD::Model_t;

        if (dict["simulation"]["AMR"]["refinement"].contains("tagging"))
        {
            if (dict["simulation"]["AMR"]["refinement"]["tagging"]["method"]
                    .template to<std::string>()
                != "none")
            {
                auto mhdTagger_ = amr::TaggerFactory<MHDModel>::make(
                    dict["simulation"]["AMR"]["refinement"]["tagging"]);
                multiphysInteg_->registerTagger(0, maxMHDLevel_ - 1, std::move(mhdTagger_));
            }
        }
    }
}


template<auto opts>
void Simulator<opts>::finalize_init(initializer::PHAREDict const& dict)
{
    auto const lb_info = amr::LoadBalancerDetails::FROM(dict["simulation"]["AMR"]["loadbalancing"]);
    auto lbm_          = std::make_unique<amr::LoadBalancerManager<opts.dimension>>(dict);
    auto loadBalancer_db = std::make_shared<SAMRAI::tbox::MemoryDatabase>("LoadBalancerDB");
    loadBalancer_db->putDouble("flexible_load_tolerance", lb_info.tolerance);
    auto loadBalancer = std::make_shared<SAMRAI::mesh::CascadePartitioner>(
        SAMRAI::tbox::Dimension{dimension}, "LoadBalancer", loadBalancer_db);

    // Load balancers break with refinement boxes - only tagging supported
    /*
      P=0000000:Program abort called in file ``/.../SAMRAI/xfer/RefineSchedule.cpp'' at line
      369 P=0000000:ERROR MESSAGE: P=0000000:RefineSchedule:RefineSchedule error: We are not
      currently P=0000000:supporting RefineSchedules with the source level finer
      P=0000000:than the destination level
    */
    if (lb_info.active and !dict["simulation"]["AMR"]["refinement"].contains("tagging"))
        throw std::runtime_error("LoadBalancer is not supported with refinement boxes");

    if (lb_info.active)
    {
        if constexpr (has_hybrid_v<opts>)
        {
            if (find_model("HybridModel"))
            {
                auto lbe_ = std::make_shared<amr::LoadBalancerEstimatorHybrid<PHARETypes>>(
                    lb_info.mode, lbm_->getId());
                lbm_->addLoadBalancerEstimator(maxMHDLevel_, maxLevelNumber_ - 1, std::move(lbe_));
            }
        }
        if constexpr (has_mhd_v<opts>)
        {
            if (find_model("MHDModel"))
            {
                auto lbe_
                    = std::make_shared<amr::LoadBalancerEstimatorMHD<PHARETypes>>(lbm_->getId());
                lbm_->addLoadBalancerEstimator(0, maxMHDLevel_ - 1, std::move(lbe_));
            }
        }
        lbm_->setLoadBalancer(loadBalancer);
    }
    auto const lbm_id = lbm_->getId(); // moved on next line
    multiphysInteg_->setLoadBalancerManager(std::move(lbm_));

    startTime_ = restart_time(dict);
    integrator_
        = std::make_unique<Integrator>(dict, hierarchy_, multiphysInteg_, multiphysInteg_,
                                       loadBalancer, startTime_, finalTime_, lb_info, lbm_id);
    timeStamper = core::TimeStamperFactory::create(dict["simulation"]);
}

} // namespace PHARE

#endif /*PHARE_SIMULATOR_SIMULATOR_INTEGRATOR_HPP*/
