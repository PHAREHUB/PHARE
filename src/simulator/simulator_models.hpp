#ifndef PHARE_SIMULATOR_SIMULATOR_MODELS_HPP
#define PHARE_SIMULATOR_SIMULATOR_MODELS_HPP

#include "simulator/simulator_decl.hpp"

#include <memory>
#include <sstream>

namespace PHARE
{

template<auto opts>
void Simulator<opts>::hybrid_model_init(initializer::PHAREDict const& dict)
{
    if constexpr (has_hybrid_v<opts>)
    {
        using HybridModel = PHARETypes::Hybrid::Model_t;
        using ResMan_t    = HybridModel::resources_manager_type;

        hyb_.resman_ = std::make_shared<ResMan_t>();
        hyb_.model_  = std::make_shared<HybridModel>(dict["simulation"], hyb_.resman_);
        hyb_.model_->registerResources();
    }
}


template<auto opts>
void Simulator<opts>::mhd_model_init(initializer::PHAREDict const& dict)
{
    if constexpr (has_mhd_v<opts>)
    {
        using MHDModel = PHARETypes::MHD::Model_t;
        using ResMan_t = MHDModel::resources_manager_type;

        mhd_.resman_ = std::make_shared<ResMan_t>();
        mhd_.model_  = std::make_shared<MHDModel>(dict["simulation"], mhd_.resman_);
        mhd_.model_->registerResources();
    }
}


template<auto opts>
std::string Simulator<opts>::models_to_str() const
{
    std::stringstream ss;
    if constexpr (has_hybrid_v<opts>)
    {
        if (hyb_.model_)
            ss << core::to_str(hyb_.model_->state);
    }
    if constexpr (has_mhd_v<opts>)
    {
        if (mhd_.model_)
            ss << core::to_str(mhd_.model_->state);
    }
    return ss.str();
}


template<auto opts>
std::string Simulator<opts>::summary() const
{
    std::stringstream ss;

    if constexpr (has_hybrid_v<opts>)
    {
        if (hyb_.model_)
            ss << hyb_.model_->summarize(*hierarchy_);
    }

    return ss.str();
}

} // namespace PHARE

#endif /*PHARE_SIMULATOR_SIMULATOR_MODELS_HPP*/
