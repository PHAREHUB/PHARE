
#ifndef PHARE_HYBRID_WORKLOAD_STRATEGY_HPP
#define PHARE_HYBRID_WORKLOAD_STRATEGY_HPP


#include <string>

#include "amr/types/amr_types.hpp"


namespace PHARE::core
{
template<typename HybridModel>
class HybridWorkLoadEstimatorStrategy
{
protected:
    using amr_t           = PHARE::amr::SAMRAI_Types;
    using gridlayout_type = typename HybridModel::gridlayout_type;

public:
    virtual void estimate(double*, HybridModel const&, gridlayout_type const&) = 0;
    virtual std::string name()                                                 = 0;
};
} // namespace PHARE::core

#endif
