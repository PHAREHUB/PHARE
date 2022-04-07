
#ifndef PHARE_HYBRID_WORKLOAD_STRATEGY_HPP
#define PHARE_HYBRID_WORKLOAD_STRATEGY_HPP



namespace PHARE::core
{
template<typename HybridModel>
class HybridWorkLoadEstimatorStrategy
{
protected:
    using amr_t = PHARE::amr::SAMRAI_Types;

public :
    virtual void estimate(double*, HybridModel const&) = 0;
    virtual std::string name() = 0;

};
} // namespace PHARE::core

#endif
