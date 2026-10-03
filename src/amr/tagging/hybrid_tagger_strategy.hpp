#ifndef HYBRID_TAGGER_STRATEGY_HPP
#define HYBRID_TAGGER_STRATEGY_HPP

namespace PHARE::amr
{

template<typename HybridModel>
class HybridTaggerStrategy
{
    using gridlayout_type = typename HybridModel::gridlayout_type;

public:
    virtual void tag(HybridModel& model, gridlayout_type const& layout, int* tags) const = 0;
    virtual ~HybridTaggerStrategy()                                                      = 0;
};

template<typename HybridModel>
HybridTaggerStrategy<HybridModel>::~HybridTaggerStrategy()
{
}
} // namespace PHARE::amr

#endif // HYBRID_TAGGER_STRATEGY_HPP
