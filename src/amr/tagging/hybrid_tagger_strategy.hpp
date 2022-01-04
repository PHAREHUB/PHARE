#ifndef HYBRID_TAGGER_STRATEGY_H
#define HYBRID_TAGGER_STRATEGY_H

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
}

#endif // HYBRID_TAGGER_STRATEGY_H
