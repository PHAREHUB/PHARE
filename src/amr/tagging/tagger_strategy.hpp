#ifndef TAGGER_STRATEGY_HPP
#define TAGGER_STRATEGY_HPP

namespace PHARE::amr
{

template<typename Model>
class TaggerStrategy
{
    using gridlayout_type = typename Model::gridlayout_type;

public:
    virtual void tag(Model& model, gridlayout_type const& layout, int* tags) const = 0;
    virtual ~TaggerStrategy()                                                      = 0;
};

template<typename Model>
TaggerStrategy<Model>::~TaggerStrategy()
{
}
} // namespace PHARE::amr

#endif // HYBRID_TAGGER_STRATEGY_HPP
