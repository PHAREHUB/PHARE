#ifndef CONST_HYBRID_TAGGER_STRATEGY_H
#define CONST_HYBRID_TAGGER_STRATEGY_H

#include "hybrid_tagger_strategy.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include <cstddef>


namespace PHARE::amr
{
template<typename HybridModel>
class ConstHybridTaggerStrategy : public HybridTaggerStrategy<HybridModel>
{
    using gridlayout_type           = typename HybridModel::gridlayout_type;
    static auto constexpr dimension = HybridModel::dimension;

public:
    void tag(HybridModel& model, gridlayout_type const& layout, int* tags) const override;
};

template<typename HybridModel>
void ConstHybridTaggerStrategy<HybridModel>::tag(HybridModel& /*model*/,
                                                 gridlayout_type const& layout, int* tags) const
{
    // SAMRAI tags int* buffer is FORTRAN ordering so we set false to the view
    bool constexpr c_ordering = false;
    auto tagsv = core::NdArrayView<dimension, int, int*, c_ordering>(tags, layout.nbrCells());

    std::fill_n(tagsv.data(), tagsv.size(), 1);
}
} // namespace PHARE::amr

#endif // CONST_HYBRID_TAGGER_STRATEGY_H
