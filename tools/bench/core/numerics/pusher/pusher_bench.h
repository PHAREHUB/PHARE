#ifndef PHARE_CORE_PUSHER_BENCH_H
#define PHARE_CORE_PUSHER_BENCH_H

#include "phare_core.h"
#include "core/numerics/pusher/boris.h"
#include "core/numerics/ion_updater/ion_updater.h"

namespace PHARE::core::bench
{
template<std::size_t dim>
using Field = PHARE::core::Field<PHARE::core::NdArrayVector<dim>,
                                 typename PHARE::core::HybridQuantity::Scalar>;
template<std::size_t dim>
using VecField
    = PHARE::core::VecField<PHARE::core::NdArrayVector<dim>, typename PHARE::core::HybridQuantity>;


template<std::size_t dim>
PHARE::core::Particle<dim> particle(int icell = 5)
{
    return {//
            /*.weight = */ 0,
            /*.charge = */ 1,
            /*.iCell  = */ PHARE::core::ConstArray<int, dim>(icell),
            /*.delta  = */ PHARE::core::ConstArray<float, dim>(.5),
            /*.v      = */ {{0, 10., 0}}};
}

template<typename GridLayout, typename Quantity, std::size_t dim = GridLayout::dimension>
Field<dim> field(std::string key, Quantity type, GridLayout const& layout)
{
    Field<dim> feeld{key, type, layout.allocSize(type)};
    std::fill(feeld.begin(), feeld.end(), 1);
    return feeld;
}


template<typename GridLayout>
auto E(GridLayout const& layout)
{
    return std::make_tuple(field("Ex", PHARE::core::HybridQuantity::Scalar::Ex, layout),
                           field("Ey", PHARE::core::HybridQuantity::Scalar::Ey, layout),
                           field("Ez", PHARE::core::HybridQuantity::Scalar::Ez, layout));
}

template<typename GridLayout>
auto B(GridLayout const& layout)
{
    return std::make_tuple(field("Bx", PHARE::core::HybridQuantity::Scalar::Bx, layout),
                           field("By", PHARE::core::HybridQuantity::Scalar::By, layout),
                           field("Bz", PHARE::core::HybridQuantity::Scalar::Bz, layout));
}

template<typename GridLayout, std::size_t dim = GridLayout::dimension>
auto EM(GridLayout const& layout)
{
    return std::make_tuple(E(layout), B(layout));
}

template<typename GridLayout, typename VecFieldT>
class Electromag : public PHARE::core::Electromag<VecFieldT>
{
public:
    using Super = PHARE::core::Electromag<VecFieldT>;

    Electromag(GridLayout const& layout)
        : Super{"EM"}
        , emFields{EM(layout)}
    {
        auto& [E, B]       = emFields;
        auto& [ex, ey, ez] = E;
        auto& [bx, by, bz] = B;

        Super::B.setBuffer("EM_B_x", &bx);
        Super::B.setBuffer("EM_B_y", &by);
        Super::B.setBuffer("EM_B_z", &bz);
        Super::E.setBuffer("EM_E_x", &ex);
        Super::E.setBuffer("EM_E_y", &ey);
        Super::E.setBuffer("EM_E_z", &ez);
    }

private:
    decltype(EM(*static_cast<GridLayout*>(0))) emFields;
};

} // namespace PHARE::core::bench

#endif /*PHARE_CORE_PUSHER_BENCH_H*/
