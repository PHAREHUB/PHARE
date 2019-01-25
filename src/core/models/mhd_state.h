#ifndef PHARE_MHD_STATE_H
#define PHARE_MHD_STATE_H

#include "hybrid/hybrid_quantities.h"
#include "models/physical_state.h"

namespace PHARE
{
namespace core
{
    using MHDQuantity = HybridQuantity;

    class MHDStateInitializer : public PhysicalStateInitializer
    {
    };


    template<typename VecFieldT>
    class MHDState : public IPhysicalState
    {
    public:
        /*virtual void allocate(ResourcesManager const& manager, SAMRAI::hier::Patch& patch)
        override
        {
            manager.allocate(B, patch);
            manager.allocate(V, patch);
        }*/

        VecFieldT B{"B", MHDQuantity::Vector::B};
        VecFieldT V{"V", MHDQuantity::Vector::V};
    };
} // namespace core
} // namespace PHARE



#endif // PHARE_MHD_STATE_H
