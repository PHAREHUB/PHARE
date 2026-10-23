#ifndef PHARE_CORE_GRID_DETAIL_RAJA_H
#define PHARE_CORE_GRID_DETAIL_RAJA_H

#include "core/def.hpp"

namespace PHARE::core::gpu::impl
{

struct GridLayout
{
    template<typename Fn, typename T, std::size_t S, typename... Args>
    static void evalOnBox(Fn& fn, Box<T, S> const& box, Args&&... args)
    {
#if !defined(RAJA_ENABLE_CUDA)
        RAJA::forall<RAJA::seq_exec>(RAJA::RangeSegment(0, box.size()),
                                     [=] RAJA_HOST_DEVICE(std::size_t const i) mutable {
                                         auto b = box.begin();
                                         b      = b + i;
                                         fn(*b, args...);
                                     });

#else
        RAJA::resources::Cuda res;
        exec(
            [=] RAJA_DEVICE(int i) mutable {
                auto b = box.begin();
                b      = b + i;
                fn(*b, args...);
            },
            res, box.size());

#endif // defined(RAJA_ENABLE_CUDA)
    }
};


} // namespace PHARE::core::gpu::impl

#endif /* PHARE_CORE_GRID_DETAIL_RAJA_H */
