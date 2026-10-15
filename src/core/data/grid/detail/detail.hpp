#ifndef PHARE_CORE_GRID_DETAIL_DETAIL_HPP
#define PHARE_CORE_GRID_DETAIL_DETAIL_HPP

#if PHARE_HAVE_RAJA
#include "core/data/grid/detail/raja.hpp"
#endif

#if PHARE_HAVE_MKN_GPU
#include "core/data/grid/detail/mkn_gpu.hpp"
#endif

#if PHARE_HAVE_KOKKOS
#include "core/data/grid/detail/kokkos.hpp"
#endif

#include "core/utilities/box/box.hpp"

#include <cstddef>

namespace PHARE::core::gpu
{
struct GridLayout
{
    template<typename Fn, typename T, std::size_t S, typename... Args>
    static void evalOnBox(Fn& fn, Box<T, S> const& box, Args&... args)
    {
        PHARE_LOG_LINE_SS(box.size());


        PHARE_WITH_RAJA({ //
            impl::GridLayout::evalOnBox(fn, box, args...);
        })

        PHARE_WITH_KOKKOS({ //
            impl::GridLayout::evalOnBox(fn, box, args...);
        })


        PHARE_WITH_MKN_GPU({ //
            impl::GridLayout::evalOnBox(fn, box, args...);
        })
    }
};


} // namespace PHARE::core::gpu



#endif /* PHARE_CORE_GRID_DETAIL_DETAIL_HPP */
