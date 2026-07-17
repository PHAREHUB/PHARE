#ifndef PHARE_CORE_GRID_DETAIL_MKN_GPU_HPP
#define PHARE_CORE_GRID_DETAIL_MKN_GPU_HPP

#include "core/def.hpp"

namespace PHARE::core::gpu::impl
{

struct GridLayout
{
    template<typename Fn, typename T, std::size_t S, typename... Args>
    static void evalOnBox(Fn& fn, Box<T, S> const& box, Args&... args)
    {
        PHARE_LOG_LINE_SS(box.size());
        mkn::gpu::GDLauncher{box.size()}([=] _PHARE_ALL_FN_() mutable {
            auto b = box.begin();
            b      = b + mkn::gpu::idx();
            fn(*b, args...);
        });
    }
};


} // namespace PHARE::core::gpu::impl


#endif /* PHARE_CORE_GRID_DETAIL_MKN_GPU_HPP */
