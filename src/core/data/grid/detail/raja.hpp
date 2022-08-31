
#ifndef PHARE_CORE_GRID_DETAIL_RAJA_H
#define PHARE_CORE_GRID_DETAIL_RAJA_H

namespace PHARE::core::raja
{
struct GridLayout
{
    template<typename Fn, typename Index, typename... Args>
    static void evalOnBox(Fn& fn, std::vector<Index>& indexes, Args&... args)
    {
        RAJA::resources::Cuda res;
        auto tuple     = std::make_tuple(args...); // copy for copy :(
        auto gpu_tuple = allocate_copy(res, tuple);
        auto d_indexes = allocate_copy(res, indexes);
        exec(
            [=] RAJA_DEVICE(int i) mutable {
                std::apply(
                    [=](auto&... targs) mutable { //
                        fn(d_indexes[i], targs...);
                    },
                    *gpu_tuple);
            },
            res, indexes.size());
        deallocate(res, gpu_tuple, d_indexes);
    }
};


} // namespace PHARE::core::raja
#endif /* PHARE_CORE_GRID_DETAIL_RAJA_H */
