// IWYU pragma: private, include "core/utilities/kernels/launchers.hpp"

#ifndef PHARE_CORE_UTILITIES_KERNELS_DETAIL_MKN_LAUNCHERS_HPP
#define PHARE_CORE_UTILITIES_KERNELS_DETAIL_MKN_LAUNCHERS_HPP

#include "core/def.hpp"


#if PHARE_HAVE_MKN_GPU

// #include "core/utilities/kernels.hpp"

#include "mkn/gpu.hpp"
#include "mkn/gpu/stream.hpp"

#include <vector>

namespace PHARE::core::gpu
{
// #if !defined(__HIPCC__) || !defined(__CUDACC__) // silence lsp
// inline static dim3 blockIdx;
// inline static dim3 threadIdx;
// #endif

template<typename T>
using Vec_t = std::vector<T, mkn::gpu::ManagedAllocator<T>>;

// mkn::gpu launches are asynchronous, _sync waits for the stream (or device) after each launch
//  Mode: mkn::gpu::Guess (grid/block from a problem size) or mkn::gpu::Fixed (caller provided)
template<typename Mode, bool _sync>
struct DeviceLauncher : public mkn::gpu::DLauncher<Mode>
{
    using Super = mkn::gpu::DLauncher<Mode>;
    using Super::Super;

    template<typename... Args>
    auto& operator()(Args&&... args)
    {
        Super::operator()(std::forward<Args>(args)...);
        if constexpr (_sync)
            this->sync();
        return *this;
    }

    template<typename... Args>
    auto& stream(Args&&... args)
    {
        Super::stream(std::forward<Args>(args)...);
        if constexpr (_sync)
            this->sync();
        return *this;
    }
};

// n threads per cell, 1 block per cell
template<typename Box_t, bool sync = true>
class BoxCellNLauncher : public DeviceLauncher<mkn::gpu::Guess, sync>
{
    using Super                       = DeviceLauncher<mkn::gpu::Guess, sync>;
    static constexpr std::size_t warp = 32; // 32 for nvidia. todo?
public:
    BoxCellNLauncher(Box_t const& box, std::size_t const _n, size_t const& dev = 0)
        : Super{std::size_t{1}, dev}
        , n{_n}
    {
        n_per_cell  = n % warp == 0 ? n : n + (warp - (n % warp));
        this->count = n_per_cell * box.size();
        this->b     = dim3{};
        this->g     = dim3{};
        this->b.x   = n_per_cell;
        this->g.x   = box.size();
    }

    static std::uint32_t block_idx() _PHARE_DEV_FN_ { return blockIdx.x; }
    static std::uint32_t thread_idx() _PHARE_DEV_FN_ { return threadIdx.x; }

private:
    std::size_t n          = 0;
    std::size_t n_per_cell = 0;
};

// 1 block per (tile, cell): grid y = tile, grid x = cell in tile (noop past a tile's cells)
//  block threads stride the cell's particles, so a fixed block size works for any cell size
template<bool sync = false>
class TileCellLauncher : public DeviceLauncher<mkn::gpu::Fixed, sync>
{
    using Super = DeviceLauncher<mkn::gpu::Fixed, sync>;

public:
    TileCellLauncher(std::size_t const n_tiles, std::size_t const max_cells_per_tile,
                     std::size_t const threads_per_cell)
        : Super{dim3{}, dim3{}}
    {
        this->g.x = max_cells_per_tile;
        this->g.y = n_tiles;
        this->b.x = threads_per_cell;
    }

    // largest per-cell local box over AoSPCTS tiles, the grid x size; 0 = nothing to launch
    static std::size_t max_cells(auto const& tiles)
    {
        std::size_t max = 0;
        for (std::size_t tidx = 0; tidx < tiles().size(); ++tidx)
            max = std::max<std::size_t>(max, tiles()[tidx]().local_box().size());
        return max;
    }

    static std::uint32_t tile_idx() _PHARE_DEV_FN_ { return blockIdx.y; }
    static std::uint32_t cell_idx() _PHARE_DEV_FN_ { return blockIdx.x; }
    static std::uint32_t thread_idx() _PHARE_DEV_FN_ { return threadIdx.x; }
};

template<bool sync = false>
class ChunkLauncher : public DeviceLauncher<mkn::gpu::Fixed, sync>
{
    using Super = DeviceLauncher<mkn::gpu::Fixed, sync>;
    // static constexpr std::size_t warp = 64;

public:
    ChunkLauncher(std::size_t const _n, std::size_t const ds = 0, size_t const& /*dev*/ = 0)
        : Super{dim3{}, dim3{}}
        , n{_n}
    {
        this->b   = dim3{};
        this->g   = dim3{};
        this->b.x = n;
        this->b.y = n;
        this->b.z = n;
        this->g.x = 1;
        this->ds  = ds;
    }


private:
    std::size_t n = 0;
};


template<typename Box_t, typename Strat, typename Fn, std::uint16_t impl = 0>
struct BoxStreamDeviceFunction : mkn::gpu::StreamFunction<Strat>
{
    using Super = mkn::gpu::StreamFunction<Strat>;
    using Super::strat;


    BoxStreamDeviceFunction(Strat& strat, Fn&& fn_, std::vector<Box_t> const& boxes_)
        : Super{strat, mkn::gpu::StreamFunctionMode::DEVICE_WAIT}
        , fn{fn_}
        , boxes{boxes_}
    {
    }

    void run(std::uint32_t const i) override
    {
        if constexpr (impl == 0)
        {
            using Launcher       = BoxCellNLauncher<Box_t, false>;
            std::size_t max_size = strat.datas[i]->size() == 0 ? 0 : strat.datas[i]->max_size();
            if (max_size > 0)
                Launcher{boxes[i], max_size}.stream(strat.streams[i],
                                                    [=, fn = fn] __device__() mutable { fn(i); });
        }
        else if constexpr (impl == 1)
        {
            using Launcher = DeviceLauncher<mkn::gpu::Guess, false>;
            Launcher{boxes[i].size()}.stream(strat.streams[i],
                                             [=, fn = fn] __device__() mutable { fn(i); });
        }
        else
            throw std::runtime_error("No impl");
    }


    Fn fn;
    std::vector<Box_t> const& boxes;
};



template<typename Strat, typename Fn, std::uint16_t impl = 0>
struct ChunkFunction : public mkn::gpu::StreamFunction<Strat>
{
    using Super = mkn::gpu::StreamFunction<Strat>;
    using Super::strat;


    ChunkFunction(Strat& strat, Fn&& fn_, std::size_t const n_, std::size_t const ds_)
        : Super{strat, mkn::gpu::StreamFunctionMode::DEVICE_WAIT}
        , fn{fn_}
        , n{n_}
        , ds{ds_}
    {
    }

    void run(std::uint32_t const i) override
    {
        if constexpr (impl == 0)
        {
            using Launcher = ChunkLauncher<false>;
            Launcher{n, ds}.stream(strat.streams[i], [=, fn = fn] __device__() mutable { fn(i); });
        }
        else
            throw std::runtime_error("No impl");
    }

    Fn fn;
    std::size_t n, ds;
};

template<typename Strat, typename Fn, std::uint16_t impl = 0, bool is = true>
struct ChunkIndexFunction : public ChunkFunction<Strat, Fn, impl>
{
    using Super = ChunkFunction<Strat, Fn, impl>;
    using Super::strat;


    ChunkIndexFunction(std::size_t const gs_, std::size_t const gid_, Strat& strat, Fn&& fn,
                       std::size_t const n, std::size_t const ds)
        : Super{strat, std::forward<Fn>(fn), n, ds}
        , gs{gs_}
        , gid{gid_}
    {
    }

    void run(std::uint32_t const i) override
    {
        if constexpr (is)
        {
            if (i % gs == gid)
                Super::run(i);
        }
        else if (i % gs != gid)
            Super::run(i);
    }

    std::size_t const gs;
    std::size_t const gid;
};


template<typename Box_t, typename Strat, typename Fn, std::uint16_t impl = 0, bool is = true>
struct BoxStreamDeviceGroupIndexFunction : BoxStreamDeviceFunction<Box_t, Strat, Fn, impl>
{
    using Super = BoxStreamDeviceFunction<Box_t, Strat, Fn, impl>;
    using Super::strat;


    BoxStreamDeviceGroupIndexFunction(std::size_t const gs_, std::size_t const gid_, Strat& strat,
                                      Fn&& fn, std::vector<Box_t> const& boxes)
        : Super{strat, std::forward<Fn>(fn), boxes}
        , gs{gs_}
        , gid{gid_}
    {
        assert(gs > 0); // groups imply one
    }

    void run(std::uint32_t const i) override
    {
        if constexpr (is)
        {
            if (i % gs == gid)
                Super::run(i);
        }
        else if (i % gs != gid)
            Super::run(i);
    }


    std::size_t const gs;
    std::size_t const gid;
};




template<typename Boxes, typename Vectors>
struct BoxStreamLauncher : public mkn::gpu::StreamLauncher<Vectors>
{
    using Super = mkn::gpu::StreamLauncher<Vectors>;
    using Box_t = Boxes::value_type;

    BoxStreamLauncher(Vectors& vectors, Boxes const& bxes)
        : Super{vectors}
        , boxes{bxes}
    {
    }

    template<std::uint16_t impl = 0, typename Fn>
    auto& async_dev(Fn&& fn)
    {
        using DevFn = BoxStreamDeviceFunction<Box_t, Super, Fn, impl>;
        Super::fns.emplace_back(std::make_shared<DevFn>(**this, std::forward<Fn>(fn), boxes));
        return *this;
    }

    Super& super() { return *this; }
    auto& operator*() { return super(); }

    Boxes const& boxes;
};


template<typename Strat, typename Fn, std::uint16_t impl = 0>
struct TileBlockFunction : public mkn::gpu::StreamFunction<Strat>
{
    using Super = mkn::gpu::StreamFunction<Strat>;
    using Super::strat;


    TileBlockFunction(Strat& strat, Fn&& fn_, std::size_t const ds_)
        : Super{strat, mkn::gpu::StreamFunctionMode::DEVICE_WAIT}
        , fn{fn_}
        , ds{ds_}
    {
    }

    void run(std::uint32_t const i) override
    {
        if (strat.datas[i]->size() == 0)
            return;
        if constexpr (impl == 0)
        {
            using Launcher = ChunkLauncher<false>;
            Launcher launcher{1, ds};
            launcher.b.x = kernel::warp_size();
            launcher.g.x = (*strat.datas[i])().size();
            launcher.stream(strat.streams[i], [=, fn = fn] __device__() mutable { fn(i); });
        }
        else
            throw std::runtime_error("No impl");
    }

    Fn fn;
    std::size_t ds;
};


template<typename Vectors>
struct ThreadedStreamLauncher : public mkn::gpu::ThreadedStreamLauncher<Vectors>
{
    using Super = mkn::gpu::ThreadedStreamLauncher<Vectors>;


    ThreadedStreamLauncher(Vectors& vectors, std::size_t const nthreads = 1)
        : Super{vectors, nthreads}
    {
    }


    template<std::uint16_t impl = 0, typename Fn>
    auto& async_dev_tiles(Fn&& fn, std::size_t const ds = 0)
    {
        if constexpr (impl == 0)
        {
            using DevFn = TileBlockFunction<Super, Fn, impl>;
            Super::fns.emplace_back(std::make_shared<DevFn>(**this, std::forward<Fn>(fn), ds));
        }
        else
        {
            throw std::runtime_error("no impl");
        }
        return *this;
    }

    template<std::uint16_t impl = 0, typename Fn>
    auto& async_dev_chunk(Fn&& fn, std::size_t const n, std::size_t const ds = 0)
    {
        if constexpr (impl == 0)
        {
            using DevFn = ChunkFunction<Super, Fn, impl>;
            Super::fns.emplace_back(std::make_shared<DevFn>(**this, std::forward<Fn>(fn), n, ds));
        }
        else
        {
            throw std::runtime_error("no impl");
        }
        return *this;
    }

    template<std::uint16_t impl = 0, typename Fn>
    auto& async_dev_chunk_idx(std::size_t const gs, std::size_t const gid, Fn&& fn,
                              std::size_t const n, std::size_t const ds = 0)
    {
        if constexpr (impl == 0)
        {
            using DevFn = ChunkIndexFunction<Super, Fn, impl>;
            Super::fns.emplace_back(
                std::make_shared<DevFn>(gs, gid, **this, std::forward<Fn>(fn), n, ds));
        }
        else
        {
            throw std::runtime_error("no impl");
        }
        return *this;
    }



    Super& super() { return *this; }
    auto& operator*() { return super(); }
};


template<typename Boxes, typename Vectors>
struct ThreadedBoxStreamLauncher : public ThreadedStreamLauncher<Vectors>
{
    using Super = ThreadedStreamLauncher<Vectors>;
    using Box_t = Boxes::value_type;

    ThreadedBoxStreamLauncher(Vectors& vectors, Boxes const& bxes, std::size_t const nthreads = 1)
        : Super{vectors, nthreads}
        , boxes{bxes}
    {
    }


    template<std::uint16_t impl = 0, typename Fn>
    auto& async_dev(Fn&& fn)
    {
        if constexpr (impl == 0)
        {
            using DevFn = BoxStreamDeviceFunction<Box_t, Super, Fn, impl>;
            Super::fns.emplace_back(std::make_shared<DevFn>(**this, std::forward<Fn>(fn), boxes));
        }
        else
        {
            throw std::runtime_error("no impl");
        }
        return *this;
    }

    template<std::uint16_t impl = 0, typename Fn>
    auto& async_dev_idx(std::size_t const gs, std::size_t const gid, Fn&& fn)
    {
        using DevFn = BoxStreamDeviceGroupIndexFunction<Box_t, Super, Fn, impl>;
        Super::fns.emplace_back(
            std::make_shared<DevFn>(gs, gid, *this, std::forward<Fn>(fn), boxes));
        return *this;
    }

    template<std::uint16_t impl = 0, typename Fn>
    auto& async_dev_not_idx(std::size_t const gs, std::size_t const gid, Fn&& fn)
    {
        if constexpr (impl < 2)
        {
            using DevFn = BoxStreamDeviceGroupIndexFunction<Box_t, Super, Fn, impl, false>;
            Super::fns.emplace_back(
                std::make_shared<DevFn>(gs, gid, *this, std::forward<Fn>(fn), boxes));
        }
        else if (impl == 2)
        {
            using DevFn = mkn::gpu::StreamDeviceGroupIndexFunction<Super, Fn, false>;
            Super::fns.emplace_back(std::make_shared<DevFn>(gs, gid, *this, std::forward<Fn>(fn)));
        }
        else
            throw std::runtime_error("No impl");
        return *this;
    }

    Super& super() { return *this; }
    auto& operator*() { return super(); }

    Boxes const& boxes;
};

// one thread per cell
template<typename Box_t, bool sync = true>
class BoxCellLauncher : public DeviceLauncher<mkn::gpu::Guess, sync>
{
    using Super = DeviceLauncher<mkn::gpu::Guess, sync>;

public:
    BoxCellLauncher(Box_t box, size_t dev = 0)
        : Super{box.size(), dev}
    {
    }
};



} // namespace PHARE::core::gpu

#endif // PHARE_HAVE_MKN_GPU

#endif /* PHARE_CORE_UTILITIES_KERNELS_DETAIL_MKN_LAUNCHERS_HPP */
