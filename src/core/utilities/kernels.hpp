#ifndef PHARE_CORE_UTILITIES_KERNELS_HPP
#define PHARE_CORE_UTILITIES_KERNELS_HPP

#include "core/def/phare_config.hpp"
#include "core/utilities/kernels/launchers.hpp"

#if PHARE_HAVE_GPU || PHARE_HAVE_MKN_GPU

#include "core/utilities/box/box.hpp"

namespace PHARE::kernel
{

template<std::size_t dim, typename... Args>
void launch(core::Box<std::uint32_t, dim> const& box, std::size_t const& size, Args&&... args)
{
    if constexpr (CompileOptions::WithRAJA)
    {
        throw std::runtime_error("launch::kernel NO IMPL");
    }
    else if (CompileOptions::WithMknGpu)
    {
        PHARE_WITH_MKN_GPU(
            { core::gpu::BoxCellNLauncher<core::Box<std::uint32_t, dim>>{box, size}(args...); });
    }
    else
        throw std::runtime_error("launch::kernel NO ALTERNATIVE");
}

template<typename... Args>
void launch(std::size_t const& size, Args&&... args)
{
    if constexpr (CompileOptions::WithRAJA)
    {
        throw std::runtime_error("launch::kernel NO IMPL");
    }
    else if (CompileOptions::WithMknGpu)
    {
        PHARE_WITH_MKN_GPU({ core::gpu::DeviceLauncher<mkn::gpu::Guess, true>{size}(args...); });
    }
    else
        throw std::runtime_error("launch::kernel NO ALTERNATIVE");
}

// same as above, on a stream, waits for that stream only
template<typename Stream, typename... Args>
    requires(!std::is_convertible_v<Stream const&, std::size_t>)
void launch(Stream& stream, std::size_t const& size, Args&&... args)
{
    if constexpr (CompileOptions::WithRAJA)
    {
        throw std::runtime_error("launch::kernel NO IMPL");
    }
    else if (CompileOptions::WithMknGpu)
    {
        PHARE_WITH_MKN_GPU({
            core::gpu::DeviceLauncher<mkn::gpu::Guess, true>{size}.stream(stream, args...);
        });
    }
    else
        throw std::runtime_error("launch::kernel NO ALTERNATIVE");
}


} // namespace PHARE::kernel

#else

namespace PHARE::kernel
{
template<typename... Args>
void launch(Args&&... args) // noop
{
    throw std::runtime_error("nooop");
}

} // namespace PHARE::kernel

#endif // PHARE_HAVE_GPU || PHARE_HAVE_MKN_GPU

#endif /* PHARE_CORE_UTILITIES_KERNELS_HPP */
