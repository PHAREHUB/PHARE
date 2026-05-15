#ifndef PHARE_OPTIONS_MHD_OPTIONS_HPP
#define PHARE_OPTIONS_MHD_OPTIONS_HPP


#include <cstddef>

namespace PHARE
{

template<auto opts>
struct MHDFieldOptions
{
    auto static constexpr dimension          = opts.dimension;
    std::size_t static constexpr ghost_width = 4;
};


} // namespace PHARE

#endif // PHARE_OPTIONS_MHD_OPTIONS_HPP
