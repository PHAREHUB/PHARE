#ifndef PHARE_CORE_DATA_TILES_TILING_OPTIONS_HPP
#define PHARE_CORE_DATA_TILES_TILING_OPTIONS_HPP

#include "core/utilities/types.hpp"

#include <string>
#include <algorithm>
#include <cstddef>
#include <stdexcept>

namespace PHARE::core
{

// how a patch box is cut into tiles, see tile_set_mapper.hpp
// tiles are max_tile_size wherever possible, a smaller remainder is rebalanced with the previous
// tile when below min_tile_size (min is soft: two tiles can still fall below it, max is not)
struct TilingOptions
{
    static inline std::string const min_tile_key     = "PHARE_TILING_MIN_TILE_SIZE";
    static inline std::string const max_tile_key     = "PHARE_TILING_MAX_TILE_SIZE";
    static inline std::string const min_before_split = "PHARE_TILING_MIN_BEFORE_SPLIT";

    std::size_t min_tile_size               = 4; // 4 is required for gpu dynamic shared memory
    std::size_t max_tile_size               = 10;
    std::size_t min_patch_size_before_split = 2 * min_tile_size;

    // unset env vars keep the defaults above; max defaults to at least min, so setting only the
    // min can't leave max < min
    NO_DISCARD static TilingOptions from_env()
    {
        TilingOptions opts;
        opts.min_tile_size = get_env_as(min_tile_key, opts.min_tile_size);
        opts.max_tile_size
            = get_env_as(max_tile_key, std::max(opts.min_tile_size, opts.max_tile_size));
        opts.min_patch_size_before_split = get_env_as(min_before_split, 2 * opts.min_tile_size);
        return opts.validate();
    }

    TilingOptions const& validate() const
    {
        if (min_tile_size == 0)
            throw std::runtime_error("TilingOptions: min_tile_size must be > 0");
        if (max_tile_size < min_tile_size)
            throw std::runtime_error("TilingOptions: max_tile_size < min_tile_size");
        return *this;
    }
};

} // namespace PHARE::core

#endif /*PHARE_CORE_DATA_TILES_TILING_OPTIONS_HPP*/
