#ifndef PHARE_CORE_DATA_PARTICLES_BISECTION_RANGE_MAPPER_HPP
#define PHARE_CORE_DATA_PARTICLES_BISECTION_RANGE_MAPPER_HPP

#include <cstddef>

#include "core/logger.hpp"
#include "core/vector.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"

namespace PHARE::core
{
template<typename Box_t>
struct CellFlattener
{
    template<typename Icell>
    std::size_t operator()(Icell const& icell) const
    {
        if constexpr (Box_t::dimension == 3)
            return icell[2] + icell[1] * box.upper[2] + icell[0] * box.upper[1] * box.upper[0];

        return 1;
    }

    Box_t const& box;
};

template<typename ParticleArray, typename Box_t>
class BisectionRangeMapper
{
    static constexpr auto dim = ParticleArray::dimension;
    using Ranges              = std::vector<std::pair<std::size_t*, std::size_t*>>;

public:
    BisectionRangeMapper(ParticleArray const& particles_, Box_t const& box_)
        : particles{particles_}
        , box{box_}
    {
    }

    auto map() const;

private:
    template<typename Icell>
    std::size_t flat_cell(Icell const& icell) const
    {
        return flattener(icell);
    }

    std::size_t find_lowest_of(std::size_t const&, std::size_t const&, std::size_t const&,
                               std::size_t const&, std::size_t const&) const;
    std::size_t find_highest_of(std::size_t const&, std::size_t const&, std::size_t const&,
                                std::size_t const&, std::size_t const&) const;

    ParticleArray const& particles;
    Box_t const& box;
    CellFlattener<Box_t> const flattener{box};

    std::size_t const lo_flat{flat_cell(box.lower())};
    std::size_t const up_flat{flat_cell(box.upper())};
};



template<typename Ps, typename Box_t>
std::size_t BisectionRangeMapper<Ps, Box_t>::find_highest_of(std::size_t const& flell,
                                                             std::size_t const& lo, //
                                                             std::size_t const& up,
                                                             std::size_t const& lidx,
                                                             std::size_t const& uidx) const
{
    std::size_t midx = (uidx - lidx) / 2;
    auto const& fmid = flat_cell(particles.iCell(midx + lidx));

    std::stringstream ss;
    ss << " " << flell << " " << lidx << " " << uidx << " " << midx << " " << fmid;
    PHARE_LOG_LINE_STR(ss.str());

    if (flell == fmid)
    {
        PHARE_LOG_LINE_STR("==");
        if (up == up_flat and flat_cell(particles.iCell(particles.size() - 1)) == flell)
            return uidx;
        return find_highest_of(flell, lo, up, lidx, uidx); // bug
    }

    if (flell < fmid)
    {
        PHARE_LOG_LINE_STR("<");
        return find_highest_of(flell, lo, up, lidx, midx);
    }
    PHARE_LOG_LINE_STR(">");
    return find_highest_of(flell, lo, up, lidx + midx, uidx);
}

template<typename Ps, typename Box_t>
std::size_t BisectionRangeMapper<Ps, Box_t>::find_lowest_of(std::size_t const& flell,
                                                            std::size_t const& lo, //
                                                            std::size_t const& up,
                                                            std::size_t const& lidx,
                                                            std::size_t const& uidx) const
{
    std::size_t midx = (uidx - lidx) / 2;
    auto const& fmid = flat_cell(particles.iCell(midx + lidx));

    std::stringstream ss;
    ss << " " << flell << " " << lidx << " " << uidx << " " << midx;
    PHARE_LOG_LINE_STR(ss.str());

    if (flell == fmid)
    {
        PHARE_LOG_LINE_STR("==");
        if (lo == lo_flat and flat_cell(particles.iCell(0)) == flell)
            return lidx;
        return find_lowest_of(flell, lo, up, lidx, uidx); // bug
    }

    if (flell < fmid)
    {
        PHARE_LOG_LINE_STR("<");
        return find_lowest_of(flell, lo, up, lidx, midx);
    }
    PHARE_LOG_LINE_STR(">");
    return find_lowest_of(flell, lo, up, lidx + midx, uidx);
}



template<typename Ps, typename Box_t>
auto BisectionRangeMapper<Ps, Box_t>::map() const
{
    static_assert(dim == 3); // for now

    auto lowest = find_lowest_of(flat_cell(box.lower()), lo_flat, up_flat, 0, particles.size() - 1);
    PHARE_LOG_LINE_STR(lowest);

    auto highest
        = find_highest_of(flat_cell(box.upper()), lo_flat, up_flat, lowest, particles.size() - 1);
    PHARE_LOG_LINE_STR(highest);

    // ranges are between lowest and highest
    Ranges ranges;

    return ranges;
}

} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_BISECTION_RANGE_MAPPER_HPP */
