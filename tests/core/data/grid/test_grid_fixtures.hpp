#ifndef PHARE_TEST_CORE_DATA_TEST_GRID_FIXTURES_HPP
#define PHARE_TEST_CORE_DATA_TEST_GRID_FIXTURES_HPP

#include "core/data/grid/grid.hpp"
#include "core/utilities/equality.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/data/field/field_box.hpp"

#include "tests/core/data/field/test_field_fixtures.hpp"

namespace PHARE::core
{


template<typename... T0s, typename... T1s>
EqualityReport compare_reduced_fields(GridTileSet<T0s...> const& ref, Grid<T1s...> const& cmp,
                                      double const diff = 1e-15)
{
    if (!ref.size() != cmp.size())
        return EqualityReport{false, "Grid/GridTileset shape/size mismatch"};
    auto tmp = cmp;
    tmp.zero();
    reduce_into(tmp, ref);
    return compare_fields(*tmp, *cmp, diff);
}


template<typename... T0s, typename... T1s>
EqualityReport compare_reduced_fields(Grid<T0s...> const& ref, GridTileSet<T1s...> const& cmp,
                                      double const diff = 1e-15)
{
    if (ref.size() != cmp.size())
        return EqualityReport{false, "Grid/GridTileset shape/size mismatch"};
    auto tmp = ref;
    tmp.zero();
    reduce_into(tmp, cmp);
    return compare_fields(*ref, *tmp, diff);
}


template<typename... T0s, typename... T1s>
EqualityReport compare_fields(GridTileSet<T0s...> const& ref, Grid<T1s...> const& cmp,
                              double const diff = 1e-15)
{
    if (ref.size() != cmp.size())
        return EqualityReport{false, "Grid/GridTileset shape/size mismatch"};
    auto tmp = cmp;
    tmp.zero();
    reduce_single(tmp, ref);
    return compare_fields(*tmp, *cmp, diff);
}


template<typename... T0s, typename... T1s>
EqualityReport compare_fields(Grid<T0s...> const& ref, GridTileSet<T1s...> const& cmp,
                              double const diff = 1e-15)
{
    if (ref.size() != cmp.size())
        return EqualityReport{false, "Grid/GridTileset shape/size mismatch"};
    auto tmp = ref;
    tmp.zero();
    reduce_single(tmp, cmp);
    return compare_fields(*ref, *tmp, diff);
}


template<typename... T0s, typename... T1s>
EqualityReport compare_reduced_fields(Grid<T0s...> const& ref, Grid<T1s...> const& cmp,
                                      double const diff = 1e-15)
{
    if (ref.size() != cmp.size())
        return EqualityReport{false, "Grid shape/size mismatch"};
    return compare_fields(*ref, *cmp, diff);
}


template<typename GridLayout_t, auto opts>
void zero_ghost_layer(basic::Field<opts>& f, GridLayout_t const& layout)
{
    for (auto const& ghost : layout.ghostBoxFor(f).remove(layout.domainBoxFor(f)))
        for (auto const& bix : ghost)
            f(bix) = 0;
}

template<typename GridLayout_t, typename Grid_t, typename Field_t>
void zero_ghost_layer(FieldTileSet<GridLayout_t, Grid_t, Field_t>& ts, GridLayout_t const& layout)
{
    auto const pq = ts.physicalQuantity();
    for (auto const& ghost : layout.ghostBoxFor(pq).remove(layout.domainBoxFor(pq)))
        FieldBox{ts, layout, ghost}.op(0);
}

} // namespace PHARE::core


#endif /*PHARE_TEST_CORE_DATA_TEST_GRID_FIXTURES_HPP*/
