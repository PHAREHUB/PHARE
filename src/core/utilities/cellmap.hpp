#ifndef PHARE_CELLMAP_H
#define PHARE_CELLMAP_H
#include <cstddef>
#include <string>
#include <iterator>
#include <array>
#include <type_traits>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <numeric>
#include <optional>

#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/indexer.hpp"
#include "core/logger.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "core/utilities/range/range.hpp"
#include "core/def.hpp"


namespace PHARE::core
{
template<std::size_t dim, typename cell_index_t = int>
class CellMap
{
private:
    using cell_t = std::array<cell_index_t, dim>;
    using box_t  = Box<cell_index_t, dim>;


public:
    CellMap(Box<cell_index_t, dim> box)
        : box_{box}
        , cellIndexes_{box.shape().template toArray<std::uint32_t>()}
    {
    }

    CellMap(CellMap const& from)            = default;
    CellMap(CellMap&& from)                 = default;
    CellMap& operator=(CellMap const& from) = default;
    CellMap& operator=(CellMap&& from)      = default;

    NO_DISCARD auto nbr_cells() const { return cellIndexes_.size(); }


    NO_DISCARD bool check_unique() const
    {
        std::unordered_map<std::size_t, std::size_t> counts;
        for (auto const& cell : box_)
        {
            auto& itemIndexes = cellIndexes_(local_(cell));
            for (auto const& index : itemIndexes)
            {
                if (counts.count(index) == 0)
                    counts[index] = 0;
                counts[index]++;
                if (counts[index] > 1)
                {
                    std::cout << "Oops " << cell << " " << index
                              << " appears mapped more than once\n";
                }
            }
        }
        bool unique = true;
        for (auto const& [index, count] : counts)
            if (count > 1)
            {
                std::cout << "CellMap CHECKUNIQUE : " << index << "appears " << count << " times\n";
                unique = false;
            }
        return unique;
    }


    // add a single index to the cellmap with the specified cell
    template<typename CellIndex>
    void addToCell(CellIndex const& cell, std::size_t itemIndex);

    static auto constexpr default_extractor = [](auto const& item) -> auto& { return item.iCell; };
    using DefaultExtractor                  = decltype(default_extractor);


    // same as above but cell is found with the CellExtractor
    template<typename Array, typename CellExtractor = DefaultExtractor,
             typename = std::enable_if_t<is_iterable_v<Array>>>
    void add(Array const& items, std::size_t itemIndex, CellExtractor extract = default_extractor);



    // add all given items indexes to the cellmap.
    // items will be added to the cell they are in, given the CellExtractor
    template<typename Array, typename CellExtractor = DefaultExtractor,
             typename = std::enable_if_t<is_iterable_v<Array>, void>>
    void add(Array const& items, CellExtractor extract = default_extractor);



    // same as above but for indexes within the given range
    template<typename Array, typename CellExtractor = DefaultExtractor,
             typename = std::enable_if_t<is_iterable_v<Array>, void>>
    void add(Array const& items, std::size_t first, std::size_t last,
             CellExtractor extract = default_extractor);



    // number of indexes stored in that cell of the cellmap
    NO_DISCARD std::size_t size(cell_t cell) const { return cellIndexes_(local_(cell)).size(); }
    NO_DISCARD std::size_t size(cell_t cell) { return cellIndexes_(local_(cell)).size(); }

    // total number of mapped indexes
    NO_DISCARD auto size() const;

    // number of indexes mapped in the given box
    NO_DISCARD auto size(box_t const& box) const;

    // total capacity over all cells
    NO_DISCARD auto capacity() const;

    // remove all indexes
    void clear()
    {
        for (auto& cell : cellIndexes_)
            cell.clear();
    }

    void empty();

    NO_DISCARD bool is_empty() const { return size() == 0; }


    NO_DISCARD float used_mem_ratio() const { return static_cast<float>(size()) / capacity(); }


    // export from 'from' into 'dest' items indexed in the map found withing 'box'
    template<typename Src, typename Dst>
    void export_to(box_t const& box, Src const& from, Dst& dest) const;


    // same as previous but applies a transformation to the items before exporting to 'dest'
    template<typename Src, typename Dst, typename Transformation>
    void export_to(box_t const& box, Src const& from, Dst& dest, Transformation&& Fn) const;


    // export items satisfying Predicate in 'from' into 'dest'
    template<typename Src, typename Dst, typename Predicate>
    void export_if(Src const& from, Dst& dest, Predicate&& pred) const;




    // item at itemIndex in items is registered at 'oldCell' but is now in a
    // different one. This method changes the cellmap to move the index to its new cell
    // new cell is obtained from the item
    template<typename Array, typename CellIndex, typename CellExtractor = DefaultExtractor>
    void update(Array& items, std::size_t itemIndex, CellIndex const& oldCell,
                CellExtractor extract = default_extractor);


    // re-orders the array so that elements satisfying the predicate are found first
    // and element not satisfying after. Returns the pivot index, which is the first
    // element not to satisfy the predicate
    // Ensures the cellmap and the re-ordered array are still consistent.
    template<typename Range, typename Predicate, typename CellExtractor = DefaultExtractor>
    auto partition(Range range, Predicate&& pred, CellExtractor = default_extractor);



    // erase all items indexed in the given range from both the cellmap and the
    // array the range is for.
    template<typename Range>
    void erase(Range&& range);


    // erase items indexes from the cellmap
    template<typename Array, typename CellExtractor = DefaultExtractor,
             typename = std::enable_if_t<is_iterable_v<Array>, void>>
    void erase(Array const& items, std::size_t itemIndex,
               CellExtractor extract = default_extractor);


    // sort all cell indexes
    void sort();

    template<typename CellIndex>
    void print(CellIndex const& cell) const;

    NO_DISCARD auto& box() { return box_; }
    NO_DISCARD auto const& box() const { return box_; }

    NO_DISCARD auto begin() { return cellIndexes_.begin(); }
    NO_DISCARD auto end() { return cellIndexes_.end(); }

    template<typename Cell>
    NO_DISCARD auto& operator()(Cell const& cell)
    {
        return cellIndexes_(local_(cell));
    }

    template<typename Cell>
    NO_DISCARD auto const& operator()(Cell const& cell) const
    {
        return cellIndexes_(local_(cell));
    }

private:
    template<typename Cell>
    auto local_(Cell const& cell) const
    {
        auto loc{cell};
        for (std::size_t i = 0; i < loc.size(); ++i)
        {
            loc[i] -= box_.lower[i];
        }
        return loc;
    }
    Box<cell_index_t, dim> box_;
    NdArrayVector<dim, Indexer> cellIndexes_;
};



template<std::size_t dim, typename cell_index_t>
template<typename CellIndex>
inline void CellMap<dim, cell_index_t>::addToCell(CellIndex const& cell, std::size_t itemIndex)
{
    if (!box_.isEmpty() and isIn(Point{cell}, box_))
        cellIndexes_(local_(cell)).add(itemIndex);
}


template<std::size_t dim, typename cell_index_t>
template<typename Array, typename CellExtractor, typename>
inline void CellMap<dim, cell_index_t>::add(Array const& items, CellExtractor extract)
{
    for (std::size_t itemIndex = 0; itemIndex < items.size(); ++itemIndex)
    {
        addToCell(extract(items[itemIndex]), itemIndex);
    }
}

template<std::size_t dim, typename cell_index_t>
template<typename Array, typename CellExtractor, typename>
inline void CellMap<dim, cell_index_t>::add(Array const& items, std::size_t first, std::size_t last,
                                            CellExtractor extract)

{
    for (auto itemIndex = first; itemIndex <= last; ++itemIndex)
    {
        addToCell(extract(items[itemIndex]), itemIndex);
    }
}


template<std::size_t dim, typename cell_index_t>
template<typename Array, typename CellExtractor, typename>
inline void CellMap<dim, cell_index_t>::add(Array const& items, std::size_t itemIndex,
                                            CellExtractor extract)
{
    addToCell(extract(items[itemIndex]), itemIndex);
}



template<std::size_t dim, typename cell_index_t>
template<typename Array, typename CellExtractor, typename>
inline void CellMap<dim, cell_index_t>::erase(Array const& items, std::size_t itemIndex,
                                              CellExtractor extract)
{
    auto& blist = cellIndexes_(local_(extract(items[itemIndex])));
    blist.remove(itemIndex);
}




template<std::size_t dim, typename cell_index_t>
inline auto CellMap<dim, cell_index_t>::size() const
{
    PHARE_LOG_SCOPE(3, "CellMap::size()");
    std::size_t s = 0;
    for (auto const& cell : cellIndexes_)
    {
        s += cell.size();
    }
    return s;
}

template<std::size_t dim, typename cell_index_t>
inline auto CellMap<dim, cell_index_t>::size(box_t const& box) const
{
    PHARE_LOG_SCOPE(3, "CellMap::size(box)");
    std::size_t s = 0;
    for (auto const& cell : box)
    {
        auto& blist = cellIndexes_(local_(cell));
        s += blist.size();
    }
    return s;
}



template<std::size_t dim, typename cell_index_t>
template<typename CellIndex>
inline void CellMap<dim, cell_index_t>::print(CellIndex const& cell) const
{
    auto& blist = cellIndexes_(local_(cell));
    for (auto itemIndex : blist)
    {
        std::cout << itemIndex << "\n";
    }
}


template<std::size_t dim, typename cell_index_t>
template<typename Src, typename Dst>
inline void CellMap<dim, cell_index_t>::export_to(box_t const& box, Src const& from,
                                                  Dst& dest) const
{
    for (auto const& cell : box)
    {
        auto& blist = cellIndexes_(local_(cell));
        for (auto itemIndex : blist)
        {
            dest.push_back(from[itemIndex]);
        }
    }
}

template<std::size_t dim, typename cell_index_t>
template<typename Src, typename Dst, typename Transformation>
inline void CellMap<dim, cell_index_t>::export_to(box_t const& box, Src const& from, Dst& dest,
                                                  Transformation&& Fn) const
{
    for (auto const& cell : box)
    {
        auto& blist = cellIndexes_(local_(cell));
        for (auto itemIndex : blist)
        {
            dest.push_back(Fn(from[itemIndex]));
        }
    }
}

template<std::size_t dim, typename cell_index_t>
template<typename Src, typename Dst, typename Predicate>
inline void CellMap<dim, cell_index_t>::export_if(Src const& from, Dst& dest,
                                                  Predicate&& pred) const
{
    for (auto const& cell : box_)
    {
        auto& blist = cellIndexes_(local_(cell));
        if (pred(cell))
        {
            for (auto itemIndex : blist)
            {
                dest.push_back(from[itemIndex]);
            }
        }
    }
}




template<std::size_t dim, typename cell_index_t>
inline auto CellMap<dim, cell_index_t>::capacity() const
{
    std::size_t tot = 0;
    for (auto const& cell : box_)
    {
        tot += cellIndexes_(local_(cell)).capacity();
    }
    return tot;
}

template<std::size_t dim, typename cell_index_t>
inline void CellMap<dim, cell_index_t>::empty()
{
    for (auto const& cell : box_)
    {
        cellIndexes_(local_(cell)).empty();
        assert(cellIndexes_(local_(cell)).is_empty());
    }
}



template<std::size_t dim, typename cell_index_t>
template<typename Array, typename CellIndex, typename CellExtractor>
inline void CellMap<dim, cell_index_t>::update(Array& items, std::size_t itemIndex,
                                               CellIndex const& oldCell, CellExtractor /*extract*/)
{
    // we want to check first if the particle is in the map
    // already. if is, needs to remove it before inserting it again
    // with at its right cell.
    // It could not be in the map at 'oldCell', in that case we want to add it.
    auto& oldList = cellIndexes_(local_(oldCell));
    if (oldList.is_indexed(itemIndex))
        oldList.remove(itemIndex);
    add(items, itemIndex);
}




template<std::size_t dim, typename cell_index_t>
template<typename Range, typename Predicate, typename CellExtractor>
inline auto CellMap<dim, cell_index_t>::partition(Range range, Predicate&& pred,
                                                  CellExtractor extract)
{
    std::size_t toSwapIndex = range.iend() - 1;
    auto pivot              = range.iend();
    for (auto const& cell : box_)
    {
        auto& itemIndexes = cellIndexes_(local_(cell));
        if (!pred(cell))
        {
            for (auto currentIdx : itemIndexes)
            {
                // partition only indexes in range
                if (currentIdx >= range.ibegin() and currentIdx < range.iend())
                {
                    assert(pivot > 0);
                    --pivot;
                    if (currentIdx < toSwapIndex)
                    {
                        while (!pred(extract(range.array()[toSwapIndex]))
                               and toSwapIndex > currentIdx)
                        {
                            --toSwapIndex;
                        }
                        if (toSwapIndex > currentIdx)
                        {
                            assert(pred(extract(range.array()[toSwapIndex])));
                            assert(toSwapIndex >= range.ibegin());
                            itemIndexes.updateIndex(currentIdx, toSwapIndex);
                            auto& l = cellIndexes_(local_(extract(range.array()[toSwapIndex])));
                            l.updateIndex(toSwapIndex, currentIdx);
                            std::swap(range.array()[currentIdx], range.array()[toSwapIndex]);
                            --toSwapIndex;
                        }
                    }
                }
            }
        }
    }

    return makeRange(range.array(), range.ibegin(), range.ibegin() + pivot);
}




template<std::size_t dim, typename cell_index_t>
template<typename Range>
inline void CellMap<dim, cell_index_t>::erase(Range&& range)
{
    auto& items = range.array();

    // first erase indexes from the cellmap
    // then items from the array
    for (std::size_t i = range.ibegin(); i < range.iend(); ++i)
    {
        erase(items, i);
    }
    items.erase(range.begin(), range.end());
}




template<std::size_t dim, typename cell_index_t>
inline void CellMap<dim, cell_index_t>::sort()
{
    for (auto const& cell : box_)
    {
        cellIndexes_(local_(cell)).sort();
    }
}


#if 0
// we don't need it now but may become handy in the future so leave it here
template<typename CellMap, typename Box, typename ParticleArray>
void get_sorted(CellMap const& cm, Box const& box, ParticleArray& dest)
{
    std::size_t pidx = 0;
    for (auto const& cell : box)
    {
        auto bucketlist = cm.list_at(cell);
        if (bucketlist)
        {
            for (auto const& part_ptr : bucketlist->get())
            {
                dest[pidx++] = *part_ptr;
            }
        }
    }
}
#endif

} // namespace PHARE::core

#endif
