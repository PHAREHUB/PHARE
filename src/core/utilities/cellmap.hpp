
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
#include "core/utilities/bucketlist.hpp"
#include "core/logger.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "core/utilities/range/range.hpp"


namespace PHARE::core
{




template<std::size_t dim, std::size_t bucket_size, typename cell_index_t = int,
         typename key_t = std::array<cell_index_t, dim>>
class CellMap
{
private:
    using bucketlist_t = BucketList;
    using cell_t       = std::array<cell_index_t, dim>;
    using box_t        = Box<cell_index_t, dim>;


public:
    CellMap(Box<cell_index_t, dim> box)
        : box_{box}
        , bucketLists_{box.shape().template toArray<std::uint32_t>()}
    {
    }

    CellMap(CellMap const& from) = default;
    CellMap(CellMap&& from)      = default;
    CellMap& operator=(CellMap const& from) = default;
    CellMap& operator=(CellMap&& from) = default;

    auto nbr_cells() const { return bucketLists_.size(); }


    bool check_unique() const
    {
        std::unordered_map<std::size_t, std::size_t> counts;
        for (auto const& cell : box_)
        {
            auto& itemIndexes = bucketLists_(local_(cell));
            for (auto const& index : itemIndexes)
            {
                // std::cout << "Checking : " << cell << " " << index << "\n";
                if (counts.count(index) == 0)
                    counts[index] = 0;
                counts[index]++;
                if (counts[index] > 1)
                {
                    std::cout << "OUPS " << cell << " " << index << "\n";
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


    template<typename Cells>
    void reserve(Cells const& cells);

    // number of indexes stored in that cell of the cellmap
    std::size_t size(cell_t cell) const { return bucketLists_(local_(cell)).size(); }
    std::size_t size(cell_t cell) { return bucketLists_(local_(cell)).size(); }

    // total number of mapped indexes
    auto size() const;

    // number of indexes mapped in the given box
    auto size(box_t const& box) const;

    // total capacity over all bucket lists
    auto capacity() const;

    // remove all bucketlists
    void clear()
    {
        for (auto& bucket : bucketLists_)
            bucket.clear();
    }

    void empty();

    bool is_empty() const { return size() == 0; }


    float used_mem_ratio() const { return static_cast<float>(size()) / capacity(); }

    // trim all bucketlists
    // void trim(std::size_t max_empty);


    // export from 'from' into 'dest' items indexed in the map found withing 'box'
    template<typename Array>
    void export_to(box_t const& box, Array const& from, Array& dest) const;


    // same as previous but applies a transformation to the items before exporting to 'dest'
    template<typename Array, typename Transformation>
    void export_to(box_t const& box, Array const& from, Array& dest, Transformation&& Fn) const;


    // export items satisfying Predicate in 'from' into 'dest'
    template<typename Array, typename Predicate>
    void export_to(Array const& from, Array& dest, Predicate&& pred) const;




    // item at itemIndex in items is registered at 'oldCell' but is now in a
    // different one. This method changes the cellmap to move the bucketListItem to its new cell
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

    // erases the elements in the (array, cellmap) belonging to the given cell index
    template<typename Array, typename CellIndex, typename CellExtractor = DefaultExtractor>
    void erase(CellIndex index, Array& items, CellExtractor extract = default_extractor);


    // erase all items indexed in the given range from both the cellmap and the
    // array the range is for.
    template<typename Range>
    void erase(Range&& range);


    // erase items indexes from the cellmap
    template<typename Array, typename CellExtractor = DefaultExtractor,
             typename = std::enable_if_t<is_iterable_v<Array>, void>>
    void erase(Array const& items, std::size_t itemIndex,
               CellExtractor extract = default_extractor);


    // sort all bucketlists
    void sort();

    template<typename CellIndex>
    void print(CellIndex const& cell) const;

    auto& box() { return box_; }
    auto const& box() const { return box_; }

    auto begin() { return bucketLists_.begin(); }
    auto end() { return bucketLists_.end(); }

    template<typename Cell>
    auto& operator()(Cell const& cell)
    {
        return bucketLists_(local_(cell));
    }

    template<typename Cell>
    auto const& operator()(Cell const& cell) const
    {
        return bucketLists_(local_(cell));
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
    NdArrayVector<dim, bucketlist_t> bucketLists_;
};



template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename CellIndex>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::addToCell(CellIndex const& cell,
                                                                      std::size_t itemIndex)
{
    if (!box_.isEmpty())
    {
        assert(isIn(Point{cell}, box_));
        if (bucketLists_(local_(cell)).in_bucket(itemIndex))
            std::cout << "ADDING : " << itemIndex << " already in bucket in cell " << Point{cell}
                      << "\n";
        if (isIn(Point{cell}, box_))
            bucketLists_(local_(cell)).add(itemIndex);
    }
}


template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename Array, typename CellExtractor, typename>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::add(Array const& items,
                                                                CellExtractor extract)
{
    PHARE_LOG_SCOPE("CellMap::add (array)");
    for (std::size_t itemIndex = 0; itemIndex < items.size(); ++itemIndex)
    {
        addToCell(extract(items[itemIndex]), itemIndex);
    }
}

template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename Array, typename CellExtractor, typename>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::add(Array const& items,
                                                                std::size_t first, std::size_t last,
                                                                CellExtractor extract)

{
    PHARE_LOG_SCOPE("CellMap::add (iterator)");
    for (auto itemIndex = first; itemIndex <= last; ++itemIndex)
    {
        addToCell(extract(items[itemIndex]), itemIndex);
    }
}


template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename Array, typename CellExtractor, typename>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::add(Array const& items,
                                                                std::size_t itemIndex,
                                                                CellExtractor extract)
{
    PHARE_LOG_SCOPE("CellMap::add (add 1 index)");
    addToCell(extract(items[itemIndex]), itemIndex);
}



template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename Array, typename CellExtractor, typename>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::erase(Array const& items,
                                                                  std::size_t itemIndex,
                                                                  CellExtractor extract)
{
    auto& blist = bucketLists_(local_(extract(items[itemIndex])));
    blist.remove(itemIndex);
}


template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename Cells>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::reserve(Cells const& cells)
{
    for (auto const& cell : cells)
    {
        // bucketLists_.emplace(cell, bucketlist_t{});
    }
}


template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
inline auto CellMap<dim, bucket_size, cell_index_t, key_t>::size() const
{
    PHARE_LOG_SCOPE("CellMap::size()");
    std::size_t s = 0;
    for (auto const& bucketlist : bucketLists_)
    {
        s += bucketlist.size();
    }
    return s;
}

template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
inline auto CellMap<dim, bucket_size, cell_index_t, key_t>::size(box_t const& box) const
{
    PHARE_LOG_SCOPE("CellMap::size(box)");
    std::size_t s = 0;
    for (auto const& cell : box)
    {
        auto& blist = bucketLists_(local_(cell));
        s += blist.size();
    }
    return s;
}



template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename CellIndex>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::print(CellIndex const& cell) const
{
    auto& blist = bucketLists_(local_(cell));
    for (auto itemIndex : blist)
    {
        std::cout << itemIndex << "\n";
    }
}


template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename Array>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::export_to(box_t const& box,
                                                                      Array const& from,
                                                                      Array& dest) const
{
    for (auto const& cell : box)
    {
        auto& blist = bucketLists_(local_(cell));
        for (auto itemIndex : blist)
        {
            dest.push_back(from[itemIndex]);
        }
    }
}

template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename Array, typename Transformation>
inline void
CellMap<dim, bucket_size, cell_index_t, key_t>::export_to(box_t const& box, Array const& from,
                                                          Array& dest, Transformation&& Fn) const
{
    for (auto const& cell : box)
    {
        auto& blist = bucketLists_(local_(cell));
        for (auto itemIndex : blist)
        {
            dest.push_back(Fn(from[itemIndex]));
        }
    }
}

template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename Array, typename Predicate>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::export_to(Array const& from,
                                                                      Array& dest,
                                                                      Predicate&& pred) const
{
    for (auto const& cell : box_)
    {
        auto& blist = bucketLists_(local_(cell));
        if (pred(cell))
        {
            for (auto itemIndex : blist)
            {
                dest.push_back(from[itemIndex]);
            }
        }
    }
}




template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
inline auto CellMap<dim, bucket_size, cell_index_t, key_t>::capacity() const
{
    std::size_t tot = 0;
    for (auto const& cell : box_)
    {
        tot += bucketLists_(local_(cell)).capacity();
    }
    return tot;
}

template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::empty()
{
    for (auto const& cell : box_)
    {
        bucketLists_(local_(cell)).empty();
        assert(bucketLists_(local_(cell)).is_empty());
    }
}


#if 0
template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::trim(std::size_t max_empty)
{
    PHARE_LOG_SCOPE("CellMap::trim");
    for (auto& [_, blist] : bucketLists_)
    {
        blist.trim(max_empty);
    }
}
#endif

template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename Array, typename CellIndex, typename CellExtractor>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::update(Array& items,
                                                                   std::size_t itemIndex,
                                                                   CellIndex const& oldCell,
                                                                   CellExtractor extract)
{
    // we want to check first if the particle is in the map
    // already. if is, needs to remove it before inserting it again
    // with at its right cell.
    // It could not be in the map at 'oldCell', in that case we want to add it.
    auto& oldList = bucketLists_(local_(oldCell));
    if (oldList.in_bucket(itemIndex))
        oldList.remove(itemIndex);
    add(items, itemIndex);
}




template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename Range, typename Predicate, typename CellExtractor>
inline auto CellMap<dim, bucket_size, cell_index_t, key_t>::partition(Range range, Predicate&& pred,
                                                                      CellExtractor extract)
{
    std::size_t shiftIdx = range.iend() - 1;
    auto pivot           = range.iend();
    for (auto const& cell : box_)
    {
        // std::cout << "cell : " << cell << "\n";
        auto& itemIndexes = bucketLists_(local_(cell));
        if (!pred(cell))
        {
            //    std::cout << "pred false, found " << itemIndexes.size()
            //              << " indexes, begin swapping in range of size " << range.size() << "\n";
            for (auto idx : itemIndexes)
            {
                //      std::cout << "idx = " << idx << "\n";
                // partition only indexes in range
                if (idx >= range.ibegin() and idx < range.iend())
                {
                    assert(pivot > 0);
                    --pivot;
                    // std::cout << "in range, --pivot : " << pivot << "\n";
                    //  find last element that satisfy the predicate
                    //  to swap with it
                    //         std::cout << "looking for item to swap from " << shiftIdx << "\n";
                    while (!pred(key_t{extract(range.array()[shiftIdx])}) and shiftIdx > idx)
                    {
                        --shiftIdx;
                    }
                    //        std::cout << "found " << shiftIdx << " with pred true, swapping...\n";
                    // only swap if idx is not already the last index
                    // not satisfying the predicate
                    if (shiftIdx > idx)
                    {
                        assert(pred(key_t{extract(range.array()[shiftIdx])}));
                        assert(shiftIdx >= range.ibegin());
                        itemIndexes.updateIndex(idx, shiftIdx);
                        auto& l = bucketLists_(local_(extract(range.array()[shiftIdx])));
                        l.updateIndex(shiftIdx, idx);
                        std::swap(range.array()[idx], range.array()[shiftIdx]);
                        --shiftIdx;
                    }
                }
            }
        }
    }

    return makeRange(range.array(), range.ibegin(), range.ibegin() + pivot);
}




template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename Range>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::erase(Range&& range)
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




template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename Array, typename CellIndex, typename CellExtractor>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::erase(CellIndex cell, Array& items,
                                                                  CellExtractor extract)
{
    auto lastIndex     = items.size() - 1;
    auto& itemsToErase = bucketLists_(local_(cell));
    assert(false);
    //  for (auto itemIndex : itemsToErase)
    //   {
    //       auto lastItemCell = extract(items[lastIndex]);
    //       std::swap(items[itemIndex], items[lastIndex]);

    //       // item previously found at lastIndex
    //       // is now at itemIndex
    //       // we need to update the index stored in
    //       // the bucketlist it is in
    //       auto& lastItemList = bucketLists_(local_(lastItemCell));

    //       lastItemList.updateIndex(lastIndex, itemIndex);
    //       --lastIndex;
    //  }
    //  itemsToErase.empty();
    //  items.erase(std::begin(items) + lastIndex + 1, std::end(items));
}


template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::sort()
{
    for (auto const& cell : box_)
    {
        bucketLists_(local_(cell)).sort();
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
