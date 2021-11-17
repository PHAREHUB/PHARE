
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
    using bucketlist_t = BucketList<bucket_size>;
    using cell_t       = std::array<cell_index_t, dim>;
    using box_t        = Box<cell_index_t, dim>;


public:
    CellMap() = default;


    auto nbr_cells() const { return bucketsLists_.size(); }



    // add a single index to the cellmap with the specified cell
    template<typename CellIndex>
    void addToCell(CellIndex const& cell, std::size_t itemIndex);

    static auto constexpr default_extractor
        = [](auto const& item) -> decltype(auto) { return item.iCell; };
    using DefaultExtractor = decltype(default_extractor);


    // same as above but cell is found with the CellExtractor
    template<typename Array, typename CellExtractor = DefaultExtractor,
             typename = std::enable_if_t<is_iterable_v<Array>, void>>
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
    std::size_t size(cell_t cell) { return bucketsLists_[cell].size(); }

    // total number of mapped indexes
    auto size() const;

    // number of indexes mapped in the given box
    auto size(box_t const& box) const;

    // total capacity over all bucket lists
    auto capacity() const;

    // remove all bucketlists
    void clear() { bucketsLists_.clear(); }

    void empty();

    bool is_empty() const { return size() == 0; }


    float used_mem_ratio() const { return static_cast<float>(size()) / capacity(); }

    // trim all bucketlists
    void trim(std::size_t max_empty);


    // export from 'from' into 'dest' items indexed in the map found withing 'box'
    template<typename Array>
    void export_to(box_t const& box, Array const& from, Array& dest) const;


    // same as previous but applies a transformation to the items before exporting to 'dest'
    template<typename Array, typename Transformation>
    void export_to(box_t const& box, Array const& from, Array& dest, Transformation&& Fn) const;


    // export items satisfying Predicate in 'from' into 'dest'
    template<typename Array, typename Predicate>
    void export_to(Array const& from, Array& dest, Predicate&& pred) const;



    template<typename CellIndex>
    auto list_at(CellIndex const& cell);

    template<typename CellIndex>
    auto list_at(CellIndex const& cell) const;




    // item at itemIndex in items is registered at 'oldCell' but is now in a
    // different one. This method changes the cellmap to move the particle to its new cell
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



    template<typename CellIndex>
    void print(CellIndex const& cell) const;

    auto begin() { return bucketsLists_.begin(); }
    auto end() { return bucketsLists_.end(); }

private:
    struct CellHasher
    {
        template<typename CellIndex>
        std::size_t operator()(CellIndex const& cell) const
        {
            std::size_t h = 0;

            for (auto c : cell)
            {
                h ^= std::hash<int>{}(c) + 0x9e3779b9 + (h << 6) + (h >> 2);
            }
            return h;
        }
    };

    std::unordered_map<key_t, bucketlist_t, CellHasher> bucketsLists_;
};



template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename CellIndex>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::addToCell(CellIndex const& cell,
                                                                      std::size_t itemIndex)
{
    bucketsLists_[cell].add(itemIndex);
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
    auto blist = list_at(extract(items[itemIndex]));
    if (blist)
    {
        blist->get().remove(itemIndex);
    }
    else
    {
        std::cout << "item cannot be found\n";
    }
}


template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename Cells>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::reserve(Cells const& cells)
{
    for (auto const& cell : cells)
    {
        bucketsLists_.emplace(cell, bucketlist_t{});
    }
}


template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
inline auto CellMap<dim, bucket_size, cell_index_t, key_t>::size() const
{
    PHARE_LOG_SCOPE("CellMap::size()");
    auto s = 0u;
    for (auto const& [_, bucketlist] : bucketsLists_)
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
        auto blist = list_at(cell);
        if (blist)
            s += blist->get().size();
    }
    return s;
}



template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename CellIndex>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::print(CellIndex const& cell) const
{
    auto blist = list_at(cell);
    if (blist)
        for (auto itemIndex : blist->get())
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
        auto blist = list_at(cell);
        if (blist)
        {
            for (auto itemIndex : blist->get())
            {
                dest.push_back(from[itemIndex]);
            }
        }
        // TODO does not work but would be better
        // dest.insert(std::end(dest), std::begin(blist->get()), std::end(blist->get()));
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
        auto blist = list_at(cell);
        if (blist)
            for (auto itemIndex : blist->get())
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
    for (auto& [cell, blist] : bucketsLists_)
    {
        if (pred(cell))
        {
            auto blist = list_at(cell);
            if (blist)
            {
                for (auto itemIndex : blist->get())
                {
                    dest.push_back(from[itemIndex]);
                }
            }
        }
    }
}




template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename CellIndex>
inline auto CellMap<dim, bucket_size, cell_index_t, key_t>::list_at(CellIndex const& cell)
{
    // TODO consider building the particleArray from a box to have the keys
    // already and not need for the find.
    if (auto it = bucketsLists_.find(cell); it != bucketsLists_.end())
        return std::make_optional<std::reference_wrapper<bucketlist_t>>(std::ref(it->second));
    return std::optional<std::reference_wrapper<bucketlist_t>>(std::nullopt);
}


template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename CellIndex>
inline auto CellMap<dim, bucket_size, cell_index_t, key_t>::list_at(CellIndex const& cell) const
{
    // TODO consider building the particleArray from a box to have the keys
    // already and not need for the find.
    if (auto it = bucketsLists_.find(cell); it != bucketsLists_.end())
        return std::make_optional<std::reference_wrapper<const bucketlist_t>>(
            std::cref(it->second));
    return std::optional<std::reference_wrapper<const bucketlist_t>>(std::nullopt);
}



template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
inline auto CellMap<dim, bucket_size, cell_index_t, key_t>::capacity() const
{
    auto tot = 0u;
    for (auto const& [_, bucketList] : bucketsLists_)
    {
        tot += bucketList.capacity();
    }
    return tot;
}

template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::empty()
{
    for (auto& [_, bucketList] : bucketsLists_)
    {
        bucketList.empty();
    }
}



template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::trim(std::size_t max_empty)
{
    PHARE_LOG_SCOPE("CellMap::trim");
    for (auto& [_, blist] : bucketsLists_)
    {
        blist.trim(max_empty);
    }
}

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
    auto oldlist_opt = list_at(oldCell);
    if (oldlist_opt)
    {
        auto& oldList = oldlist_opt->get();
        if (oldList.in_bucket(itemIndex))
            oldList.remove(itemIndex);
    }
    add(items, itemIndex);
}




template<std::size_t dim, std::size_t bucket_size, typename cell_index_t, typename key_t>
template<typename Range, typename Predicate, typename CellExtractor>
inline auto CellMap<dim, bucket_size, cell_index_t, key_t>::partition(Range range, Predicate&& pred,
                                                                      CellExtractor extract)
{
    std::size_t shiftIdx = range.iend() - 1;
    auto pivot           = range.iend();
    for (auto& [cell, itemIndexes] : bucketsLists_)
    {
        if (!pred(cell))
        {
            for (auto idx : itemIndexes)
            {
                // partition only indexes in rangeparticles_
                if (idx >= range.ibegin() and idx < range.iend())
                {
                    assert(pivot > 0);
                    pivot--;
                    // find last element that satisfy the predicate
                    // to swap with it
                    while (!pred(key_t{extract(range.array()[shiftIdx])}) and shiftIdx > idx)
                    {
                        shiftIdx--;
                    }
                    // only swap if idx is not already the last index
                    // not satisfying the predicate
                    if (shiftIdx > idx)
                    {
                        assert(pred(key_t{extract(range.array()[shiftIdx])}));
                        assert(shiftIdx >= range.ibegin());
                        itemIndexes.updateIndex(idx, shiftIdx);
                        auto bl = list_at(extract(range.array()[shiftIdx]));
                        auto& l = bl->get();
                        l.updateIndex(shiftIdx, idx);
                        std::swap(range.array()[idx], range.array()[shiftIdx]);
                        shiftIdx--;
                    }
                }
            }
        }
    }

    auto retRange = makeRange(range.array(), range.ibegin(), range.ibegin() + pivot);
    return retRange;
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
inline void CellMap<dim, bucket_size, cell_index_t, key_t>::erase(CellIndex index, Array& items,
                                                                  CellExtractor extract)
{
    auto lastIndex    = items.size() - 1;
    auto itemsToErase = list_at(index);
    if (itemsToErase)
    {
        for (auto itemIndex : itemsToErase->get())
        {
            auto lastItemCell = extract(items[lastIndex]);
            std::swap(items[itemIndex], items[lastIndex]);

            // item previously found at lastIndex
            // is now at itemIndex
            // we need to update the index stored in
            // the bucketlist it is in
            auto lastItemList = list_at(lastItemCell);

            lastItemList->get().updateIndex(lastIndex, itemIndex);
            lastIndex--;
        }
        itemsToErase->get().empty();
        items.erase(std::begin(items) + lastIndex + 1, std::end(items));
    }
}

#if 0
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
