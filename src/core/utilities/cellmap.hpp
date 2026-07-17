#ifndef PHARE_CELLMAP_H
#define PHARE_CELLMAP_H

#include "core/data/ndarray/ndarray_view.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/def.hpp"
#include "core/def/phare_config.hpp"
#include "core/utilities/indexer.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

#include <array>
#include <cstddef>
#include <type_traits>
#include <unordered_map>


namespace PHARE::core
{

struct CellMapOptions
{
    using value_type = int;

    std::size_t dimension;
    AllocatorMode alloc_mode = AllocatorMode::CPU;
    StorageMode storage_mode = StorageMode::VECTOR;

    CellMapOptions static FROM(auto opts)
    {
        return {opts.dimension, opts.alloc_mode, opts.storage_mode};
    }
};


template<auto storage_mode, auto opts>
struct CellMapStorage;


// SPAN — non-owning view; Indexer<SPAN> objects live inside the VECTOR storage
template<auto opts>
struct CellMapStorage<StorageMode::SPAN, opts>
{
    bool constexpr static c_order = true;
    using Options                 = decltype(opts);
    using cell_index_t            = Options::value_type;

    CellMapStorage() = default;

    using idx_span_t   = Indexer<IndexerOptions{opts.alloc_mode, StorageMode::SPAN}>;
    using idx_vector_t = Indexer<IndexerOptions{opts.alloc_mode, StorageMode::VECTOR}>;

    CellMapStorage(auto& cellmap)
        : box_{cellmap.box_}
        , cellIndexes_{static_cast<idx_span_t*>(cellmap.cellIndexes_.data()),
                       cellmap.cellIndexes_.shape()}
    {
    }

    void reset_from(auto& cellmap)
    {
        cellIndexes_ = NdArrayViewSpan<opts.dimension, idx_span_t, idx_vector_t>{
            static_cast<idx_span_t*>(cellmap.cellIndexes_.data()), cellmap.cellIndexes_.shape()};
    }

    Box<cell_index_t, opts.dimension> box_;
    NdArrayViewSpan<opts.dimension, idx_span_t, idx_vector_t> cellIndexes_;
};

// VECTOR — owns the single NdArrayVector; its Indexer<VECTOR> elements inherit Indexer<SPAN>,
// so the SPAN base above can alias them directly via static_cast with no extra allocation.
template<auto opts>
struct CellMapStorage<StorageMode::VECTOR, opts>
    : CellMapStorage<StorageMode::SPAN,
                     CellMapOptions{opts.dimension, opts.alloc_mode, StorageMode::SPAN}>
{
    bool constexpr static c_order = true;
    using Options                 = decltype(opts);
    using cell_index_t            = Options::value_type;
    using SpanBase
        = CellMapStorage<StorageMode::SPAN,
                         CellMapOptions{opts.dimension, opts.alloc_mode, StorageMode::SPAN}>;

    using idx_span_t   = Indexer<IndexerOptions{opts.alloc_mode, StorageMode::SPAN}>;
    using idx_vector_t = Indexer<IndexerOptions{opts.alloc_mode, StorageMode::VECTOR}>;

    CellMapStorage(Box<cell_index_t, opts.dimension> const& box)
        : SpanBase{}
        , cellIndexes_{box.shape().template toArray<std::uint32_t>()}
    {
        SpanBase::box_ = box;
        SpanBase::cellIndexes_
            = {static_cast<idx_span_t*>(cellIndexes_.data()), cellIndexes_.shape()};
    }

    NdArrayVector<opts.dimension, idx_vector_t, c_order, opts.alloc_mode> cellIndexes_;
};


template<auto opts>
class CellMap : public CellMapStorage<opts.storage_mode, opts>
{
    using Super = CellMapStorage<opts.storage_mode, opts>;
    using Super::box_;
    using Super::cellIndexes_;
    using Options      = decltype(opts);
    using cell_index_t = Options::value_type;

    template<auto, auto>
    friend struct CellMapStorage;

private:
    using cell_t = std::array<cell_index_t, opts.dimension>;
    using box_t  = Box<cell_index_t, opts.dimension>;

public:
    CellMap(Box<cell_index_t, opts.dimension> const& box)
        requires(opts.storage_mode == StorageMode::VECTOR)
        : Super{box}
    {
    }

    CellMap(CellMap const& from)            = default;
    CellMap(CellMap&& from)                 = default;
    CellMap& operator=(CellMap const& from) = default;
    CellMap& operator=(CellMap&& from)      = default;


    template<typename... Args>
    CellMap(Args&&... args)
        requires std::is_constructible_v<Super, Args&&...>
    _PHARE_ALL_FN_ : Super{std::forward<Args>(args)...}
    {
    }

    // produce a non-owning span view of this CellMap for GPU kernel dispatch.
    // Syncs every Indexer's span before building the view so GPU sees valid pointers.
    auto view()
        requires(opts.storage_mode == StorageMode::VECTOR)
    {
        for (auto& idx : cellIndexes_)
            idx.sync();
        constexpr auto span_opts
            = CellMapOptions{opts.dimension, opts.alloc_mode, StorageMode::SPAN};
        return CellMap<span_opts>{*this};
    }

    // refresh span pointer after source CellMap's NdArrayVector may have reallocated
    void reset_from(auto& from)
        requires(opts.storage_mode == StorageMode::SPAN)
    {
        Super::reset_from(from);
    }



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

    static auto constexpr default_extractor
        = [](auto const& item) -> auto& { return item.iCell(); };
    using DefaultExtractor = decltype(default_extractor);


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


    void erase(auto& items, box_t const& box);


    // erase all items indexed in the given range from both the cellmap and the
    // array the range is for.
    template<typename Range>
    void erase(Range range);


    // erase items indexes from the cellmap
    template<typename Array, typename CellExtractor = DefaultExtractor,
             typename = std::enable_if_t<is_iterable_v<Array>, void>>
    void erase(Array const& items, std::size_t itemIndex,
               CellExtractor extract = default_extractor);


    template<typename Array, typename = std::enable_if_t<is_iterable_v<Array>, void>>
    void erase(Array const& items, cell_t const oldCell, std::size_t const itemIndex);


    // sort all cell indexes
    void sort();

    template<typename CellIndex>
    void print(CellIndex const& cell) const;

    NO_DISCARD auto& box() { return box_; }
    NO_DISCARD auto const& box() const { return box_; }

    NO_DISCARD auto begin() { return cellIndexes_.begin(); }
    NO_DISCARD auto end() { return cellIndexes_.end(); }


    template<template<typename, std::size_t> typename Arr>
    NO_DISCARD auto& operator()(Arr<std::uint32_t, opts.dimension> const& arr)
    {
        return cellIndexes_(arr);
    }

    template<template<typename, std::size_t> typename Arr>
    NO_DISCARD auto& operator()(Arr<std::uint32_t, opts.dimension> const& arr) const
    {
        return cellIndexes_(arr);
    }

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

    void swap(cell_t const a, cell_t const b, std::size_t const i0, std::size_t i1);

private:
    template<typename Cell>
    auto local_(Cell const& cell) const
    {
        return (Point{cell} - box_.lower).as_unsigned();
    }


    // Box<cell_index_t, dim> box_;
    // bool constexpr static c_order = true;
    // NdArrayVector<dim, Indexer, c_order, default_allocator_mode()> cellIndexes_;
};


template<auto opts>
void CellMap<opts>::swap(cell_t const a, cell_t const b, std::size_t const i0, std::size_t i1)
{
    assert(cellIndexes_(local_(a)).is_indexed(i0));
    assert(cellIndexes_(local_(b)).is_indexed(i1));
    cellIndexes_(local_(a)).updateIndex(i0, i1);
    cellIndexes_(local_(b)).updateIndex(i1, i0);
}


template<auto opts>
template<typename CellIndex>
inline void CellMap<opts>::addToCell(CellIndex const& cell, std::size_t itemIndex)
{
    if (isIn(Point{cell}, box_))
        cellIndexes_(local_(cell)).add(itemIndex);
}


template<auto opts>
template<typename Array, typename CellExtractor, typename>
inline void CellMap<opts>::add(Array const& items, CellExtractor extract)
{
    for (std::size_t itemIndex = 0; itemIndex < items.size(); ++itemIndex)
    {
        addToCell(extract(items[itemIndex]), itemIndex);
    }
}

template<auto opts>
template<typename Array, typename CellExtractor, typename>
inline void CellMap<opts>::add(Array const& items, std::size_t first, std::size_t last,
                               CellExtractor extract)

{
    for (auto itemIndex = first; itemIndex <= last; ++itemIndex)
    {
        addToCell(extract(items[itemIndex]), itemIndex);
    }
}


template<auto opts>
template<typename Array, typename CellExtractor, typename>
inline void CellMap<opts>::add(Array const& items, std::size_t itemIndex, CellExtractor extract)
{
    addToCell(extract(items[itemIndex]), itemIndex);
}



template<auto opts>
template<typename Array, typename CellExtractor, typename>
inline void CellMap<opts>::erase(Array const& items, std::size_t itemIndex, CellExtractor extract)
{
    auto& blist = cellIndexes_(local_(extract(items[itemIndex])));
    blist.remove(itemIndex);
}


template<auto opts>
template<typename Array, typename>
inline void CellMap<opts>::erase(Array const& items, cell_t const oldCell,
                                 std::size_t const itemIndex)
{
    auto& blist = cellIndexes_(local_(oldCell));
    blist.remove(itemIndex);
}


template<auto opts>
inline auto CellMap<opts>::size() const
{
    PHARE_LOG_SCOPE(3, "CellMap::size()");
    std::size_t s = 0;
    for (auto const& cell : cellIndexes_)
    {
        s += cell.size();
    }
    return s;
}

template<auto opts>
inline auto CellMap<opts>::size(box_t const& box) const
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



template<auto opts>
template<typename CellIndex>
inline void CellMap<opts>::print(CellIndex const& cell) const
{
    auto& blist = cellIndexes_(local_(cell));
    for (auto itemIndex : blist)
    {
        std::cout << itemIndex << "\n";
    }
}


template<auto opts>
template<typename Src, typename Dst>
inline void CellMap<opts>::export_to(box_t const& box, Src const& from, Dst& dest) const
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

template<auto opts>
template<typename Src, typename Dst, typename Transformation>
inline void CellMap<opts>::export_to(box_t const& box, Src const& from, Dst& dest,
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

template<auto opts>
template<typename Src, typename Dst, typename Predicate>
inline void CellMap<opts>::export_if(Src const& from, Dst& dest, Predicate&& pred) const
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




template<auto opts>
inline auto CellMap<opts>::capacity() const
{
    std::size_t tot = 0;
    for (auto const& cell : box_)
    {
        tot += cellIndexes_(local_(cell)).capacity();
    }
    return tot;
}

template<auto opts>
inline void CellMap<opts>::empty()
{
    for (auto const& cell : box_)
    {
        cellIndexes_(local_(cell)).empty();
        assert(cellIndexes_(local_(cell)).is_empty());
    }
}



template<auto opts>
template<typename Array, typename CellIndex, typename CellExtractor>
inline void CellMap<opts>::update(Array& items, std::size_t itemIndex, CellIndex const& oldCell,
                                  CellExtractor /*extract*/)
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




template<auto opts>
template<typename Range, typename Predicate, typename CellExtractor>
inline auto CellMap<opts>::partition(Range range, Predicate&& pred, CellExtractor extract)
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

    return makeRange(range.array(), range.ibegin(), pivot);
}


template<auto opts>
void CellMap<opts>::erase(auto& items, box_t const& box)
{
    std::size_t prev_idx = 0; // to throw if not sorted

    for (auto const& cell : box)
    {
        auto& blist = cellIndexes_(local_(cell));
        if (blist.size() == 0)
            continue;
        prev_idx = blist.back();
        for (std::size_t i = blist.size(); i-- > 0;) // reverse loop
        {
            auto const& idx = blist[i];
            PHARE_LOG_LINE_SS(i << " " << idx << " " << items.size());
            if (idx > prev_idx)
                throw std::runtime_error("Index list is not sorted");
            if (idx > items.size())
                throw std::runtime_error("Index list is invalid");
            items.erase(items.begin() + idx);
            prev_idx = idx;
        }
        blist.clear();
    }
}


template<auto opts>
template<typename Range>
inline void CellMap<opts>::erase(Range range)
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




template<auto opts>
inline void CellMap<opts>::sort()
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
