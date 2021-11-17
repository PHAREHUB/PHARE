
#ifndef PHARE_CELLMAP_H
#define PHARE_CELLMAP_H
#include <string>
#include <iterator>
#include <array>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <numeric>


#include "core/utilities/box/box.h"
#include "core/utilities/bucketlist.h"

namespace PHARE::core
{
template<std::size_t dim, typename Type, std::size_t bucket_size, typename cell_index_t = int>
class CellMap
{
private:
    using bucketlist_t = BucketList<bucket_size, Type>;
    using cell_t       = std::array<cell_index_t, dim>;
    using box_t        = Box<cell_index_t, dim>;

    template<typename CellIndex>
    auto to_string(CellIndex const& cell) const
    {
        return std::accumulate(
            std::begin(cell), std::end(cell), std::string{},
            [](auto& s, auto const& other) { return s += std::to_string(other); });
    }

public:
    CellMap() = default;


    auto nbr_cells() const { return bucketsLists_.size(); }


    template<typename CellIndex>
    void addToCell(CellIndex const& cell, Type const& obj)
    {
        bucketsLists_[to_string(cell)].add(obj);
    }

    static auto constexpr default_extractor = [](auto const& item) { return item.iCell; };
    using DefaultExtractor                  = decltype(default_extractor);

    template<typename CellExtractor = DefaultExtractor>
    void add(std::vector<Type> const& items, CellExtractor extract = default_extractor)
    {
        for (auto const& item : items)
        {
            addToCell(extract(item), item);
        }
    }

    template<typename Iterator, typename CellExtractor = DefaultExtractor>
    void add(Iterator const& begin, Iterator const& end, CellExtractor extract = default_extractor)
    {
        for (auto it = begin; it != end; ++it)
        {
            addToCell(extract(*it), *it);
        }
    }

    std::size_t size(cell_t cell) { return bucketsLists_[to_string(cell)].size(); }

    auto size() const
    {
        auto s = 0u;
        for (auto const& [_, bucketlist] : bucketsLists_)
        {
            s += bucketlist.size();
        }
        return s;
    }


    // this assumes the given box is all within the box mapped by the CellMap
    // i.e. we do not select the objects from the intersection of the given box
    // with the map box. Any cell in the given box that does not fall into the map
    // box will result in an exception thrown when trying to get value at that key
    auto select(std::vector<Type>& items, box_t const& box) const
    {
        // count all items
        std::size_t nbrTot = 0;
        for (auto const& cell : box)
        {
            auto& blist = list_at(cell);
            nbrTot += blist.size();
        }
        std::vector<Type> selection(nbrTot);

        std::size_t item_idx = 0;
        for (auto const& cell : box)
        {
            auto blist = list_at(cell);
            for (Type const* item : blist)
            {
                selection[item_idx++] = *item;
            }
        }
        return selection;
    }


    template<typename CellIndex>
    bucketlist_t& list_at(CellIndex const& cell)
    {
        auto& c = bucketsLists_[to_string(cell)];
        return c;
    }

    template<typename CellIndex>
    bucketlist_t const& list_at(CellIndex const& cell) const
    {
        auto& c = bucketsLists_.at(to_string(cell));
        return c;
    }


    auto capacity() const
    {
        auto tot = 0u;
        for (auto const& [_, bucketList] : bucketsLists_)
        {
            tot += bucketList.capacity();
        }
        return tot;
    }


    void empty()
    {
        for (auto& [_, bucketList] : bucketsLists_)
        {
            bucketList.empty();
        }
    }


    bool is_empty() const { return size() == 0; }


    float used_mem_ratio() const { return static_cast<float>(size()) / capacity(); }

    void trim(std::size_t max_empty)
    {
        for (auto& [_, blist] : bucketsLists_)
        {
            blist.trim(max_empty);
        }
    }

private:
    box_t box_;
    std::array<cell_index_t, dim> shape_;
    std::unordered_map<std::string, bucketlist_t> bucketsLists_;
};




} // namespace PHARE::core

#endif
