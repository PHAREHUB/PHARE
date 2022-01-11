
#ifndef PHARE_CELLMAP_H
#define PHARE_CELLMAP_H
#include <cstddef>
#include <string>
#include <iterator>
#include <array>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <numeric>
#include <optional>


#include "core/utilities/box/box.hpp"
#include "core/utilities/bucketlist.hpp"
#include "core/logger.hpp"
#include "core/utilities/iterators.hpp"
#include "core/utilities/meta/meta_utilities.hpp"

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


    template<typename CellIndex>
    void addToCell(CellIndex const& cell, std::size_t itemIndex)
    {
        bucketsLists_[cell].add(itemIndex);
    }

    static auto constexpr default_extractor = [](auto const& item) { return item.iCell; };
    using DefaultExtractor                  = decltype(default_extractor);

    template<typename Array, typename CellExtractor = DefaultExtractor,
             typename = std::enable_if_t<is_iterable_v<Array>, void>>
    void add(Array const& items, CellExtractor extract = default_extractor)
    {
        PHARE_LOG_SCOPE("CellMap::add (array)");
        for (std::size_t itemIndex = 0; itemIndex < items.size(); ++itemIndex)
        {
            addToCell(extract(items[itemIndex]), itemIndex);
        }
    }

    template<typename Array, typename CellExtractor = DefaultExtractor,
             typename = std::enable_if_t<is_iterable_v<Array>, void>>
    void add(Array const& items, std::size_t first, std::size_t last,
             CellExtractor extract = default_extractor)
    {
        PHARE_LOG_SCOPE("CellMap::add (iterator)");
        for (auto itemIndex = first; itemIndex <= last; ++itemIndex)
        {
            addToCell(extract(items[itemIndex]), itemIndex);
        }
    }

    template<typename Array, typename CellExtractor = DefaultExtractor,
             typename = std::enable_if_t<is_iterable_v<Array>, void>>
    void add(Array const& items, std::size_t itemIndex, CellExtractor extract = default_extractor)
    {
        PHARE_LOG_SCOPE("CellMap::add (array)");
        addToCell(extract(items[itemIndex]), itemIndex);
    }

    std::size_t size(cell_t cell) { return bucketsLists_[cell].size(); }

    template<typename Array, typename CellExtractor = DefaultExtractor,
             typename = std::enable_if_t<is_iterable_v<Array>, void>>
    void erase(Array const& items, std::size_t itemIndex, CellExtractor extract = default_extractor)
    {
        auto blist = list_at(extract(items[itemIndex]));
        if (blist)
        {
            std::cout << "removed one item\n";
            blist->get().remove(itemIndex);
        }
        else
        {
            std::cout << "item cannot be found\n";
        }
    }

    auto size() const
    {
        PHARE_LOG_SCOPE("CellMap::size()");
        auto s = 0u;
        for (auto const& [_, bucketlist] : bucketsLists_)
        {
            s += bucketlist.size();
        }
        return s;
    }

    auto size(box_t const& box) const
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


    template<typename Array>
    auto select(box_t const& box, Array const& from) const
    {
        // count all items
        std::size_t nbrTot = size(box);
        std::vector<typename Array::value_type> selection(nbrTot);
        if (nbrTot == 0)
            return selection;

        std::size_t selected_idx = 0;
        for (auto const& cell : box)
        {
            auto blist = list_at(cell);
            if (blist)
                for (auto itemIndex : blist->get())
                {
                    selection[selected_idx++] = from[itemIndex];
                }
        }
        return selection;
    }

    template<typename Array>
    void export_to(box_t const& box, Array const& from, Array& dest) const
    {
        for (auto const& cell : box)
        {
            auto blist = list_at(cell);
            if (blist)
                for (auto itemIndex : blist->get())
                {
                    dest.push_back(from[itemIndex]);
                }
            // TODO does not work but would be better
            // dest.insert(std::end(dest), std::begin(blist->get()), std::end(blist->get()));
        }
    }

    template<typename Array, typename Transformation>
    void export_to(box_t const& box, Array const& from, Array& dest, Transformation&& Fn) const
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

    template<typename CellIndex>
    auto list_at(CellIndex const& cell)
    {
        // TODO consider building the particleArray from a box to have the keys
        // already and not need for the find.
        if (auto it = bucketsLists_.find(cell); it != bucketsLists_.end())
            return std::make_optional<std::reference_wrapper<bucketlist_t>>(std::ref(it->second));
        return std::optional<std::reference_wrapper<bucketlist_t>>(std::nullopt);
    }

    template<typename CellIndex>
    auto list_at(CellIndex const& cell) const
    {
        if (auto it = bucketsLists_.find(cell); it != bucketsLists_.end())
            return std::make_optional<std::reference_wrapper<bucketlist_t const>>(
                std::ref(it->second));
        return std::optional<std::reference_wrapper<bucketlist_t const>>(std::nullopt);
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

    void clear() { bucketsLists_.clear(); }




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
        PHARE_LOG_SCOPE("CellMap::trim");
        for (auto& [_, blist] : bucketsLists_)
        {
            blist.trim(max_empty);
        }
    }

    template<typename Array, typename CellIndex, typename CellExtractor = DefaultExtractor>
    void update(Array const& items, std::size_t itemIndex, CellIndex const& oldCell,
                CellExtractor extract = default_extractor)
    {
        auto oldlist = list_at(oldCell);
        auto newlist = list_at(extract(items[itemIndex]));

        if (oldlist and newlist)
        {
            oldlist->get().remove(itemIndex);
            newlist->get().add(itemIndex);
        }
    }


    template<typename CellIndex>
    void print(CellIndex const& cell) const
    {
        auto blist = list_at(cell);
        if (blist)
            for (auto itemIndex : blist->get())
            {
                std::cout << itemIndex << "\n";
            }
    }

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

    box_t box_;
    std::array<cell_index_t, dim> shape_;
    std::unordered_map<key_t, bucketlist_t, CellHasher> bucketsLists_;
};

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
