
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

namespace PHARE::core
{
template<std::size_t dim, typename Type, std::size_t bucket_size, typename cell_index_t = int,
         typename key_t = std::array<cell_index_t, dim>>
class CellMap
{
private:
    using bucketlist_t = BucketList<bucket_size, Type>;
    using cell_t       = std::array<cell_index_t, dim>;
    using box_t        = Box<cell_index_t, dim>;


public:
    CellMap() = default;


    auto nbr_cells() const { return bucketsLists_.size(); }


    template<typename CellIndex>
    void addToCell(CellIndex const& cell, Type& obj)
    {
        bucketsLists_[cell].add(obj);
    }

    static auto constexpr default_extractor = [](auto const& item) { return item.iCell; };
    using DefaultExtractor                  = decltype(default_extractor);

    template<typename Array, typename CellExtractor = DefaultExtractor>
    void add(Array& items, CellExtractor extract = default_extractor)
    {
        PHARE_LOG_SCOPE("CellMap::add (array)");
        for (auto& item : items)
        {
            addToCell(extract(item), item);
        }
    }

    template<typename Iterator, typename CellExtractor = DefaultExtractor>
    void add(Iterator const& begin, Iterator const& end, CellExtractor extract = default_extractor)
    {
        PHARE_LOG_SCOPE("CellMap::add (iterator)");
        for (auto it = begin; it != end; ++it)
        {
            addToCell(extract(*it), *it);
        }
    }

    std::size_t size(cell_t cell) { return bucketsLists_[cell].size(); }

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


    auto select(box_t const& box) const
    {
        // count all items
        std::size_t nbrTot = size(box);
        std::vector<Type> selection(nbrTot);
        if (nbrTot == 0)
            return selection;

        std::size_t item_idx = 0;
        for (auto const& cell : box)
        {
            auto blist = list_at(cell);
            if (blist)
                for (Type const* item : blist->get())
                {
                    selection[item_idx++] = *item;
                }
        }
        return selection;
    }

    template<typename Container>
    void export_to(box_t const& box, Container& dest) const
    {
        for (auto const& cell : box)
        {
            auto blist = list_at(cell);
            if (blist)
                for (Type const* item : blist->get())
                {
                    dest.push_back(*item);
                }
            // TODO does not work but would be better
            // dest.insert(std::end(dest), std::begin(blist->get()), std::end(blist->get()));
        }
    }

    template<typename Container, typename Transformation>
    void export_to(box_t const& box, Container& dest, Transformation&& Fn) const
    {
        for (auto const& cell : box)
        {
            auto blist = list_at(cell);
            if (blist)
                for (Type const* item : blist->get())
                {
                    // TODO if particleArray this does clean_=false
                    // for each element... but cannot do insert()
                    // since there is no insert with Op
                    // a bit silly....
                    dest.push_back(Fn(*item));
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

    template<typename CellIndex, typename CellExtractor = DefaultExtractor>
    void update(Type& item, CellIndex const& oldCell, CellExtractor extract = default_extractor)
    {
        auto oldlist = list_at(oldCell);
        auto newlist = list_at(extract(item));

        if (oldlist and newlist)
        {
            oldlist->get().remove(item);
            newlist->get().add(item);
        }
    }


    template<typename CellIndex>
    void print(CellIndex const& cell) const
    {
        auto blist = list_at(cell);
        if (blist)
            for (auto const& item : blist->get())
            {
                std::cout << *item << "\n";
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


} // namespace PHARE::core

#endif
