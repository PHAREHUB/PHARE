#ifndef PHARE_INDEXER_H
#define PHARE_INDEXER_H

#include "core/utilities/types.hpp"

#include <cassert>
#include <cstdint>
#include <iostream>
#include <cstddef>
#include <iterator>
#include <array>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <tuple>

namespace PHARE::core
{


class Indexer
{
public:
    Indexer()                                = default;
    Indexer(Indexer const& other)            = default;
    Indexer(Indexer&& other)                 = default;
    Indexer& operator=(Indexer const& other) = default;
    Indexer& operator=(Indexer&& other)      = default;

    void add(std::size_t itemIndex) { indexes_.push_back(itemIndex); }
    void remove(std::size_t itemIndex)
    {
        auto it = std::find(std::begin(indexes_), std::end(indexes_), itemIndex);
        if (it != std::end(indexes_))
        {
            indexes_.erase(it);
        }
        assert(!is_indexed(itemIndex));
    }
    [[nodiscard]] bool is_indexed(std::size_t itemIndex)
    {
        return std::end(indexes_) != std::find(std::begin(indexes_), std::end(indexes_), itemIndex);
    }

    // reallocate bucketlist memory if more empty space than max_empty
    // void trim(std::size_t max_empty);

    // to use if an item in an indexed array is moved at another index
    void updateIndex(std::size_t oldIndex, std::size_t newIndex)
    {
        //
        auto it = std::find(std::begin(indexes_), std::end(indexes_), oldIndex);
        if (it != std::end(indexes_))
        {
            *it = newIndex;
        }
    }

    // empty the bucketlist, but leaves the capacity untouched
    void empty() { indexes_.resize(0); };
    [[nodiscard]] bool is_empty() const { return indexes_.size() == 0; }

    [[nodiscard]] std::size_t size() const { return indexes_.size(); }
    [[nodiscard]] std::size_t capacity() const { return indexes_.capacity(); }

    [[nodiscard]] auto begin() { return indexes_.begin(); }
    [[nodiscard]] auto begin() const { return indexes_.begin(); }
    [[nodiscard]] auto cbegin() const { return indexes_.begin(); }
    [[nodiscard]] auto end() { return indexes_.end(); }
    [[nodiscard]] auto end() const { return indexes_.end(); }
    [[nodiscard]] auto cend() const { return indexes_.end(); }
    void sort() { std::sort(indexes_.begin(), indexes_.end()); }
    void clear() { indexes_.clear(); }


private:
    std::vector<std::size_t> indexes_;
};


// ==========================================================




} // namespace PHARE::core
#endif
