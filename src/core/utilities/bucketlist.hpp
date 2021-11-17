#ifndef PHARE_BUCKETLIST_H
#define PHARE_BUCKETLIST_H

#include "core/logger.hpp"
#include "core/utilities/meta/meta_utilities.hpp"

#include <cstdint>
#include <iostream>
#include <cstddef>
#include <iterator>
#include <array>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <cassert>
#include <unordered_map>

namespace PHARE::core
{
struct BucketListIndex
{
    using Bindex                          = int; // std::uint_fast16_t;
    static constexpr Bindex default_idx_v = -1;
    Bindex bucket_idx                     = default_idx_v;
    Bindex pos                            = default_idx_v;
};



template<std::size_t bucket_size>
class BucketList
{
    // ----------------- BucketList Iterator ----------------------------------
    template<typename BucketListRef>
    class bucket_iterator : public std::iterator<std::forward_iterator_tag, std::size_t>
    {
    private:
        using Bindex = typename BucketListIndex::Bindex;

    public:
        auto& operator*() const { return bucketsList_.buckets_[curr_bucket_][curr_pos_]; }
        auto& operator*() { return bucketsList_.buckets_[curr_bucket_][curr_pos_]; }
        auto operator++();
        bool operator!=(bucket_iterator const& other) const;

        bucket_iterator(BucketListRef blist, Bindex curr_bucket = 0, Bindex curr_pos = 0)
            : curr_bucket_{curr_bucket}
            , curr_pos_{curr_pos}
            , bucketsList_{blist}
        {
        }

    private:
        typename BucketListIndex::Bindex curr_bucket_ = 0, curr_pos_ = 0;
        BucketListRef bucketsList_;
    };
    // ----------------- END BucketList Iterator ------------------------------




public:
    BucketList(std::size_t min_bucket_nbr = 1)
        : buckets_(min_bucket_nbr)
    {
        for (auto& b : buckets_)
            for (auto& v : b)
                v = BucketListIndex::default_idx_v;
    }

    using iterator       = bucket_iterator<BucketList&>;
    using const_iterator = bucket_iterator<BucketList const&>;

    void add(std::size_t itemIndex);
    void remove(std::size_t itemIndex);
    bool in_bucket(std::size_t itemIndex);

    // reallocate bucketlist memory if more empty space than max_empty
    void trim(std::size_t max_empty);

    // to use if an item in an indexed array is moved at another index
    void updateIndex(std::size_t oldIndex, std::size_t newIndex);

    // empty the bucketlist, but leaves the capacity untouched
    void empty();
    bool is_empty() const { return bucket_idx == 0 and curr == 0; }

    std::size_t size() const { return bucket_size * (bucket_idx) + curr; }
    std::size_t capacity() const { return buckets_.capacity() * bucket_size; }

    auto begin() { return iterator{*this}; }
    auto begin() const { return const_iterator{*this}; }
    auto cbegin() const { return const_iterator{*this}; }
    auto end();
    auto end() const;
    auto cend() const;



private:
    BucketListIndex lastIdx_() const;


    void decrement_()
    {
        if (curr == 0)
        {
            bucket_idx--;
            curr = bucket_size - 1;
        }
        else
        {
            auto check
                = bucket_idx >= 0
                  and bucket_idx < static_cast<typename BucketListIndex::Bindex>(buckets_.size())
                  and curr >= 0 and (curr - 1 < static_cast<BucketListIndex::Bindex>(bucket_size));
            if (!check)
            {
                std::cout << curr << " " << bucket_idx << "\n";
            }
            assert(check);
            buckets_[bucket_idx][curr - 1] = BucketListIndex::default_idx_v;
            curr--;
        }
    }

    using bucket_t                              = std::array<std::size_t, bucket_size>;
    typename BucketListIndex::Bindex bucket_idx = 0;
    typename BucketListIndex::Bindex curr       = 0;
    std::vector<bucket_t> buckets_;
};


// ----------------- BucketList Iterator ----------------------------------
template<std::size_t bucket_size>
template<typename BucketListRef>
inline auto BucketList<bucket_size>::bucket_iterator<BucketListRef>::operator++()
{
    auto val = *this;
    curr_pos_++;
    if (curr_pos_ == bucket_size)
    {
        curr_bucket_++;
        curr_pos_ = 0;
    }
    return val;
}


template<std::size_t bucket_size>
template<typename BucketListRef>
inline bool BucketList<bucket_size>::bucket_iterator<BucketListRef>::operator!=(
    bucket_iterator const& other) const
{
    return (other.curr_bucket_ != curr_bucket_ or other.curr_pos_ != curr_pos_)
           or &bucketsList_ != &other.bucketsList_;
}

// ----------------- BucketList Iterator ----------------------------------



template<std::size_t bucket_size>
void BucketList<bucket_size>::add(std::size_t itemIndex)
{
    if (curr != bucket_size)
    {
        buckets_[bucket_idx][curr++] = itemIndex;
    }
    else
    {
        if (size() == capacity())
        {
            buckets_.reserve(buckets_.size() + 1);
            buckets_.push_back(bucket_t{});
        }
        bucket_idx++;
        curr = 0;
        add(itemIndex);
    }
}

template<std::size_t bucket_size>
BucketListIndex BucketList<bucket_size>::lastIdx_() const
{
    if (curr == 0)
        return {bucket_idx - 1, bucket_size - 1};
    else
        return {bucket_idx, curr - 1};
}


template<std::size_t bucket_size>
void BucketList<bucket_size>::empty()
{
    if (capacity() > 3 * size())
    {
        trim(1);
    }
    bucket_idx = 0;
    curr       = 0;
}


template<std::size_t bucket_size>
bool BucketList<bucket_size>::in_bucket(std::size_t itemIndex)
{
    return std::find(std::begin(*this), std::end(*this), itemIndex) != std::end(*this);
}


template<std::size_t bucket_size>
void BucketList<bucket_size>::updateIndex(std::size_t oldIndex, std::size_t newIndex)
{
    auto it = std::find(std::begin(*this), std::end(*this), oldIndex);
    if (it != std::end(*this))
        *it = newIndex;
}



template<std::size_t bucket_size>
void BucketList<bucket_size>::remove(std::size_t itemIndex)
{
    iterator it = std::find(std::begin(*this), std::end(*this), itemIndex);
    if (it != std::end(*this))
    {
        auto lastIdxInBucket = lastIdx_();
        auto last_idx        = buckets_[lastIdxInBucket.bucket_idx][lastIdxInBucket.pos];
        (*it)                = last_idx;
        decrement_();
    }
}



template<std::size_t bucket_size>
auto BucketList<bucket_size>::end()
{
    // if the current cursor position is equal to bucket_size
    // it means we really are positioned on the next
    // bucket at cursor 0
    if (bucket_idx == 0 and curr == 0) // never added
        return iterator{*this, bucket_idx, curr};

    else if (curr != bucket_size)
    {
        auto it = iterator{*this, bucket_idx, curr};
        return it;
    }
    else
    {
        auto it = iterator{*this, bucket_idx + 1, 0};
        return it;
    }
}

template<std::size_t bucket_size>
auto BucketList<bucket_size>::end() const
{
    // if the current cursor position is equal to bucket_size
    // it means we really are positioned on the next
    // bucket at cursor 0
    if (curr != bucket_size)
    {
        auto it = const_iterator{*this, bucket_idx, curr};
        return it;
    }
    else
    {
        auto it = const_iterator{*this, bucket_idx + 1, 0};
        return it;
    }
}

template<std::size_t bucket_size>
auto BucketList<bucket_size>::cend() const
{
    // if the current cursor position is equal to bucket_size
    // it means we really are positioned on the next
    // bucket at cursor 0
    if (curr != bucket_size)
    {
        auto it = const_iterator{*this, bucket_idx, curr};
        return it;
    }
    else
    {
        auto it = const_iterator{*this, bucket_idx + 1, 0};
        return it;
    }
}

template<std::size_t bucket_size>
void BucketList<bucket_size>::trim(std::size_t max_empty)
{
    auto nbr_elem         = size();
    auto nbr_full_buckets = nbr_elem / bucket_size;
    auto nbr_buckets      = buckets_.size();

    auto occupied_buckets = nbr_full_buckets + ((curr == 0 or curr == bucket_size) ? 0 : 1);
    auto curr_empty       = nbr_buckets - occupied_buckets;
    if (curr_empty > max_empty)
    {
        std::vector<bucket_t> new_buckets_(std::begin(buckets_),
                                           std::begin(buckets_) + occupied_buckets + max_empty);
        buckets_.clear();
        buckets_.swap(new_buckets_);
    }
}


} // namespace PHARE::core
#endif
