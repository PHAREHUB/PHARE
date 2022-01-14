#ifndef PHARE_BUCKETLIST_H
#define PHARE_BUCKETLIST_H

#include "core/utilities/meta/meta_utilities.hpp"
#include <cstdint>
#include <iostream>
#include <cstddef>
#include <iterator>
#include <array>
#include <stdexcept>
#include <type_traits>
#include <vector>
namespace PHARE::core
{
struct BucketListIndex
{
    using Bindex                          = std::uint_fast16_t;
    static constexpr Bindex default_idx_v = 0;
    Bindex bucket_idx                     = default_idx_v;
    Bindex pos                            = default_idx_v;
};


struct BucketListItem
{
    BucketListIndex index;
};




template<std::size_t bucket_size>
class BucketList
{
    class iterator : public std::iterator<std::forward_iterator_tag, std::size_t>
    {
    private:
        using Bindex = typename BucketListIndex::Bindex;

    public:
        auto operator*() const { return bucketsList_.buckets_[curr_bucket_][curr_pos_]; }
        auto operator*() { return bucketsList_.buckets_[curr_bucket_][curr_pos_]; }
        iterator operator++();
        bool operator!=(iterator const& other) const;

        iterator(BucketList const& blist, Bindex curr_bucket = 0, Bindex curr_pos = 0)
            : curr_bucket_{curr_bucket}
            , curr_pos_{curr_pos}
            , bucketsList_{blist}
        {
        }

    private:
        std::size_t curr_bucket_ = 0, curr_pos_ = 0;
        BucketList<bucket_size> const& bucketsList_;
    };

public:
    BucketList(std::size_t min_bucket_nbr = 1)
        : buckets_(min_bucket_nbr)
    {
    }

    void add(std::size_t itemIndex);

    void remove(std::size_t itemIndex);

    std::size_t size() const { return bucket_size * (bucket_idx) + curr; }

    auto begin() { return iterator{*this}; }

    auto begin() const { return iterator{*this}; }

    auto end();

    auto end() const;

    BucketListIndex last_idx() const;

    void empty()
    {
        if (capacity() > 3 * size())
        {
            trim(1);
        }
        bucket_idx = 0;
        curr       = 0;
    }


    // TODO
    bool in_bucket(std::size_t itemIndex)
    {
        auto const& index = indexes_[itemIndex];
        if (index.bucket_idx < bucket_idx)
            return index.pos < bucket_size;
        else if (index.bucket_idx == bucket_idx)
            return index.pos < curr;
        else
            return false;
    }

    bool is_empty() const { return bucket_idx == 0 and curr == 0; }

    std::size_t capacity() const { return buckets_.capacity() * bucket_size; }

    void trim(std::size_t max_empty);

    void updateIndex(std::size_t oldIndex, std::size_t newIndex)
    {
        auto bidx                           = indexes_.at(oldIndex);
        buckets_[bidx.bucket_idx][bidx.pos] = newIndex;
        indexes_[newIndex]                  = bidx;
        // should I erase the key oldIndex from indexes_?
    }

private:
    void decrement_()
    {
        if (curr == 0)
        {
            bucket_idx--;
            curr = bucket_size - 1;
        }
        else
        {
            curr--;
        }
    }

    using bucket_t                              = std::array<std::size_t, bucket_size>;
    typename BucketListIndex::Bindex bucket_idx = 0;
    typename BucketListIndex::Bindex curr       = 0;
    std::vector<bucket_t> buckets_;
    std::unordered_map<std::size_t, BucketListIndex> indexes_;
};


template<std::size_t bucket_size>
inline typename BucketList<bucket_size>::iterator BucketList<bucket_size>::iterator::operator++()
{
    curr_pos_++;
    if (curr_pos_ == bucket_size)
    {
        curr_bucket_++;
        curr_pos_ = 0;
    }
    return *this;
}


template<std::size_t bucket_size>
inline bool BucketList<bucket_size>::iterator::operator!=(iterator const& other) const
{
    return (other.curr_bucket_ != curr_bucket_ or other.curr_pos_ != curr_pos_)
           or &bucketsList_ != &other.bucketsList_;
}




template<std::size_t bucket_size>
void BucketList<bucket_size>::add(std::size_t itemIndex)
{
    if (curr != bucket_size)
    {
        indexes_[itemIndex].bucket_idx = bucket_idx;
        indexes_[itemIndex].pos        = curr;
        buckets_[bucket_idx][curr++]   = itemIndex;
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
BucketListIndex BucketList<bucket_size>::last_idx() const
{
    if (curr == 0)
        return {bucket_idx - 1, bucket_size};
    else
        return {bucket_idx, curr - 1};
}




template<std::size_t bucket_size>
void BucketList<bucket_size>::remove(std::size_t itemIndex)
{
    if (!in_bucket(itemIndex))
        throw std::runtime_error("not in bucket");

    auto& to_remove = buckets_[indexes_[itemIndex].bucket_idx][indexes_[itemIndex].pos];
    auto idx        = last_idx();
    auto last       = buckets_[idx.bucket_idx][idx.pos];
    to_remove       = last;
    decrement_();
}



template<std::size_t bucket_size>
auto BucketList<bucket_size>::end()
{
    // if the current cursor position is equal to bucket_size
    // it means we really are positioned on the next
    // bucket at cursor 0
    if (bucket_idx == 0 and curr == 0) // never added
        return iterator{*this, bucket_idx, curr + 1};

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
