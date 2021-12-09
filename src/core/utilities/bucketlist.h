#ifndef PHARE_BUCKETLIST_H
#define PHARE_BUCKETLIST_H

#include "core/utilities/meta/meta_utilities.h"
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
    std::uint_fast16_t bucket_idx = 0;
    std::uint_fast16_t pos        = 0;
};


struct BucketListItem
{
    BucketListIndex index;
};




template<std::size_t bucket_size, typename T>
class BucketList
{
    class iterator : public std::iterator<std::forward_iterator_tag, T>
    {
    public:
        auto operator*() const { return bucketsList_.buckets_[curr_bucket_][curr_pos_]; }
        auto operator*() { return bucketsList_.buckets_[curr_bucket_][curr_pos_]; }
        iterator operator++();
        bool operator!=(iterator const& other) const;

        iterator(BucketList const& blist, std::size_t curr_bucket = 0, std::size_t curr_pos = 0)
            : curr_bucket_{curr_bucket}
            , curr_pos_{curr_pos}
            , bucketsList_{blist}
        {
        }

    private:
        std::size_t curr_bucket_ = 0, curr_pos_ = 0;
        BucketList<bucket_size, T> const& bucketsList_;
    };

public:
    BucketList(std::size_t min_bucket_nbr = 1)
        : buckets_(min_bucket_nbr)
    {
    }

    void add(T& t);

    void remove(T& t);

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


    bool in_bucket(BucketListItem const& t)
    {
        auto const& index = t.index;
        return index.bucket_idx <= bucket_idx and index.pos <= curr;
    }

    bool is_empty() const { return bucket_idx == 0 and curr == 0; }

    std::size_t capacity() const { return buckets_.capacity() * bucket_size; }

    void trim(std::size_t max_empty);




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

    using bucket_t         = std::array<const T*, bucket_size>;
    std::size_t bucket_idx = 0;
    std::size_t curr       = 0;
    std::vector<bucket_t> buckets_;
};


template<std::size_t bucket_size, typename T>
inline typename BucketList<bucket_size, T>::iterator
BucketList<bucket_size, T>::iterator::operator++()
{
    curr_pos_++;
    if (curr_pos_ == bucket_size)
    {
        curr_bucket_++;
        curr_pos_ = 0;
    }
    return *this;
}


template<std::size_t bucket_size, typename T>
inline bool BucketList<bucket_size, T>::iterator::operator!=(iterator const& other) const
{
    return (other.curr_bucket_ != curr_bucket_ or other.curr_pos_ != curr_pos_)
           or &bucketsList_ != &other.bucketsList_;
}


template<typename BucketListItem, typename Attempt = void>
struct has_index : std::false_type
{
};

template<typename BucketListItem>
struct has_index<BucketListItem,
                 core::tryToInstanciate<decltype(std::declval<BucketListItem>().index)>>
    : std::true_type
{
};


template<typename BucketListItem>
static auto constexpr inline has_index_v = has_index<BucketListItem>::value;


template<std::size_t bucket_size, typename T>
void BucketList<bucket_size, T>::add(T& t)
{
    if (curr != bucket_size)
    {
        if constexpr (has_index_v<T>)
        {
            t.index.bucket_idx = bucket_idx;
            t.index.pos        = curr;
        }
        buckets_[bucket_idx][curr++] = &t;
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
        add(t);
    }
}

template<std::size_t bucket_size, typename T>
BucketListIndex BucketList<bucket_size, T>::last_idx() const
{
    if (curr == 0)
        return {bucket_idx - 1, bucket_size};
    else
        return {bucket_idx, curr - 1};
}




template<std::size_t bucket_size, typename T>
void BucketList<bucket_size, T>::remove(T& t)
{
    if constexpr (has_index_v<T>)
    {
        if (!in_bucket(t))
            throw std::runtime_error("not in bucket " + std::to_string(t.index.bucket_idx) + " "
                                     + std::to_string(t.index.pos));

        auto& to_remove = buckets_[t.index.bucket_idx][t.index.pos];
        auto idx        = last_idx();
        auto last       = buckets_[idx.bucket_idx][idx.pos];
        to_remove       = last;
        decrement_();
    }
    else
        throw std::runtime_error("cannot remove this type");
}



template<std::size_t bucket_size, typename T>
auto BucketList<bucket_size, T>::end()
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


template<std::size_t bucket_size, typename T>
auto BucketList<bucket_size, T>::end() const
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

template<std::size_t bucket_size, typename T>
void BucketList<bucket_size, T>::trim(std::size_t max_empty)
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
