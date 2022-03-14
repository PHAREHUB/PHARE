#ifndef PHARE_BUCKETLIST_H
#define PHARE_BUCKETLIST_H

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

namespace PHARE::core
{

template<std::size_t bucket_size>
struct BucketListIndex
{
    using Bindex = int;
    struct BucketIndex
    {
        Bindex val;
    };
    struct CursorIndex
    {
        Bindex val;
    };
    BucketListIndex(Bindex curr_bucket, Bindex curr_pos)
        : pos{curr_pos}
        , bucket_idx{curr_bucket}
    {
    }
    bool operator==(BucketListIndex const& other) const
    {
        return (bucket_idx == other.bucket_idx) and (pos == other.pos);
    }

    bool operator!=(BucketListIndex const& other) const { return !(*this == other); }

    static constexpr Bindex default_idx_v = -1;
    Bindex pos                            = default_idx_v;
    Bindex bucket_idx                     = default_idx_v;
    auto globalCurrPos() const { return pos + bucket_idx * bucket_size; }
};



template<std::size_t bucket_size>
class BucketList
{
    using Bindex = typename BucketListIndex<bucket_size>::Bindex;
    // ----------------- BucketList Iterator ----------------------------------
    template<typename BucketListPtr>
    class bucket_iterator : public std::iterator<std::random_access_iterator_tag, std::size_t>
    {
    public:
        auto& operator*() const { return bucketsList_->buckets_[index_.bucket_idx][index_.pos]; }
        auto& operator*() { return bucketsList_->buckets_[index_.bucket_idx][index_.pos]; }

        auto& operator++();
        auto& operator--();

        auto operator-(bucket_iterator const& other) const;
        auto operator+(int n) const;
        auto operator-(int n) const;

        auto operator<(bucket_iterator const& that) const;
        auto operator>=(bucket_iterator const& that) const;

        auto operator>(bucket_iterator const& that) const;
        auto operator<=(bucket_iterator const& that) const;

        bool operator!=(bucket_iterator const& other) const;
        auto operator==(bucket_iterator const& that) const;

        auto& operator+=(int n);

        bucket_iterator(BucketListPtr blist, Bindex curr_bucket = 0, Bindex curr_pos = 0)
            : index_{curr_bucket, curr_pos}
            , bucketsList_{blist}
        {
        }

        bucket_iterator& operator=(bucket_iterator const& that) = default;
        bucket_iterator& operator=(bucket_iterator&& that) = default;
        bucket_iterator(bucket_iterator const& that)       = default;
        bucket_iterator(bucket_iterator&& that)            = default;

    private:
        BucketListIndex<bucket_size> index_;
        BucketListPtr bucketsList_;
    };
    // ----------------- END BucketList Iterator ------------------------------




public:
    BucketList(std::size_t min_bucket_nbr = 1)
        : index_{0, 0}
        , buckets_(min_bucket_nbr, ConstArray<std::size_t, bucket_size>(
                                       BucketListIndex<bucket_size>::default_idx_v))
    {
    }

    using iterator       = bucket_iterator<BucketList*>;
    using const_iterator = bucket_iterator<BucketList const* const>;

    void add(std::size_t itemIndex);
    void remove(std::size_t itemIndex);
    bool in_bucket(std::size_t itemIndex);

    // reallocate bucketlist memory if more empty space than max_empty
    void trim(std::size_t max_empty);

    // to use if an item in an indexed array is moved at another index
    void updateIndex(std::size_t oldIndex, std::size_t newIndex);

    // empty the bucketlist, but leaves the capacity untouched
    void empty();
    bool is_empty() const { return index_.bucket_idx == 0 and index_.pos == 0; }

    std::size_t size() const { return index_.globalCurrPos(); }
    std::size_t capacity() const { return buckets_.capacity() * bucket_size; }

    auto begin() { return iterator{this}; }
    auto begin() const { return const_iterator{this}; }
    auto cbegin() const { return const_iterator{this}; }
    auto end();
    auto end() const;
    auto cend() const;

    void sort();



private:
    BucketListIndex<bucket_size> lastIdx_() const;


    void decrement_()
    {
        if (index_.pos == 0)
        {
            index_.bucket_idx--;
            index_.pos = bucket_size - 1;
        }
        else
        {
            auto bs    = static_cast<Bindex>(bucket_size);
            auto check = index_.bucket_idx >= 0 and index_.bucket_idx < bs and (index_.pos >= 0)
                         and ((index_.pos - 1) < bs);
            if (!check)
            {
                throw std::runtime_error("BucketListError : bucket_idx ("
                                         + std::to_string(index_.bucket_idx) + "), curr ("
                                         + std::to_string(index_.pos) + " bucket_size("
                                         + std::to_string(bucket_size) + ")");
            }
            buckets_[index_.bucket_idx][index_.pos - 1]
                = BucketListIndex<bucket_size>::default_idx_v;
            index_.pos--;
        }
    }

    using bucket_t = std::array<std::size_t, bucket_size>;
    BucketListIndex<bucket_size> index_;
    std::vector<bucket_t> buckets_;
};


// ----------------- BucketList Iterator ----------------------------------
template<std::size_t bucket_size>
template<typename BucketListPtr>
inline auto& BucketList<bucket_size>::bucket_iterator<BucketListPtr>::operator++()
{
    index_.pos++;
    if (index_.pos == bucket_size)
    {
        index_.bucket_idx++;
        index_.pos = 0;
    }
    return *this;
}

template<std::size_t bucket_size>
template<typename BucketListPtr>
inline auto& BucketList<bucket_size>::bucket_iterator<BucketListPtr>::operator--()
{
    if (index_.pos == 0)
    {
        index_.pos = bucket_size - 1;
        --index_.bucket_idx;
    }
    else
        --index_.pos;
    return *this;
}

template<std::size_t bucket_size>
template<typename BucketListPtr>
inline bool BucketList<bucket_size>::bucket_iterator<BucketListPtr>::operator!=(
    bucket_iterator const& other) const
{
    return (this->index_ != other.index_) or (bucketsList_ != other.bucketsList_);
}

template<std::size_t bucket_size>
template<typename BucketListPtr>
inline auto& BucketList<bucket_size>::bucket_iterator<BucketListPtr>::operator+=(int n)
{
    auto globPos      = index_.globalCurrPos();
    auto newGloPos    = globPos + n;
    index_.bucket_idx = newGloPos / static_cast<int>(bucket_size);
    index_.pos        = newGloPos - index_.bucket_idx * bucket_size;
    return *this;
}

template<std::size_t bucket_size>
template<typename BucketListPtr>
inline auto BucketList<bucket_size>::bucket_iterator<BucketListPtr>::operator+(int n) const
{
    auto copy{*this};
    copy += n;
    return copy;
}


template<std::size_t bucket_size>
template<typename BucketListPtr>
inline auto BucketList<bucket_size>::bucket_iterator<BucketListPtr>::operator-(int n) const
{
    auto copy{*this};
    auto globPos           = copy.index_.globalCurrPos();
    auto newGloPos         = globPos - n;
    copy.index_.bucket_idx = newGloPos / static_cast<int>(bucket_size);
    copy.index_.pos        = newGloPos - copy.index_.bucket_idx * bucket_size;
    assert(copy.index_.pos < static_cast<int>(bucket_size));
    return copy;
}


template<std::size_t bucket_size>
template<typename BucketListPtr>
inline auto BucketList<bucket_size>::bucket_iterator<BucketListPtr>::operator-(
    bucket_iterator const& other) const
{
    return index_.globalCurrPos() - other.index_.globalCurrPos();
}



template<std::size_t bucket_size>
template<typename BucketListPtr>
inline auto BucketList<bucket_size>::bucket_iterator<BucketListPtr>::operator<(
    bucket_iterator const& that) const
{
    return index_.globalCurrPos() < that.index_.globalCurrPos();
}


template<std::size_t bucket_size>
template<typename BucketListPtr>
inline auto BucketList<bucket_size>::bucket_iterator<BucketListPtr>::operator>=(
    bucket_iterator const& that) const
{
    return !(*this < that);
}

template<std::size_t bucket_size>
template<typename BucketListPtr>
inline auto BucketList<bucket_size>::bucket_iterator<BucketListPtr>::operator>(
    bucket_iterator const& that) const
{
    return index_.globalCurrPos() > that.index_.globalCurrPos();
}


template<std::size_t bucket_size>
template<typename BucketListPtr>
inline auto BucketList<bucket_size>::bucket_iterator<BucketListPtr>::operator<=(
    bucket_iterator const& that) const
{
    return !(*this > that);
}

template<std::size_t bucket_size>
template<typename BucketListPtr>
inline auto BucketList<bucket_size>::bucket_iterator<BucketListPtr>::operator==(
    bucket_iterator const& that) const
{
    return index_ == that.index_ and bucketsList_ == that.bucketsList_;
}

// ----------------- BucketList Iterator ----------------------------------



template<std::size_t bucket_size>
void BucketList<bucket_size>::add(std::size_t itemIndex)
{
    if (index_.pos != bucket_size)
    {
        buckets_[index_.bucket_idx][index_.pos++] = itemIndex;
    }
    else
    {
        if (size() == capacity())
        {
            buckets_.reserve(buckets_.size() + 1);
            buckets_.push_back(bucket_t{});
        }
        index_.bucket_idx++;
        index_.pos = 0;
        add(itemIndex);
    }
}

template<std::size_t bucket_size>
BucketListIndex<bucket_size> BucketList<bucket_size>::lastIdx_() const
{
    if (index_.pos == 0)
        return {index_.bucket_idx - 1, bucket_size - 1};
    else
    {
        return {index_.bucket_idx, index_.pos - 1};
    }
}


template<std::size_t bucket_size>
void BucketList<bucket_size>::empty()
{
    constexpr std::size_t arbitraryTrimThresh = 3;
    if (capacity() > arbitraryTrimThresh * size())
    {
        trim(1);
    }
    index_.bucket_idx = 0;
    index_.pos        = 0;
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
    if (index_.bucket_idx == 0 and index_.pos == 0) // never added
        return iterator{this, index_.bucket_idx, index_.pos};

    else if (index_.pos != bucket_size)
    {
        auto it = iterator{this, index_.bucket_idx, index_.pos};
        return it;
    }
    else
    {
        auto it = iterator{this, index_.bucket_idx + 1, 0};
        return it;
    }
}

template<std::size_t bucket_size>
auto BucketList<bucket_size>::end() const
{
    // if the current cursor position is equal to bucket_size
    // it means we really are positioned on the next
    // bucket at cursor 0
    if (index_.pos != bucket_size)
    {
        auto it = const_iterator{this, index_.bucket_idx, index_.pos};
        return it;
    }
    else
    {
        auto it = const_iterator{this, index_.bucket_idx + 1, 0};
        return it;
    }
}

template<std::size_t bucket_size>
auto BucketList<bucket_size>::cend() const
{
    // if the current cursor position is equal to bucket_size
    // it means we really are positioned on the next
    // bucket at cursor 0
    if (index_.pos != bucket_size)
    {
        auto it = const_iterator{this, index_.bucket_idx, index_.pos};
        return it;
    }
    else
    {
        auto it = const_iterator{this, index_.bucket_idx + 1, 0};
        return it;
    }
}

template<std::size_t bucket_size>
void BucketList<bucket_size>::trim(std::size_t max_empty)
{
    auto nbr_elem         = size();
    auto nbr_full_buckets = nbr_elem / bucket_size;
    auto nbr_buckets      = buckets_.size();

    auto occupied_buckets
        = nbr_full_buckets + ((index_.pos == 0 or index_.pos == bucket_size) ? 0 : 1);
    auto curr_empty = nbr_buckets - occupied_buckets;
    if (curr_empty > max_empty)
    {
        std::vector<bucket_t> new_buckets_(std::begin(buckets_),
                                           std::begin(buckets_) + occupied_buckets + max_empty);
        buckets_.clear();
        buckets_.swap(new_buckets_);
    }
}



template<std::size_t bucket_size>
void BucketList<bucket_size>::sort()
{
    std::sort(std::begin(*this), std::end(*this));
}



} // namespace PHARE::core
#endif
