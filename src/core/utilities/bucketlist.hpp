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
#include <tuple>

namespace PHARE::core
{
using Bindex = int;

class BucketIndex
{
public:
    explicit BucketIndex(Bindex val)
        : val_{val}
    {
    }

    auto& get() { return val_; }
    auto const& get() const { return val_; }

private:
    Bindex val_;
};

class CursorIndex
{
public:
    explicit CursorIndex(Bindex val)
        : val_{val}
    {
    }

    auto& get() { return val_; }
    auto const& get() const { return val_; }

private:
    Bindex val_;
};

template<std::size_t bucket_size>
class BucketListIndex
{
public:
    BucketListIndex(BucketIndex curr_bucket, CursorIndex curr_pos)
        : pos_{curr_pos}
        , bucket_idx_{curr_bucket}
    {
    }

    bool operator==(BucketListIndex const& other) const
    {
        return (bucket() == other.bucket()) and (cursor() == other.cursor());
    }

    bool operator!=(BucketListIndex const& other) const { return !(*this == other); }


    auto& operator++()
    {
        cursor()++;
        if (cursor() == bucket_size)
        {
            bucket()++;
            cursor() = 0;
        }
        return *this;
    }


    auto& operator--()
    {
        if (cursor() == 0)
        {
            cursor() = bucket_size - 1;
            --bucket();
        }
        else
            --cursor();
        return *this;
    }


    auto operator+=(int n)
    {
        auto globPos   = globalCurrPos();
        auto newGloPos = globPos + n;
        bucket()       = newGloPos / static_cast<int>(bucket_size);
        cursor()       = newGloPos - bucket() * bucket_size;
    }


    auto operator-(int n) const
    {
        auto copy{*this};
        auto globPos   = copy.globalCurrPos();
        auto newGloPos = globPos - n;
        copy.bucket()  = newGloPos / static_cast<int>(bucket_size);
        copy.cursor()  = newGloPos - copy.bucket() * bucket_size;
        assert(copy.cursor() < static_cast<int>(bucket_size));
        return copy;
    }


    auto operator-(BucketListIndex const& other) const
    {
        return globalCurrPos() - other.globalCurrPos();
    }


    auto operator<(BucketListIndex const& other) const
    {
        return globalCurrPos() < other.globalCurrPos();
    }

    auto operator>(BucketListIndex const& other) const
    {
        return globalCurrPos() > other.globalCurrPos();
    }

    auto& bucket() { return bucket_idx_.get(); }
    auto const& bucket() const { return bucket_idx_.get(); }
    auto& cursor() { return pos_.get(); }
    auto const& cursor() const { return pos_.get(); }

    auto isBucketEnd() const { return cursor() == bucket_size; }

    void reset()
    {
        cursor() = 0;
        bucket() = 0;
    }
    auto globalCurrPos() const { return cursor() + bucket() * bucket_size; }
    static constexpr Bindex default_idx_v = -1;

    template<std::size_t s>
    friend auto& operator<<(std::ostream& os, BucketListIndex<s> const& index);


private:
    CursorIndex pos_{default_idx_v};
    BucketIndex bucket_idx_{default_idx_v};
};

template<std::size_t bucket_size>
auto& operator<<(std::ostream& os, BucketListIndex<bucket_size> const& index)
{
    os << "bucket : " << index.bucket() << " cursor " << index.cursor();
    return os;
}


template<std::size_t bucket_size>
class Buckets
{
private:
    using bucket_t = std::array<std::size_t, bucket_size>;

public:
    Buckets(std::size_t min_bucket_nbr)
        : buckets_(min_bucket_nbr, ConstArray<std::size_t, bucket_size>(
                                       BucketListIndex<bucket_size>::default_idx_v))
    {
    }

    auto& operator()(BucketListIndex<bucket_size> const& index)
    {
        assert(index.cursor() < static_cast<int>(bucket_size));
        return buckets_[index.bucket()][index.cursor()];
    }

    auto const& operator()(BucketListIndex<bucket_size> const& index) const
    {
        return buckets_[index.bucket()][index.cursor()];
    }

    auto& operator[](std::size_t i) { return buckets_[i]; }
    auto const& operator[](std::size_t i) const { return buckets_[i]; }

    void push_back(bucket_t&& b) { buckets_.push_back(b); }
    void push_back(bucket_t const& b) { buckets_.push_back(b); }
    void reserve(std::size_t size) { buckets_.reserve(size); }
    auto size() const { return buckets_.size(); }
    void clear() { buckets_.clear(); }
    void swap(std::vector<bucket_t>& other) { buckets_.swap(other); }
    auto capacity() const { return buckets_.capacity(); }
    auto begin() { return buckets_.begin(); }
    auto end() { return buckets_.end(); }
    void add_bucket()
    {
        buckets_.reserve(buckets_.size() + 1);
        buckets_.push_back(bucket_t{});
    }
    bool isFull() const { return size() == capacity(); }

private:
    std::vector<bucket_t> buckets_;
};

template<std::size_t bucket_size>
class BucketList
{
    // ----------------- BucketList Iterator ----------------------------------
    template<typename BucketListPtr>
    class bucket_iterator : public std::iterator<std::random_access_iterator_tag, std::size_t>
    {
    public:
        auto const& operator*() const { return bucketsList_->buckets_(index_); }
        auto& operator*() { return bucketsList_->buckets_(index_); }

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
            : index_{BucketIndex{curr_bucket}, CursorIndex{curr_pos}}
            , bucketsList_{blist}
        {
        }

        bucket_iterator(BucketListPtr blist, BucketListIndex<bucket_size> index)
            : index_{index}
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
        : index_{BucketIndex{0}, CursorIndex{0}}
        , buckets_(min_bucket_nbr)
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
    bool is_empty() const { return index_.bucket() == 0 and index_.cursor() == 0; }

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
        if (index_.cursor() == 0)
        {
            index_.bucket()--;
            index_.cursor() = bucket_size - 1;
        }
        else
        {
            auto bs    = static_cast<Bindex>(bucket_size);
            auto check = index_.bucket() >= 0 and index_.bucket() < bs and (index_.cursor() >= 0)
                         and ((index_.cursor() - 1) < bs);
            if (!check)
            {
                throw std::runtime_error("BucketListError : bucket_idx ("
                                         + std::to_string(index_.bucket()) + "), curr ("
                                         + std::to_string(index_.cursor()) + " bucket_size("
                                         + std::to_string(bucket_size) + ")");
            }
            buckets_(index_) = BucketListIndex<bucket_size>::default_idx_v;
            index_.cursor()--;
        }
    }

    using bucket_t = std::array<std::size_t, bucket_size>;
    BucketListIndex<bucket_size> index_;
    Buckets<bucket_size> buckets_;
};


// ----------------- BucketList Iterator ----------------------------------
template<std::size_t bucket_size>
template<typename BucketListPtr>
inline auto& BucketList<bucket_size>::bucket_iterator<BucketListPtr>::operator++()
{
    ++index_;
    return *this;
}

template<std::size_t bucket_size>
template<typename BucketListPtr>
inline auto& BucketList<bucket_size>::bucket_iterator<BucketListPtr>::operator--()
{
    --index_;
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
    index_ += n;
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
    return bucket_iterator{bucketsList_, BucketListIndex<bucket_size>{this->index_ - n}};
}


template<std::size_t bucket_size>
template<typename BucketListPtr>
inline auto BucketList<bucket_size>::bucket_iterator<BucketListPtr>::operator-(
    bucket_iterator const& other) const
{
    return index_ - other.index_;
}



template<std::size_t bucket_size>
template<typename BucketListPtr>
inline auto BucketList<bucket_size>::bucket_iterator<BucketListPtr>::operator<(
    bucket_iterator const& that) const
{
    return index_ < that.index_;
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
    return index_ > that.index_;
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
    if (size() == capacity())
        buckets_.add_bucket();
    buckets_(index_) = itemIndex;
    ++index_;
}


template<std::size_t bucket_size>
BucketListIndex<bucket_size> BucketList<bucket_size>::lastIdx_() const
{
    if (index_.cursor() == 0)
        return {BucketIndex{index_.bucket() - 1}, CursorIndex{bucket_size - 1}};
    else
    {
        return {BucketIndex{index_.bucket()}, CursorIndex{index_.cursor() - 1}};
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
    index_.reset();
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
        *it = buckets_(lastIdx_());
        decrement_();
    }
}



template<std::size_t bucket_size>
auto BucketList<bucket_size>::end()
{
    if (is_empty())
        return iterator{this, index_};

    else if (index_.cursor() != bucket_size)
    {
        auto it = iterator{this, index_};
        return it;
    }
    else
    {
        // if the current cursor position is equal to bucket_size
        // it means we really are positioned on the next
        // bucket at cursor 0
        auto it = iterator{this, index_.bucket() + 1, 0};
        return it;
    }
}

template<std::size_t bucket_size>
auto BucketList<bucket_size>::end() const
{
    // if the current cursor position is equal to bucket_size
    // it means we really are positioned on the next
    // bucket at cursor 0
    if (!index_.isBucketEnd())
    {
        auto it = const_iterator{this, index_};
        return it;
    }
    else
    {
        auto it = const_iterator{this, index_.bucket() + 1, 0};
        return it;
    }
}

template<std::size_t bucket_size>
auto BucketList<bucket_size>::cend() const
{
    // if the current cursor position is equal to bucket_size
    // it means we really are positioned on the next
    // bucket at cursor 0
    if (index_.cursor() != bucket_size)
    {
        auto it = const_iterator{this, index_};
        return it;
    }
    else
    {
        auto it = const_iterator{this, index_.bucket() + 1, 0};
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
        = nbr_full_buckets + ((index_.cursor() == 0 or index_.cursor() == bucket_size) ? 0 : 1);
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
