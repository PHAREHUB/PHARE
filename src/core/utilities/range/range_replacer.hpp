#ifndef PHARE_CORE_UTILITITES_RANGE_REPLACER_HPP
#define PHARE_CORE_UTILITITES_RANGE_REPLACER_HPP


#include <vector>
#include <iterator>
#include <cstddef>
#include <iterator>
#include <algorithm>


#include "core/utilities/range/ranges.hpp"


namespace PHARE::core
{
template<typename Container, typename RangeType>
struct SyncInfo
{
    SyncInfo(Container& container_, RangeType&& range_)
        : container{container_}
        , range{range_}
    {
    }

    Container& container;
    RangeType range;
    std::size_t offset{};
};

template<typename Container>
class RangeSynchrotron
{
    using iterator         = typename Container::iterator;
    using value_type       = typename Container::value_type;
    using RangeSyncInfo_t  = SyncInfo<Container, Range<iterator>>;
    using VectorSyncInfo_t = SyncInfo<Container, Container>;
    using RangeSyncVec_t   = std::vector<std::unique_ptr<RangeSyncInfo_t>>;
    using VectorSyncVec_t  = std::vector<std::unique_ptr<VectorSyncInfo_t>>;

public:
    ~RangeSynchrotron() { sync(); }

    void sync()
    {
        _erase(to_erase_list_per_thread);
        _add(to_add_list_per_thread);
        _add(to_add_copy_list_per_thread);
        _erase(to_erase_src_list_per_thread);
    }

    RangeSynchrotron(std::uint16_t n_threads_ = 1)
        : n_threads{n_threads_}
        , to_add_list_per_thread(n_threads)
        , to_add_copy_list_per_thread(n_threads)
        , to_erase_list_per_thread(n_threads)
        , to_erase_src_list_per_thread(n_threads)
    {
    }

    void add(std::uint16_t const thread_idx, Container& container, Range<iterator>&& range,
             bool copy = false)
    {
        if (copy)
            to_add_copy_list_per_thread[thread_idx].emplace_back(std::make_unique<VectorSyncInfo_t>(
                container, Container{range.begin(), range.end()}));
        else
            to_add_list_per_thread[thread_idx].emplace_back(
                std::make_unique<RangeSyncInfo_t>(container, std::forward<Range<iterator>>(range)));
    }

    void erase(std::uint16_t const thread_idx, Container& container, std::size_t const offset,
               Range<iterator>&& range, bool src = false)
    {
        auto& list_per_thread = src ? to_erase_src_list_per_thread : to_erase_list_per_thread;

        list_per_thread[thread_idx]
            .emplace_back(
                std::make_unique<RangeSyncInfo_t>(container, std::forward<Range<iterator>>(range)))
            ->offset
            = offset;
    }

    std::uint16_t const n_threads = 1;

private:
    template<typename AddType>
    void static _add(AddType& list_per_thread)
    {
        for (auto& to_add_list : list_per_thread)
        {
            for (auto& to_add : to_add_list)
                std::copy(to_add->range.begin(), to_add->range.end(),
                          std::back_inserter(to_add->container));

            to_add_list.clear();
        }
    }

    // erases must be in reverse, i.e. container.end() to container.begin()
    void static _erase(std::vector<RangeSyncVec_t>& list_per_thread)
    { // assumption is higher threads have higher offsets
        for (auto it = list_per_thread.rbegin(); it < list_per_thread.rend(); ++it)
        {
            auto& to_erase_list = *it;

            std::sort(to_erase_list.begin(), to_erase_list.end(), // highest offset first
                      [](auto const& a, auto const& b) { return a->offset > b->offset; });

            for (auto& to_erase : to_erase_list)
                to_erase->container.erase(to_erase->container.begin() + to_erase->offset,
                                          to_erase->container.begin() + to_erase->offset
                                              + to_erase->range.size());

            to_erase_list.clear();
        }
    }

    std::vector<RangeSyncVec_t> to_add_list_per_thread;
    std::vector<VectorSyncVec_t> to_add_copy_list_per_thread;
    std::vector<RangeSyncVec_t> to_erase_list_per_thread;
    std::vector<RangeSyncVec_t> to_erase_src_list_per_thread;
};

template<typename Container>
class RangeReplacer
{
public:
    using RangeSynchrotron_t = RangeSynchrotron<Container>;
    using iterator           = typename Container::iterator;
    using Range_t            = Range<iterator>;
    using PartitionInfo      = std::vector<std::pair<std::size_t, std::size_t>>;

    RangeReplacer(Container& dst_, RangeSynchrotron_t& synchrotron,
                  std::uint16_t const thread_idx_ = 0)
        : dst{dst_}
        , synchrotron_{synchrotron}
        , thread_idx{thread_idx_}
    {
    }

    auto static make_unique(Container& dst, RangeSynchrotron_t& synchrotron,
                            std::uint16_t const thread_idx = 0)
    {
        return std::make_unique<RangeReplacer<Container>>(dst, synchrotron, thread_idx);
    }

    void replace(Container& src, Range_t src_range, bool copy_src = false, bool erase_src = false)
    {
        auto src_size  = src_range.size();
        auto src_begin = src_range.begin();

        for (auto& partition_pair : replaceable)
        {
            if (src_size == 0)
                break;

            auto& [offset, n_replaceable] = partition_pair;

            if (n_replaceable > 0)
            {
                auto replaceable = src_size > n_replaceable ? n_replaceable : src_size;

                std::copy(src_begin, src_begin + replaceable, dst.begin() + offset);
                src_size -= replaceable;
                n_replaceable -= replaceable;
                src_begin += replaceable;
                offset += replaceable;
            }
        }

        if (src_size > 0)
        {
            synchrotron_.add(thread_idx, dst, makeRange(src_begin, src_range.end()), copy_src);
        }

        if (erase_src)
        {
            std::size_t offset = std::distance(src.begin(), src_range.begin());
            synchrotron_.erase(thread_idx, src, offset,
                               makeRange(src_range.begin(), src_range.end()), true);
        }
    }


    void erase()
    {
        std::sort(replaceable.begin(), replaceable.end(),
                  [](auto const& a, auto const& b) { return a > b; });

        for (auto& [offset, n_replaceable] : replaceable)
        {
            if (n_replaceable > 0)
            {
                if (threads_ == 1)
                    dst.erase(dst.begin() + offset, dst.begin() + offset + n_replaceable);
                else
                    synchrotron_.erase(
                        thread_idx, dst, offset,
                        makeRange(dst.begin() + offset, dst.begin() + offset + n_replaceable));
            }
        }
        replaceable.clear();
    }

    void add_replaceable(Range_t range, iterator newEnd)
    {
        auto offset = std::distance(dst.begin(), newEnd);
        auto size   = std::distance(newEnd, range.end());
        replaceable.push_back({offset, size});
    }

    template<typename Fn>
    void add_replaceable(Range_t range, Fn fn)
    {
        add_replaceable(range, std::partition(range.begin(), range.end(), fn));
    }

    void replace_with(Container& src, Range_t range, iterator newEnd)
    {
        replace(src, makeRange(newEnd, range.end()));
    }

    template<typename Fn>
    void replace_with(Container& src, Range_t range, Fn fn)
    {
        replace_with(src, range, std::partition(range.begin(), range.end(), fn));
    }

private:
    Container& dst;
    RangeSynchrotron_t& synchrotron_;
    std::uint16_t const thread_idx = -1;
    PartitionInfo replaceable{};
    std::uint16_t threads_ = synchrotron_.n_threads;
};
} // namespace PHARE::core

#endif // PHARE_CORE_UTILITITES_RANGE_REPLACER_HPP
