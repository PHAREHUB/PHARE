#ifndef RANGE_HPP
#define RANGE_HPP

#include <iterator>
#include <type_traits>

namespace PHARE
{
namespace core
{
    template<typename Iterator>
    struct Range
    {
        using iterator_category = typename Iterator::iterator_category;
        using value_type        = typename Iterator::value_type;
        using difference_type   = typename Iterator::difference_type;
        using reference         = typename Iterator::reference;
        using pointer           = typename Iterator::pointer;

        template<class Container>
        explicit Range(Container const& c)
            : first_{std::begin(c)}
            , end_{std::end(c)}
        {
        }

        Range(Iterator begin, Iterator end)
            : first_{begin}
            , end_{end}
        {
        }
        Iterator begin() const { return first_; }
        Iterator end() const { return end_; }
        std::size_t size() const { return std::distance(first_, end_); }

    private:
        Iterator first_;
        Iterator end_;
    };

    template<typename Array, typename Index = std::size_t,
             typename = std::enable_if_t<std::is_integral_v<Index>>>
    class IndexRange
    {
    public:
        using array_t = Array;

        IndexRange(Array& arr, Index first, Index end)
            : first_{first}
            , end_{end}
            , array_{&arr}
        {
        }

        IndexRange(Array& arr)
            : first_{0}
            , end_{arr.size()}
            , array_{&arr}
        {
        }

        IndexRange& operator=(IndexRange const& from) = default;
        IndexRange& operator=(IndexRange&& from) = default;
        IndexRange(IndexRange&& from)            = default;
        IndexRange(IndexRange const& from)       = default;

        auto size() const { return end_ - first_; }
        auto ibegin() const { return first_; }
        auto iend() const { return end_; }
        auto begin() const { return std::begin(*array_) + first_; }
        auto end() const { return std::begin(*array_) + end_; }
        auto& array() { return *array_; }
        auto const& array() const { return *array_; }

        // these ones are not used for now... they give access to the array element
        // via a range index (i.e. not an array index, but relative to the array index
        // at the start of the range).
        // it is commented out since it may be used now by mistake by a function
        // that'd take an Array or a Range as template arg... resulting in wrong results
        // more thinking is needed if they become useful...
        // auto& operator[](std::size_t idx) { return (*array_)[first_ + idx]; }
        // auto const& operator[](std::size_t idx) const { return (*array_)[first_ + idx]; }


    private:
        Index first_;
        Index end_;
        Array* array_;
    };

    template<typename Array>
    auto makeRange(Array& arr, std::size_t begin, std::size_t end)
    {
        return IndexRange{arr, begin, end};
    }

    template<typename Iterator>
    Range<Iterator> makeRange(Iterator&& begin, Iterator&& end)
    {
        return Range{std::forward<Iterator>(begin), std::forward<Iterator>(end)};
    }


    template<typename Array, typename Index>
    auto makeRange(IndexRange<Array, Index> irange)
    {
        auto& arr  = irange.array();
        auto begin = std::begin(arr) + irange.ibegin();
        auto end   = std::begin(arr) + irange.iend();
        return makeRange(begin, end);
    }

    template<typename Container>
    auto makeIndexRange(Container& container)
    {
        return IndexRange<Container>{container, 0, container.size()};
    }

    template<typename Container>
    auto makeRange(Container& container)
    {
        return makeRange(std::begin(container), std::end(container));
    }

} // namespace core
} // namespace PHARE

#endif // RANGE_HPP
