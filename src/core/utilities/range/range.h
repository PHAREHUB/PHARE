#ifndef RANGE_H
#define RANGE_H

#include <iterator>

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
            , last_{std::end(c)}
        {
        }

        Range(Iterator begin, Iterator end)
            : first_{begin}
            , last_{end}
        {
        }


        Iterator begin() const { return first_; }

        Iterator end() const { return last_; }

        std::size_t size() const { return std::distance(first_, last_); }

        auto& operator[](std::size_t i) { return *(first_ + i); }
        auto& operator[](std::size_t i) const { return *(first_ + i); }

    private:
        Iterator first_;
        Iterator last_;
    };

    template<typename Iterator>
    Range<Iterator> makeRange(Iterator&& begin, Iterator&& last)
    {
        return Range<Iterator>{std::forward<Iterator>(begin), std::forward<Iterator>(last)};
    }

    template<typename Container>
    auto makeRange(Container& container)
    {
        return makeRange(std::begin(container), std::end(container));
    }

} // namespace core
} // namespace PHARE

#endif // RANGE_H
