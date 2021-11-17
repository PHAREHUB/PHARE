
#ifndef PHARE_CORE_UTILITIES_ITERATORS_H
#define PHARE_CORE_UTILITIES_ITERATORS_H

#include <type_traits>

namespace PHARE::core
{
template<typename T, typename Vector>
struct wrapped_iterator
{
    auto constexpr static is_const = std::is_const_v<T>;
    using value_type               = typename Vector::value_type;
    using outer_type               = T;
    using iterator_t
        = std::conditional_t<is_const, typename Vector::const_iterator, typename Vector::iterator>;

    using iterator_category = typename iterator_t::iterator_category;
    using difference_type   = typename iterator_t::difference_type;
    using reference         = typename iterator_t::reference;
    using pointer           = typename iterator_t::pointer;

    auto operator-(wrapped_iterator const& other) { return it - other.it; }

    auto& operator+=(std::size_t i)
    {
        it += i;
        return *this;
    }


    auto& operator-=(wrapped_iterator const& that)
    {
        it -= that;
        return *this;
    }


    auto& operator++()
    {
        ++it;
        return *this;
    }

    auto& operator--()
    {
        --it;
        return *this;
    }

    bool operator!=(wrapped_iterator const& that) const { return this->it != that.it; }
    wrapped_iterator& operator=(wrapped_iterator const& src) = default;
    bool operator==(wrapped_iterator const& src) const
    {
        return (container == src.container and it == src.it);
    }

    auto operator->() { return &(*it); }
    auto operator->() const { return &(*it); }

    auto& operator*() { return *it; }
    auto& operator*() const { return *it; }

    auto& operator()() const { return *container; }
    auto& operator()() { return *container; }


    iterator_t it;

    T* container; // might be const
};

} // namespace PHARE::core


namespace std
{
using namespace PHARE::core;

template<typename T, typename V>
auto distance(wrapped_iterator<T, V> const& a, wrapped_iterator<T, V> const& b)
{
    return std::distance(a.it, b.it);
}

} // namespace std


#endif /* PHARE_CORE_UTILITIES_ITERATORS_H */
