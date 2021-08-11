// https://stackoverflow.com/a/51311730/795574

#ifndef PHARE_CORE_RANGE_HPP
#define PHARE_CORE_RANGE_HPP


#include <string>
#include <sstream>
#include <iostream>
#include <chrono>
#include <vector>
#include <utility>
#include <memory>


namespace PHARE
{
template<class F>
struct generator_iterator
{
    F f;
    std::size_t i = 0;

    using self = generator_iterator;
    friend bool operator==(self const& lhs, self const& rhs) { return lhs.i == rhs.i; }
    friend bool operator!=(self const& lhs, self const& rhs) { return lhs.i != rhs.i; }
    using reference         = std::result_of_t<F const&(std::size_t const&)>;
    using value_type        = std::decay_t<reference>;
    using difference_type   = std::ptrdiff_t;
    using pointer           = value_type*;
    using iterator_category = std::input_iterator_tag;

    self& operator++()
    {
        ++i;
        return *this;
    }
    self operator++(int)
    {
        auto tmp = *this;
        ++*this;
        return tmp;
    }
    reference operator*() const { return f(i); }
    pointer operator->() const { return std::addressof(f(i)); }

    friend difference_type operator-(self const& lhs, self const& rhs) { return lhs.i - rhs.i; }
    friend difference_type operator-(self const& lhs, difference_type rhs) { return lhs.i - rhs; }
    friend difference_type operator+(self const& lhs, difference_type rhs) { return lhs.i + rhs; }
};

template<class It>
struct range_t
{
    It b, e;
    It begin() const { return b; }
    It end() const { return e; }

private:
    struct container_maker
    {
        range_t const* self;

        template<class C>
        operator C() &&
        {
            return {self->begin(), self->end()};
        }
    };

public:
    container_maker make_container() const { return {this}; }
    // C is optional
    template<class C>
    C make_container() const
    {
        return make_container();
    }

    std::vector<typename It::value_type> operator()() { return make_container(); }
};

template<class It>
range_t<It> range(It s, It f)
{
    return {std::move(s), std::move(f)};
}

template<class F>
auto generate(F&& f, std::size_t count)
{
    generator_iterator<std::decay_t<F>> e{f, count};
    generator_iterator<std::decay_t<F>> b{std::forward<F>(f)};
    return range(std::move(b), std::move(e));
}

} // namespace PHARE

#endif /*PHARE_CORE_RANGE_HPP*/
