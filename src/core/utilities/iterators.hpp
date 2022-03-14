
#ifndef PHARE_CORE_UTILITIES_ITERATORS_HPP
#define PHARE_CORE_UTILITIES_ITERATORS_HPP

namespace PHARE::core
{
template<typename Vector, bool is_const>
using wrapped_iterator_base
    = std::conditional_t<is_const, typename Vector::const_iterator, typename Vector::iterator>;


template<typename T, typename Vector, bool is_const = std::is_const_v<T>>
struct wrapped_iterator : public wrapped_iterator_base<Vector, is_const>
{
    using outer_type = std::decay_t<T>;
    using value_type = typename Vector::value_type;
    using Super      = wrapped_iterator_base<Vector, is_const>;

    wrapped_iterator operator+(std::int64_t i)
    {
        wrapped_iterator copy = *this;
        static_cast<Super&>(copy) += i;
        return copy;
    }

    auto& operator()() { return *container; }
    auto& operator()() const { return *container; }


    auto idx() const
    {
        auto val = std::distance(static_cast<Super const&>(container->begin()),
                                 static_cast<Super const&>(*this));
        if (val < 0)
            throw std::runtime_error("Invalid iterator distance");
        return static_cast<std::size_t>(val);
    }


    T* container = nullptr; // might be const
};

} // namespace PHARE::core



#endif /* PHARE_CORE_UTILITIES_ITERATORS_HPP */
