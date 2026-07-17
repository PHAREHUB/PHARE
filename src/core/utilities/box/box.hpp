#ifndef PHARE_CORE_UTILITIES_BOX_BOX_HPP
#define PHARE_CORE_UTILITIES_BOX_BOX_HPP


#include "core/def.hpp"
#include "core/logger.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/span.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/meta/meta_utilities.hpp"

#include <tuple>
#include <limits>
#include <vector>
#include <cstddef>
#include <iostream>
#include <optional>
#include <algorithm>
#include <unordered_map>


namespace PHARE::core
{
template<typename Type, std::size_t dim>
class box_iterator;

namespace
{
    template<std::size_t dim, typename Point>
    void constexpr verify(Point const& lower, Point const& upper)
    {
        for (std::uint16_t i = 0; i < dim; ++i)
            if (lower[i] > upper[i])
            {
                PHARE_LOG_LINE_STR("Invalid box " + std::to_string(lower) + "-"
                                   + std::to_string(upper));
                std::abort();
            }
    }
} // namespace

/** Represents a 1D, 2D or 3D box of integer or floating point
 * points.
 */
template<typename Type, std::size_t dim>
struct Box
{
    using value_type                = Type;
    using iterator                  = box_iterator<Type, dim>;
    static constexpr auto dimension = dim;

    using Point_t = Point<Type, dim>;


    Point_t lower{ConstArray<Type, dim>()};
    Point_t upper{ConstArray<Type, dim>()};


    Box()                      = default;
    Box(Box const&)            = default;
    Box(Box&&)                 = default;
    Box& operator=(Box const&) = default;
    Box& operator=(Box&&)      = default;


    constexpr Box(std::array<Type, dim> const& _lower, std::array<Type, dim> const& _upper)
        : lower{_lower}
        , upper{_upper}
    {
        PHARE_DEBUG_DO(verify<dim>(lower, upper));
    }

    template<typename T>
    Box(Point<T, dim> const& _lower, Point<T, dim> const& _upper)
        : lower{_lower}
        , upper{_upper}
    {
        PHARE_DEBUG_DO(verify<dim>(lower, upper));
    }

    template<typename T2>
    NO_DISCARD bool operator==(Box<T2, dim> const& box) const
    {
        return box.lower == lower && box.upper == upper;
    }
    NO_DISCARD bool operator!=(Box const& other) const { return !(*this == other); }

    NO_DISCARD auto operator*(Box const& other) const
    {
        Box intersection{other};
        for (auto idim = 0u; idim < dim; ++idim)
        {
            intersection.lower[idim] = std::max(lower[idim], other.lower[idim]);
            intersection.upper[idim] = std::min(upper[idim], other.upper[idim]);
            if (intersection.lower[idim] > intersection.upper[idim])
                return std::optional<Box>{std::nullopt};
        }
        return std::optional<Box>{intersection};
    }


    NO_DISCARD auto unsafe_intersection(Box const& other) const _PHARE_ALL_FN_ // no optional on gpu
    {
        Box intersection{other};
        for (auto idim = 0u; idim < dim; ++idim)
        {
            intersection.lower[idim] = std::max(lower[idim], other.lower[idim]);
            intersection.upper[idim] = std::min(upper[idim], other.upper[idim]);
            assert(intersection.lower[idim] < intersection.upper[idim]);
        }
        return intersection;
    }



    template<typename Modifier>
    auto& grow(Modifier const& size) _PHARE_ALL_FN_
    {
        lower -= size;
        upper += size;
        return *this;
    }


    template<typename Size>
    auto& shrink(std::array<Size, dim> const& size) _PHARE_ALL_FN_
    {
        return shrink(Point_t{size});
    }

    template<typename Modifier>
    auto& shrink(Modifier const& size) _PHARE_ALL_FN_
    {
        assert(size >= 0);
        lower += size;
        upper -= size;
        return *this;
    }


    NO_DISCARD auto shape() const _PHARE_ALL_FN_ { return upper - lower + 1; }
    NO_DISCARD std::size_t size() const _PHARE_ALL_FN_ { return core::product(shape()); }


    NO_DISCARD auto begin() _PHARE_ALL_FN_ { return iterator{this, lower}; }

    //   // since the 1D scan of the multidimensional box is done assuming C ordering
    //   // the end (in the sense of container.end()) is one beyond last for the last
    //   // direction only, previous dimensions have not reached the end.
    NO_DISCARD auto begin() const _PHARE_ALL_FN_ { return iterator{this, lower}; }
    NO_DISCARD auto end() _PHARE_ALL_FN_
    {
        static_assert(dim <= 3 and dim > 0);
        // following could maybe be a one liner?
        if constexpr (dim == 1)
        {
            return iterator{this, {upper[0] + 1}};
        }
        else if constexpr (dim == 2)
        {
            return iterator{this, {upper[0] + 1, upper[1] + 1}};
        }
        else
        {
            return iterator{this, {upper[0] + 1, upper[1] + 1, upper[2] + 1}};
        }
    }

    NO_DISCARD auto end() const _PHARE_ALL_FN_
    {
        static_assert(dim <= 3 and dim > 0);
        if constexpr (dim == 1)
        {
            return iterator{this, {upper[0] + 1}};
        }
        else if constexpr (dim == 2)
        {
            return iterator{this, {upper[0] + 1, upper[1] + 1}};
        }
        else
        {
            return iterator{this, {upper[0] + 1, upper[1] + 1, upper[2] + 1}};
        }
    }


    NO_DISCARD constexpr static std::size_t nbrRemainBoxes()
    {
        if constexpr (dim == 1)
        {
            return 2;
        }
        else if constexpr (dim == 2)
        {
            return 4;
        }
        else
            return 6;
    }


    std::vector<Box> remove(Box const& that) const;
    std::vector<Box> remove_all(auto&&... them) const;
    Box merge(Box const& that) const;
    Box remove_edge(Box const& that) const;
};



template<typename Type, std::size_t dim>
class box_iterator
{
public:
    box_iterator(Box<Type, dim> const* box,
                 Point<Type, dim> index = Point<Type, dim>{}) _PHARE_ALL_FN_ : box_{box},
                                                                               index_{index}
    {
    }

public:
    NO_DISCARD auto& operator*() const _PHARE_ALL_FN_ { return index_; }
    // NO_DISCARD auto operator->() const _PHARE_ALL_FN_ { return &index_; }


    void increment(std::size_t idim) _PHARE_ALL_FN_
    {
        index_[idim]++;
        if (idim == 0)
            return;
        if (index_[idim] == box_->upper[idim] + 1)
        {
            increment(idim - 1);
            if (index_[idim - 1] <= box_->upper[idim - 1])
                index_[idim] = box_->lower[idim];
        }
    }

    auto& operator++() _PHARE_ALL_FN_
    {
        increment(dim - 1);
        return *this;
    }


    auto& operator+(std::uint32_t s) _PHARE_ALL_FN_
    {
        auto lo = box_->upper;
        for (std::uint16_t d = 0; d < dim; ++d)
            lo[d] += 1 - box_->lower[d];
        if constexpr (dim == 3)
        {
            auto const prod = lo[0] * lo[1];
            auto const div  = s / prod;
            s -= div * prod;
            index_[2] += div;
        }
        if constexpr (dim > 1)
        {
            auto const div = s / lo[0];
            s -= div * lo[0];
            index_[1] += div;
        }
        index_[0] += s;
        return *this;
    }

    bool operator!=(box_iterator const& other) const _PHARE_ALL_FN_
    {
        return box_ != other.box_ or index_ != other.index_;
    }


private:
    Box<Type, dim> const* box_;
    Point<Type, dim> index_;
};


template<typename T, std::size_t s>
Box(Point<T, s> lower, Point<T, s> upper) -> Box<T, s>;

template<typename... Boxes>
struct boxes_iterator
{
    auto constexpr static N = sizeof...(Boxes);

    boxes_iterator(Boxes const&... boxes)
        : boxes{std::make_tuple(boxes...)}
    {
    }

    struct iterator
    {
        using Tuple_t = std::tuple<typename Boxes::iterator...>;
        static_assert(N == std::tuple_size_v<Tuple_t>);

        iterator(std::tuple<typename Boxes::iterator...> iterators)
            : its{iterators}
        {
        }

        iterator& operator++()
        {
            for_N<N>([&](auto i) { ++std::get<i>(its); });
            return *this;
        }

        auto operator*()
        {
            return for_N<N>([&](auto i) { return *std::get<i>(its); });
        }

        auto operator!=(iterator const& that) const
        {
            return for_N_any<N>([&](auto i) { return std::get<i>(its) != std::get<i>(that.its); });
        }

        Tuple_t its;
    };



    auto begin()
    {
        return iterator{for_N<N>([&](auto i) { return std::get<i>(boxes).begin(); })};
    }
    auto end()
    {
        return iterator{for_N<N>([&](auto i) { return std::get<i>(boxes).end(); })};
    }

    std::tuple<Boxes...> boxes;
};

template<typename Particle>
NO_DISCARD auto _PHARE_ALL_FN_ isIn(Particle const& particle,
                                    Box<int, Particle::dimension> const& box)
    -> decltype(isIn(particle.iCell(), box), bool())
{
    return isIn(particle.iCell(), box);
}

template<template<typename, std::size_t> typename Point, typename Type, std::size_t SIZE>
NO_DISCARD bool isIn(Point<Type, SIZE> const& point, Box<Type, SIZE> const& box) _PHARE_ALL_FN_
{
    for (std::size_t iDim = 0; iDim < SIZE; ++iDim)
        if (point[iDim] < box.lower[iDim] || point[iDim] > box.upper[iDim])
            return false;
    return true;
}


/** this overload of isIn takes a Point and a Container of boxes
 * and returns true if the Point is at least in one of the boxes.
 * Returns occurs at the first box the point is in.
 */
template<typename Boxes>
NO_DISCARD bool _PHARE_ALL_FN_ isIn(auto const& point, Boxes const& boxes)
    requires(is_span_like_v<Boxes>)
{
    for (auto const& box : boxes)
        if (isIn(point, box))
            return true;
    return false;
}



template<typename Type, std::size_t dim, typename OType>
Box<Type, dim> grow(Box<Type, dim> const& box, std::array<Type, dim> const& by) _PHARE_ALL_FN_
{
    auto copy{box};
    copy.grow(by);
    return copy;
}


template<typename Type, std::size_t dim, typename OType>
Box<Type, dim> grow(Box<Type, dim> const& box, OType const& size) _PHARE_ALL_FN_
{
    auto copy{box};
    copy.grow(size);
    return copy;
}


template<typename Type, std::size_t dim, typename T2>
NO_DISCARD Box<Type, dim> shrink(Box<Type, dim> const& box, T2 const& size) _PHARE_ALL_FN_
{
    auto copy{box};
    copy.shrink(size);
    return copy;
}



template<typename Type, std::size_t dim, typename Shifter>
NO_DISCARD Box<Type, dim> shift(Box<Type, dim> const& box, Shifter const& offset)
{
    auto copy{box};
    copy.lower += offset;
    copy.upper += offset;
    return copy;
}



template<std::uint8_t idx, typename Type, std::size_t dim>
NO_DISCARD Box<Type, dim> shift_idx(Box<Type, dim> const& box, Type const& offset)
{
    auto copy{box};
    copy.lower[idx] += offset;
    copy.upper[idx] += offset;
    return copy;
}




template<typename Type, std::size_t dim>
auto& operator<<(std::ostream& os, Box<Type, dim> const& box)
{
    os << "Box<Type," << dim << "> : ( ";
    for (auto& c : box.lower)
        os << c << " ";
    os << ")-->( ";
    for (auto& c : box.upper)
        os << c << " ";
    os << ")";
    return os;
}



template<typename Type, std::size_t dim>
std::vector<Box<Type, dim>> Box<Type, dim>::remove(Box<Type, dim> const& to_remove) const
{
    using box_t = Box<Type, dim>;
    using _m    = std::unordered_map<std::uint16_t, std::uint32_t>;


    auto const box = *this;
    auto overlap   = box * to_remove;

    if (not overlap)
        return std::vector{*this};


    auto copy = [](auto cpy, auto const& replace) {
        for (auto const& [i, v] : replace)
            cpy[i] = v;
        return cpy;
    };

    auto intersection = *overlap;


    std::unordered_map<std::string, box_t> boxes;

    if (intersection.lower[0] > box.lower[0])
        boxes["left"] = Box(box.lower, copy(box.upper, _m{{0, intersection.lower[0] - 1}}));
    if (intersection.upper[0] < box.upper[0])
        boxes["right"] = box_t{copy(box.lower, _m{{0, intersection.upper[0] + 1}}), box.upper};

    [[maybe_unused]] Type minx = 0, maxx = 0;
    if constexpr (dim > 1)
    {
        minx = boxes.count("left") > 0 ? intersection.lower[0] : box.lower[0];
        maxx = boxes.count("right") > 0 ? intersection.upper[0] : box.upper[0];

        if (intersection.lower[1] > box.lower[1])
            boxes["down"] = box_t{copy(box.lower, _m{{0, minx}}),
                                  copy(box.upper, _m{{0, maxx}, {1, intersection.lower[1] - 1}})};

        if (intersection.upper[1] < box.upper[1])
            boxes["up"] = Box(copy(box.lower, _m{{0, minx}, {1, intersection.upper[1] + 1}}),
                              copy(box.upper, _m{{0, maxx}}));
    }

    if constexpr (dim > 2)
    {
        Type miny = boxes.count("down") > 0 ? intersection.lower[1] : box.lower[1];
        Type maxy = boxes.count("up") > 0 ? intersection.upper[1] : box.upper[1];

        if (intersection.lower[2] > box.lower[2])
            boxes["back"] = Box(copy(box.lower, _m{{0, minx}, {1, miny}}),
                                copy(intersection.lower - 1, _m{{0, maxx}, {1, maxy}}));
        if (intersection.upper[2] < box.upper[2])
            boxes["front"] = Box(copy(intersection.upper + 1, _m{{0, minx}, {1, miny}}),
                                 copy(box.upper, _m{{0, maxx}, {1, maxy}}));
    }

    std::vector<box_t> remaining;
    for (auto const& [key, val] : boxes)
        remaining.emplace_back(val);
    return remaining;
}


template<typename Type, std::size_t dim>
std::vector<Box<Type, dim>> Box<Type, dim>::remove_all(auto&&... them) const
{
    auto const add_to_vec = [](auto& vec, auto&& res) {
        for (auto const& r : res)
            vec.emplace_back(r);
    };

    auto const&& arr = std::array{them...};

    std::vector<Box<Type, dim>> remaining = this->remove(arr[0]);
    for (std::uint16_t a = 1; a < arr.size(); ++a)
    {
        auto const& to_rm = arr[a];
        std::vector<Box<Type, dim>> tmp;
        std::vector<std::size_t> remove;

        for (std::uint16_t r = 0; r < remaining.size(); ++r)
        {
            auto const& rem = remaining[r];
            if (rem * to_rm)
            {
                remove.emplace_back(r);
                add_to_vec(tmp, rem.remove(to_rm));
            }
        }

        for (auto it = remove.rbegin(); it != remove.rend(); ++it)
            remaining.erase(remaining.begin() + *it);

        add_to_vec(remaining, tmp);
    }

    return remaining;
}


template<typename Type, std::size_t dim>
Box<Type, dim> Box<Type, dim>::merge(Box<Type, dim> const& to_merge) const
{
    Box<Type, dim> merged{*this};
    for (std::size_t i = 0; i < dim; ++i)
    {
        merged.lower[i] = std::min(merged.lower[i], to_merge.lower[i]);
        merged.upper[i] = std::max(merged.upper[i], to_merge.upper[i]);
    }
    return merged;
}

template<typename Type, std::size_t dim>
Box<Type, dim> Box<Type, dim>::remove_edge(Box<Type, dim> const& to_remove) const
{
    auto const remaining = this->remove(to_remove);

    if (remaining.size() == 0)
        return *this;
    if (remaining.size() == 1)
        return remaining[0];

    Box<Type, dim> leftover{ConstArray<Type, dim>(std::numeric_limits<Type>::max()),
                            ConstArray<Type, dim>(std::numeric_limits<Type>::min())};

    for (auto const& rem : remaining)
        leftover = leftover.merge(rem);

    return leftover;
}

template<typename Boxes>
bool any_overlaps(Boxes const& boxes)
{
    for (std::size_t i = 0; i < boxes.size() - 1; ++i)
        for (std::size_t j = i + 1; j < boxes.size(); ++j)
            if (auto overlap = boxes[i] * boxes[j])
                return true;
    return false;
}

template<typename Boxes>
bool any_overlaps(Boxes const& boxes, typename Boxes::value_type const& box)
{
    for (std::size_t i = 0; i < boxes.size(); ++i)
        if (auto overlap = boxes[i] * box)
            return true;
    return false;
}


template<typename BoxHavers, typename Accessor>
bool any_overlaps_in(BoxHavers const& havers, Accessor&& fn)
{
    if (havers.size())
        for (std::size_t i = 0; i < havers.size() - 1; ++i)
            for (std::size_t j = i + 1; j < havers.size(); ++j)
                if (auto overlap = fn(havers[i]) * fn(havers[j]))
                    return true;
    return false;
}




template<typename Boxes, typename Box>
bool all_overlaps(Boxes const& boxes, Box const& box)
{
    std::uint32_t overlaps = 0;
    for (std::size_t i = 0; i < boxes.size(); ++i)
        if (auto overlap = boxes[i] * box)
            ++overlaps;
    return overlaps == boxes.size();
}

template<typename Boxes>
Boxes distinct_overlaps(Boxes& boxes)
{
    for (std::size_t i = 0; i < boxes.size() - 1; ++i)
    {
        auto const& a = boxes[i];

        for (std::size_t j = i + 1; j < boxes.size(); ++j)
        {
            auto const& b = boxes[j];

            if (auto overlap = a * b)
            {
                auto remaining0 = a.remove(b);
                auto remaining1 = b.remove(a);

                boxes.insert(boxes.end(), remaining0.begin(), remaining0.end());
                boxes.insert(boxes.end(), remaining1.begin(), remaining1.end());

                boxes[i] = *overlap;
                boxes.erase(boxes.begin() + j);

                return distinct_overlaps(boxes); // dragons
            }
        }
    }

    return boxes;
}

template<typename Boxes, typename Box>
auto distinct_overlaps(Boxes const& boxes, Box const& box)
{
    Boxes overlaps;

    for (std::size_t i = 0; i < boxes.size(); ++i)
        if (auto overlap = boxes[i] * box)
            overlaps.emplace_back(*overlap);

    return distinct_overlaps(overlaps);
}


template<template<typename, std::size_t> typename Point_t, typename Type, std::size_t SIZE>
auto box_from_zero_to_upper(Point_t<Type, SIZE> const& upper)
{
    return core::Box<std::uint32_t, SIZE>{
        core::Point<std::uint32_t, SIZE>{core::ConstArray<std::uint32_t, SIZE>()},
        core::Point<std::uint32_t, SIZE>{upper}};
}

template<template<typename, std::size_t> typename Point_t, typename Type, std::size_t SIZE>
auto box_from_zero_to_upper_minus_one(Point_t<Type, SIZE> const& upper) _PHARE_ALL_FN_
{
    return core::Box<std::uint32_t, SIZE>{
        core::Point<std::uint32_t, SIZE>{core::ConstArray<std::uint32_t, SIZE>()},
        core::Point<std::uint32_t, SIZE>{
            core::generate_from([](auto& u) { return u - 1; }, upper)}};
}

template<std::size_t SIZE>
NO_DISCARD auto as_unsigned(Box<int, SIZE> const& box) _PHARE_ALL_FN_
{
    return Box<std::uint32_t, SIZE>{box.lower.as_unsigned(), box.upper.as_unsigned()};
}

template<template<typename, std::size_t> typename Point_t, typename Type, std::size_t SIZE>
NO_DISCARD auto asBox(Point_t<Type, SIZE> const& point) _PHARE_ALL_FN_
{
    return Box<Type, SIZE>{point, point};
}



} // namespace PHARE::core



#endif
