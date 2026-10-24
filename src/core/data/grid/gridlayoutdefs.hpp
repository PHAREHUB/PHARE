#ifndef PHARE_CORE_GRID_GRIDLAYOUTDEFS_HPP
#define PHARE_CORE_GRID_GRIDLAYOUTDEFS_HPP

#include "core/def.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <type_traits>
#include <utility>

namespace PHARE
{
namespace core
{
    enum class Direction { X, Y, Z };

    enum class Side { Lower = -1, Upper = 1 };

    enum class QtyCentering : std::uint16_t { primal = 0, dual = 1 };


    template<std::size_t dim>
    struct WeightPoint
    {
        Point<int, dim> indexes{};
        double coef{0.0};
    };

    // using LinearCombination = std::vector<WeightPoint>;

    enum class Layout { Yee };


    template<typename T, std::size_t dim>
    NO_DISCARD Box<T, dim> toCenteredBox(Box<T, dim> box,
                                         std::array<QtyCentering, dim> const& centerings)
    {
        for (std::size_t i = 0; i < dim; ++i)
            if (centerings[i] == QtyCentering::primal)
                ++box.upper[i];

        return box;
    }

    enum class InterpDir { DualToPrimal = 0, PrimalToDual = 1 };


    NO_DISCARD constexpr std::array<int, 2> linearRecenteringOffsets(QtyCentering source)
    {
        if (source == QtyCentering::primal)
            return {0, 1};
        return {-1, 0};
    }


    template<std::size_t dim, auto dir, InterpDir interpDir, std::size_t order = 2>
    NO_DISCARD consteval auto directionalInterp()
    {
        if constexpr (dir >= dim)
        {
            return std::array{WeightPoint<dim>{Point<int, dim>{}, 1.0}};
        }
        else
        {
            constexpr int baseidx = (interpDir == InterpDir::PrimalToDual);

            auto make_p = [](auto offset) {
                Point<int, dim> p{};
                p[dir] = offset;
                return p;
            };

            if constexpr (order == 2)
            {
                constexpr auto offsets = linearRecenteringOffsets(
                    interpDir == InterpDir::PrimalToDual ? QtyCentering::primal
                                                         : QtyCentering::dual);
                auto w = 0.5;
                return std::array{WeightPoint<dim>{make_p(offsets[0]), w},
                                  WeightPoint<dim>{make_p(offsets[1]), w}};
            }
            else if constexpr (order == 4)
            {
                auto w1 = -1.0 / 16.0;
                auto w2 = 9.0 / 16.0;

                return std::array{WeightPoint<dim>{make_p(baseidx - 2), w1},
                                  WeightPoint<dim>{make_p(baseidx - 1), w2},
                                  WeightPoint<dim>{make_p(baseidx), w2},
                                  WeightPoint<dim>{make_p(baseidx + 1), w1}};
            }
            else if constexpr (order == 6)
            {
                auto w1 = 3.0 / 256.0;
                auto w2 = -25.0 / 256.0;
                auto w3 = 150.0 / 256.0;

                return std::array{WeightPoint<dim>{make_p(baseidx - 3), w1},
                                  WeightPoint<dim>{make_p(baseidx - 2), w2},
                                  WeightPoint<dim>{make_p(baseidx - 1), w3},
                                  WeightPoint<dim>{make_p(baseidx), w3},
                                  WeightPoint<dim>{make_p(baseidx + 1), w2},
                                  WeightPoint<dim>{make_p(baseidx + 2), w1}};
            }
            else
                static_assert(order == 2 || order == 4 || order == 6, "unsupported order");
        }
    }


    namespace detail
    {
        template<std::size_t dim, std::size_t n>
        constexpr std::size_t weightPointDim(std::array<WeightPoint<dim>, n> const&)
        {
            return dim;
        }
    } // namespace detail


    template<auto... dirs>
    NO_DISCARD consteval auto tensorProduct(auto const&... rows)
    {
        static_assert(sizeof...(dirs) == sizeof...(rows) && sizeof...(rows) > 0);
        static_assert(
            []() {
                constexpr std::array<std::size_t, sizeof...(dirs)> ds{
                    static_cast<std::size_t>(dirs)...};
                for (std::size_t i = 0; i < ds.size(); ++i)
                    for (std::size_t j = i + 1; j < ds.size(); ++j)
                        if (ds[i] == ds[j])
                            return false;
                return true;
            }(),
            "tensorProduct directions must be pairwise distinct");

        constexpr std::size_t dim
            = detail::weightPointDim(std::get<0>(std::tuple<std::decay_t<decltype(rows)>...>{}));
        constexpr std::size_t total = (std::tuple_size_v<std::decay_t<decltype(rows)>> * ...);

        std::array<WeightPoint<dim>, total> result{};
        for (std::size_t k = 0; k < total; ++k)
        {
            Point<int, dim> pt{};
            double coef        = 1.0;
            std::size_t stride = total;

            auto accumulate = [&](auto const& row, auto dirConstant) {
                stride /= row.size();
                auto const& wp = row[(k / stride) % row.size()];
                coef *= wp.coef;
                if constexpr (decltype(dirConstant)::value < dim)
                    pt[decltype(dirConstant)::value] = wp.indexes[decltype(dirConstant)::value];
            };

            (accumulate(rows, std::integral_constant<std::decay_t<decltype(dirs)>, dirs>{}), ...);
            result[k] = {pt, coef};
        }
        return result;
    }


    namespace detail
    {
        template<auto src, auto dst, std::size_t order, auto diff, std::size_t... k>
        NO_DISCARD consteval auto recenteringStencilImpl(std::index_sequence<k...>)
        {
            constexpr std::size_t dim = src.size();
            return tensorProduct<diff[k]...>(
                directionalInterp<dim, diff[k],
                                  (src[diff[k]] == QtyCentering::primal ? InterpDir::PrimalToDual
                                                                        : InterpDir::DualToPrimal),
                                  order>()...);
        }

        template<auto src, auto dst>
        NO_DISCARD consteval std::size_t differingDirections()
        {
            std::size_t count = 0;
            for (std::size_t d = 0; d < src.size(); ++d)
                count += (src[d] != dst[d]);
            return count;
        }

        template<auto src, auto dst, std::size_t p>
        NO_DISCARD consteval auto differingDirectionsArray()
        {
            std::array<std::size_t, p> diff{};
            std::size_t k = 0;
            for (std::size_t d = 0; d < src.size(); ++d)
                if (src[d] != dst[d])
                    diff[k++] = d;
            return diff;
        }
    } // namespace detail


    template<auto src, auto dst, std::size_t order = 2>
    NO_DISCARD consteval auto recenteringStencil()
    {
        static_assert(src.size() == dst.size());
        constexpr std::size_t dim = src.size();
        constexpr std::size_t p   = detail::differingDirections<src, dst>();

        if constexpr (p == 0)
        {
            return std::array{WeightPoint<dim>{Point<int, dim>{}, 1.0}};
        }
        else
        {
            constexpr auto diff = detail::differingDirectionsArray<src, dst, p>();
            return detail::recenteringStencilImpl<src, dst, order, diff>(
                std::make_index_sequence<p>{});
        }
    }


    template<auto src, auto dst, std::size_t order = 2, typename Field>
    NO_DISCARD auto project(Field const& field, Point<std::size_t, Field::dimension> const& index)
    {
        constexpr auto stencil = recenteringStencil<src, dst, order>();

        typename Field::type result = 0.;
        for (auto const& wp : stencil)
            result += wp.coef * field(index + wp.indexes);
        return result;
    }


} // namespace core
} // namespace PHARE

#endif // GRIDLAYOUTDEFS_HPP
