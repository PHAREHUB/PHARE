#ifndef PHARE_CORE_DATA_NDARRAY_NDARRAY_BASE_HPP
#define PHARE_CORE_DATA_NDARRAY_NDARRAY_BASE_HPP

#include "core/def.hpp"
#include "core/utilities/types.hpp"

#include <cstdint>

namespace PHARE::core
{

template<std::size_t dim, bool c_ordering = true>
struct NdArrayViewer
{
    template<typename NCells, template<typename, std::size_t> typename Indexes, typename Index>
    static inline std::size_t idx(NCells const& nCells, Indexes<Index, dim> const& indexes)
    {
        return for_N_fold<dim>([&](auto... d) { return idx(nCells, indexes[d]...); });
    }

    static inline std::size_t idx(auto const& /*nCells*/, std::size_t const i) { return i; }

    static inline std::size_t idx(auto const& nCells, std::size_t const i, std::size_t const j)
    {
        if constexpr (c_ordering)
            return j + i * nCells[1];
        else
            return i + j * nCells[0];
    }

    static inline std::size_t idx(auto const& nCells, std::size_t const i, std::size_t const j,
                                  std::size_t const k)
    {
        if constexpr (c_ordering)
            return k + j * nCells[2] + i * nCells[1] * nCells[2];
        else
            return i + j * nCells[0] + k * nCells[1] * nCells[0];
    }

    template<template<typename, std::size_t> typename Indexes, typename Index>
    NO_DISCARD static inline auto& at(auto* data, auto const& nCells,
                                      Indexes<Index, dim> const& indexes)
    {
        auto const& i = idx(nCells, indexes);
        assert(i < product(nCells, std::size_t{1}));
        return data[i];
    }

    static inline auto& at(auto* data, auto const& nCells, auto const... indexes)
    {
        auto const& i = idx(nCells, indexes...);
        assert(i < product(nCells, std::size_t{1}));
        return data[i];
    }
};

} // namespace PHARE::core

#endif // PHARE_CORE_DATA_NDARRAY_NDARRAY_BASE_HPP
