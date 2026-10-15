#ifndef PHARE_CORE_DATA_NDARRAY_NDARRAY_BASE_HPP
#define PHARE_CORE_DATA_NDARRAY_NDARRAY_BASE_HPP


#include "core/def.hpp"

#include <cstdint>


namespace PHARE::core
{
template<std::size_t dim, bool c_ordering = true>
struct NdArrayViewer
{
    template<typename NCells, template<typename, std::size_t> typename Indexes, typename Index>
    static inline std::uint32_t idx(NCells const& nCells,
                                    Indexes<Index, dim> const& indexes) _PHARE_ALL_FN_
    {
        if constexpr (dim == 1)
            return idx(nCells, indexes[0]);

        else if constexpr (dim == 2)
            return idx(nCells, indexes[0], indexes[1]);

        else if constexpr (dim == 3)
            return idx(nCells, indexes[0], indexes[1], indexes[2]);
    }

    static inline std::uint32_t idx(auto const /*nCells*/, std::uint32_t const i) _PHARE_ALL_FN_
    {
        return i;
    }


    static inline std::uint32_t idx(auto const nCells, std::uint32_t const i,
                                    std::uint32_t const j) _PHARE_ALL_FN_
    {
        if constexpr (c_ordering)
            return j + i * nCells[1];
        else
            return i + j * nCells[0];
    }
    static inline std::uint32_t idx(auto const nCells, std::uint32_t const i, std::uint32_t const j,
                                    std::uint32_t const k) _PHARE_ALL_FN_
    {
        if constexpr (c_ordering)
            return k + j * nCells[2] + i * nCells[1] * nCells[2];
        else
            return i + j * nCells[0] + k * nCells[1] * nCells[0];
    }



    template<template<typename, std::size_t> typename Indexes, typename Index>
    NO_DISCARD static inline auto& at(auto* data, auto const& nCells,
                                      Indexes<Index, dim> const& indexes) _PHARE_ALL_FN_

    {
        auto const& i = idx(nCells, indexes);
        assert(i < product(nCells, std::uint32_t{1}));
        return data[i];
    }

    static inline auto& at(auto* data, auto const nCells, auto const... indexes) _PHARE_ALL_FN_
    {
        auto const& i = idx(nCells, indexes...);
        assert(i < product(nCells, std::uint32_t{1}));
        return data[i];
    }
};




} // namespace PHARE::core

#endif // PHARE_CORE_DATA_NDARRAY_NDARRAY_BASE_HPP
