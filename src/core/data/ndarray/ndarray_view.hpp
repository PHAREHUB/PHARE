#ifndef PHARE_CORE_DATA_NDARRAY_NDARRAY_VIEW_HPP
#define PHARE_CORE_DATA_NDARRAY_NDARRAY_VIEW_HPP

#include "core/def.hpp"
#include "core/utilities/types.hpp"
#include "core/data/ndarray/ndarray_base.hpp"
#include "core/data/ndarray/ndarray_mask.hpp"

#include <array>
#include <vector>
#include <cstdint>
#include <stdexcept>

namespace PHARE::core
{

template<std::size_t dim, typename DataType = double, bool c_ordering = true>
class NdArrayView
{
    using viewer = NdArrayViewer<dim, c_ordering>;

public:
    std::size_t static const dimension = dim;
    using type                         = DataType;
    using value_type                   = DataType;
    using pointer_type                 = DataType*;

    NdArrayView() = default;

    NdArrayView(pointer_type ptr, std::array<std::uint32_t, dim> const nCells)
        : ptr_{ptr}
        , size_{core::product(nCells)}
        , nCells_{nCells}
    {
    }

    NdArrayView(NdArrayView const&)            = default;
    NdArrayView& operator=(NdArrayView const&) = default;


    template<typename Index>
    NO_DISCARD inline auto& operator[](std::array<Index, dim> const& indexes)
    {
        return viewer::at(ptr_, nCells_, indexes);
    }
    template<typename Index>
    NO_DISCARD inline auto& operator[](std::array<Index, dim> const& indexes) const
    {
        return viewer::at(ptr_, nCells_, indexes);
    }

    template<typename Index>
    NO_DISCARD inline auto const& operator()(std::array<Index, dim> const& indexes) const
    {
        return viewer::at(ptr_, nCells_, indexes);
    }

    template<typename Index>
    NO_DISCARD inline auto& operator()(std::array<Index, dim> const& indexes)
    {
        return const_cast<DataType&>(static_cast<NdArrayView const&>(*this)(indexes));
    }

    NO_DISCARD inline auto const& operator()(auto const... indexes) const
    {
        return viewer::at(ptr_, nCells_, indexes...);
    }

    inline auto& operator()(auto const... indexes) { return viewer::at(ptr_, nCells_, indexes...); }


    NO_DISCARD auto& data() const { return ptr_; }
    NO_DISCARD auto& data() { return ptr_; }

    NO_DISCARD auto& size() const { return size_; }
    NO_DISCARD auto& shape() const { return nCells_; }

    NO_DISCARD auto begin() const { return ptr_; }
    NO_DISCARD auto begin() { return ptr_; }

    NO_DISCARD auto end() const { return ptr_ + size_; }
    NO_DISCARD auto end() { return ptr_ + size_; }

    void zero() { fill(0); }
    auto zeros() const
    {
        return sum_from(*this, [](auto const e) { return e == 0 ? 1 : 0; });
    }

    auto& fill(DataType const& v)
    {
        std::fill(begin(), end(), v);
        return *this;
    }

    bool isclose(NdArrayView const& that, double diff = 1e-13) const
    {
        if (this->size() != that.size())
            return false;
        for (std::size_t i = 0; i < this->size(); ++i)
            if (!float_equals(this->data()[i], that.data()[i], diff))
                return false;
        return true;
    }

    bool operator==(NdArrayView const& that) const
    {
        if (this->size() != that.size())
            return false;
        for (std::size_t i = 0; i < this->size(); ++i)
            if (this->data()[i] != that.data()[i])
                return false;
        return true;
    }

    template<typename View>
    void reset(View& view)
    {
        this->ptr_    = view.data();
        this->size_   = view.size();
        this->nCells_ = view.nCells_;
    }
    template<typename Vec>
    void reset(Vec& vec, std::array<std::uint32_t, dim> const& nCells)
    {
        this->ptr_    = vec.data();
        this->size_   = vec.size();
        this->nCells_ = nCells;
    }

    void fill_from(NdArrayView const& that)
    {
        if (for_N_any<dim>([&](auto i) { return shape()[i] != that.shape()[i]; }))
            throw std::runtime_error("ArrayView::fill_from: Incompatible input shape");
        std::copy(that.data(), that.data() + size(), data());
    }

    void setBuffer(pointer_type ptr) { ptr_ = ptr; }
    void setShape(std::array<std::uint32_t, dim> const nCells)
    {
        nCells_ = nCells;
        size_   = core::product(nCells);
    }

    auto& reshape(auto const& shape)
    {
        setShape(shape);
        return *this;
    }

    auto size_address() { return &size_; }

    NO_DISCARD auto operator[](NdArrayMask&& mask)
    {
        return MaskedView{*this, std::forward<NdArrayMask>(mask)};
    }
    NO_DISCARD auto operator[](NdArrayMask const& mask) { return MaskedView{*this, mask}; }

private:
    pointer_type ptr_                      = nullptr;
    std::size_t size_                      = 0;
    std::array<std::uint32_t, dim> nCells_ = {};
};


// Like NdArrayView but steps by sizeof(T) rather than sizeof(V) when indexing.
// V must be at offset 0 within T (i.e. V is the first base of T).
// This mirrors the ViewSpan::hax pattern used in GridTileSet.
template<std::size_t dim, typename V, typename T, bool c_ordering = true>
class NdArrayViewSpan
{
    using viewer = NdArrayViewer<dim, c_ordering>;

    template<typename V2, typename T2>
    static V2& hax(V2* data, std::uint32_t const i)
    {
        return *reinterpret_cast<V2*>(reinterpret_cast<T2*>(data) + i);
    }

public:
    std::size_t static const dimension = dim;
    using type                         = V;
    using value_type                   = V;

    NdArrayViewSpan() = default;

    NdArrayViewSpan(V* ptr, std::array<std::uint32_t, dim> const& nCells)
        : ptr_{ptr}
        , size_{core::product(nCells)}
        , nCells_{nCells}
    {
    }

    NdArrayViewSpan(NdArrayViewSpan const&)            = default;
    NdArrayViewSpan& operator=(NdArrayViewSpan const&) = default;

    template<typename Index>
    NO_DISCARD inline auto& operator()(std::array<Index, dim> const& indexes)
    {
        return hax<V, T>(ptr_, viewer::idx(nCells_, indexes));
    }
    template<typename Index>
    NO_DISCARD inline auto const& operator()(std::array<Index, dim> const& indexes) const
    {
        return hax<V const, T const>(ptr_, viewer::idx(nCells_, indexes));
    }
    NO_DISCARD inline auto& operator()(auto const... indexes)
    {
        return hax<V, T>(ptr_, viewer::idx(nCells_, indexes...));
    }
    NO_DISCARD inline auto const& operator()(auto const... indexes) const
    {
        return hax<V const, T const>(ptr_, viewer::idx(nCells_, indexes...));
    }

    NO_DISCARD auto data() { return ptr_; }
    NO_DISCARD auto data() const { return ptr_; }
    NO_DISCARD auto& size() const { return size_; }
    NO_DISCARD auto& shape() const { return nCells_; }

    // Iteration yields T& (full derived type) — correct step size, full CPU interface.
    NO_DISCARD auto begin() { return reinterpret_cast<T*>(ptr_); }
    NO_DISCARD auto begin() const { return reinterpret_cast<T const*>(ptr_); }
    NO_DISCARD auto end() { return reinterpret_cast<T*>(ptr_) + size_; }
    NO_DISCARD auto end() const { return reinterpret_cast<T const*>(ptr_) + size_; }

private:
    V* ptr_                                = nullptr;
    std::size_t size_                      = 0;
    std::array<std::uint32_t, dim> nCells_ = {};
};


template<bool c_ordering = true, typename DataType, std::size_t dim>
auto make_array_view(DataType* data, std::array<std::uint32_t, dim> const shape)
{
    return NdArrayView<dim, DataType, c_ordering>{data, shape};
}

template<bool c_ordering = true, typename DataType, std::size_t dim>
auto make_array_view(DataType const* const data, std::array<std::uint32_t, dim> const shape)
{
    return NdArrayView<dim, DataType const, c_ordering>{data, shape};
}

template<typename DataType, std::size_t dim>
auto make_array_view(std::vector<DataType>& vec, std::array<std::uint32_t, dim> const shape)
{
    return NdArrayView<dim, DataType>{vec.data(), shape};
}


template<typename DataType, std::size_t dim>
auto make_array_view(std::vector<DataType> const& vec, std::array<std::uint32_t, dim> const shape)
{
    return NdArrayView<dim, DataType const>{vec.data(), shape};
}


template<typename T>
concept is_ndarray_c
    = requires(T* p) { []<std::size_t dim, typename D>(NdArrayView<dim, D> const*) {}(p); };

template<typename T>
inline constexpr bool is_ndarray_v = is_ndarray_c<T>;


} // namespace PHARE::core

#endif // PHARE_CORE_DATA_NDARRAY_NDARRAY_VIEW_HPP
