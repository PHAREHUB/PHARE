#ifndef PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_H
#define PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_H

#include <stdexcept>
#include <array>
#include <cstdint>
#include <vector>
#include <tuple>
#include <numeric>


namespace PHARE::core
{
template<std::size_t dim, typename DataType = double>
struct NdArrayViewer
{
    template<typename NCells, typename... Indexes>
    static DataType const& at(DataType const* data, NCells const& nCells, Indexes const&... indexes)
    {
        auto params = std::forward_as_tuple(indexes...);
        static_assert(sizeof...(Indexes) == dim);
        // static_assert((... && std::is_unsigned_v<decltype(indexes)>)); TODO : manage later if
        // this test should be included

        if constexpr (dim == 1)
        {
            auto i = std::get<0>(params);

            return data[i];
        }

        if constexpr (dim == 2)
        {
            auto i = std::get<0>(params);
            auto j = std::get<1>(params);

            return data[j + i * nCells[1]];
        }

        if constexpr (dim == 3)
        {
            auto i = std::get<0>(params);
            auto j = std::get<1>(params);
            auto k = std::get<2>(params);

            return data[k + j * nCells[2] + i * nCells[1] * nCells[2]];
        }
    }

    template<typename NCells, typename Index>
    static DataType const& at(DataType const* data, NCells const& nCells,
                              std::array<Index, dim> const& indexes)

    {
        if constexpr (dim == 1)
            return data[indexes[0]];

        else if constexpr (dim == 2)
            return data[indexes[1] + indexes[0] * nCells[1]];

        else if constexpr (dim == 3)
            return data[indexes[2] + indexes[1] * nCells[2] + indexes[0] * nCells[1] * nCells[2]];
    }
};


template<std::size_t dim, typename DataType = double, typename Pointer = DataType const*>
class NdArrayView : NdArrayViewer<dim, DataType>
{
public:
    explicit NdArrayView(Pointer ptr, std::array<std::uint32_t, dim> const& nCells)
        : ptr_{ptr}
        , nCells_{nCells}
    {
    }

    explicit NdArrayView(std::vector<DataType> const& v,
                         std::array<std::uint32_t, dim> const& nbCell)
        : NdArrayView{v.data(), nbCell}
    {
    }

    template<typename... Indexes>
    DataType const& operator()(Indexes... indexes) const
    {
        return NdArrayViewer<dim, DataType>::at(ptr_, nCells_, indexes...);
    }

    template<typename... Indexes>
    DataType& operator()(Indexes... indexes)
    {
        return const_cast<DataType&>(static_cast<NdArrayView const&>(*this)(indexes...));
    }

    template<typename Index>
    DataType const& operator()(std::array<Index, dim> const& indexes) const
    {
        return NdArrayViewer<dim, DataType>::at(ptr_, nCells_, indexes);
    }

    template<typename Index>
    DataType& operator()(std::array<Index, dim> const& indexes)
    {
        return const_cast<DataType&>(static_cast<NdArrayView const&>(*this)(indexes));
    }

private:
    Pointer ptr_ = nullptr;
    std::array<std::uint32_t, dim> nCells_;
};

template<std::size_t dim, typename DataType = double>
class NdArrayVector
{
public:
    NdArrayVector() = delete;

    template<typename... Nodes>
    NdArrayVector(Nodes... nodes)
        : nCells_{nodes...}
        , data_((... * nodes))
    {
    }

    explicit NdArrayVector(std::array<std::uint32_t, dim> const& ncells)
        : nCells_{ncells}
        , data_(std::accumulate(ncells.begin(), ncells.end(), 1, std::multiplies<int>()))
    {
    }


    NdArrayVector(NdArrayVector const& source) = default;
    NdArrayVector(NdArrayVector&& source)      = default;

    static constexpr bool is_contiguous = 1;

    auto data() const { return data_.data(); }

    auto size() const { return data_.size(); }

    auto begin() const { return std::begin(data_); }
    auto begin() { return std::begin(data_); }

    auto end() const { return std::end(data_); }
    auto end() { return std::end(data_); }

    void zero() { data_ = std::vector<DataType>(data_.size(), {0}); }


    NdArrayVector& operator=(NdArrayVector const& source)
    {
        if (nCells_ != source.nCells_)
        {
            throw std::runtime_error("Error NdArrayVector cannot be assigned, incompatible sizes");
        }

        this->data_ = source.data_;
        return *this;
    }

    NdArrayVector& operator=(NdArrayVector&& source)
    {
        if (nCells_ != source.nCells_)
        {
            throw std::runtime_error("Error NdArrayVector cannot be assigned, incompatible sizes");
        }

        this->data_ = std::move(source.data_);
        return *this;
    }

    template<typename... Indexes>
    DataType const& operator()(Indexes... indexes) const
    {
        return NdArrayViewer<dim, DataType>::at(data_.data(), nCells_, indexes...);
    }

    template<typename... Indexes>
    DataType& operator()(Indexes... indexes)
    {
        return const_cast<DataType&>(static_cast<NdArrayVector const&>(*this)(indexes...));
    }

    template<typename Index>
    DataType const& operator()(std::array<Index, dim> const& indexes) const
    {
        return NdArrayViewer<dim, DataType>::at(data_.data(), nCells_, indexes);
    }

    template<typename Index>
    DataType& operator()(std::array<Index, dim> const& indexes)
    {
        return const_cast<DataType&>(static_cast<NdArrayVector const&>(*this)(indexes));
    }

    static const std::size_t dimension = dim;
    using type                         = DataType;

private:
    std::array<std::uint32_t, dim> nCells_;
    std::vector<DataType> data_;
};


} // namespace PHARE::core

#endif // PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_H
