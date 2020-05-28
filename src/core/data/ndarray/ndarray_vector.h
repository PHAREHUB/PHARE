#ifndef PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_H
#define PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_H

#include <stdexcept>
#include <array>
#include <cstdint>
#include <vector>
#include <tuple>
#include <numeric>


namespace PHARE
{
namespace core
{
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
                throw std::runtime_error(
                    "Error NdArrayVector cannot be assigned, incompatible sizes");
            }

            this->data_ = source.data_;
            return *this;
        }

        NdArrayVector& operator=(NdArrayVector&& source)
        {
            if (nCells_ != source.nCells_)
            {
                throw std::runtime_error(
                    "Error NdArrayVector cannot be assigned, incompatible sizes");
            }

            this->data_ = std::move(source.data_);
            return *this;
        }


        template<typename... Indexes>
        DataType const& operator()(Indexes... indexes) const
        {
            auto params = std::tuple<Indexes...>{indexes...};
            static_assert(sizeof...(Indexes) == dim);
            // static_assert((... && std::is_unsigned_v<decltype(indexes)>)); TODO : manage later if
            // this test should be included

            if constexpr (dim == 1)
            {
                auto i = std::get<0>(params);

                return this->data_[i];
            }

            if constexpr (dim == 2)
            {
                auto i = std::get<0>(params);
                auto j = std::get<1>(params);

                return this->data_[j + i * nCells_[1]];
            }

            if constexpr (dim == 3)
            {
                auto i = std::get<0>(params);
                auto j = std::get<1>(params);
                auto k = std::get<2>(params);

                return this->data_[k + j * nCells_[2] + i * nCells_[1] * nCells_[2]];
            }

            if constexpr (dim != 1 && dim != 2 && dim != 3)
            {
                throw std::runtime_error(
                    "Error NdArrayVector cannot be accessed, incompatible sizes");
            }
        }

        template<typename... Indexes>
        DataType& operator()(Indexes... indexes)
        {
            return const_cast<DataType&>(static_cast<NdArrayVector const&>(*this)(indexes...));
        }

        static const std::size_t dimension = dim;
        using type                         = DataType;

    private:
        std::array<std::uint32_t, dim> nCells_;
        std::vector<DataType> data_;
    };



    template<std::size_t dim>
    auto makeNdArray(std::array<std::uint32_t, dim> sizes)
    {
        if constexpr (dim == 1)
            return NdArrayVector<1>{sizes[0]};
        if constexpr (dim == 2)
            return NdArrayVector<2>{sizes[0], sizes[1]};
        if constexpr (dim == 3)
            return NdArrayVector<3>{sizes[0], sizes[1], sizes[2]};
    }


} // namespace core
} // namespace PHARE
#endif // PHARE_CORE_DATA_NDARRAY_NDARRAY_VECTOR_H
