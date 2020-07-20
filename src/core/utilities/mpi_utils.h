
#ifndef PHARE_CORE_UTILITIES_MPI_H
#define PHARE_CORE_UTILITIES_MPI_H

#include <vector>
#include <string>
#include <cstring>

// clang-format off
#include "initializer/pragma_disable.h"
DISABLE_WARNING(cast-function-type, bad-function-cast, 42)
#include "mpi.h"
ENABLE_WARNING(cast-function-type, bad-function-cast, 42)
// clang-format on

#include "core/utilities/types.h"
#include "kul/tuple.hpp"

namespace PHARE::core::mpi
{
template<typename Data>
std::vector<Data> collect(Data const& data, int mpi_size = 0);

size_t max(size_t const local, int mpi_size = 0);

int size();

template<typename Data, typename GatherFunc>
void _gather(GatherFunc const&& gather)
{
    if constexpr (std::is_same_v<double, Data>)
        gather(MPI_DOUBLE);
    else if constexpr (std::is_same_v<float, Data>)
        gather(MPI_FLOAT);
    else if constexpr (std::is_same_v<int, Data>)
        gather(MPI_INT);
    else if constexpr (std::is_same_v<std::uint32_t, Data>)
        gather(MPI_UNSIGNED);
    else if constexpr (std::is_same_v<uint8_t, Data>)
        gather(MPI_UNSIGNED_SHORT);
    else if constexpr (std::is_same_v<std::size_t, Data>)
        gather(MPI_UINT64_T);
    else if constexpr (std::is_same_v<char, Data>)
        gather(MPI_CHAR);
    else
        throw std::runtime_error(std::string("Unhandled MPI data type collection: ")
                                 + typeid(Data).name());
}

template<typename Data>
void _collect(Data const* const sendbuf, std::vector<Data>& rcvBuff,
              std::size_t const sendcount = 1, std::size_t const recvcount = 1)
{
    _gather<Data>([&](auto mpi_type) {
        MPI_Allgather(      // MPI_Allgather
            sendbuf,        //   void         *sendbuf,
            sendcount,      //   int          sendcount,
            mpi_type,       //   MPI_Datatype sendtype,
            rcvBuff.data(), //   void         *recvbuf,
            recvcount,      //   int          recvcount,
            mpi_type,       //   MPI_Datatype recvtype,
            MPI_COMM_WORLD  //   MPI_Comm     comm

        );
    });
}



template<typename Data, typename SendBuff, typename RcvBuff>
void _collect_vector(SendBuff const& sendBuff, RcvBuff& rcvBuff, std::vector<int> const& recvcounts,
                     int const mpi_size)
{
    _gather<Data>([&](auto const mpi_type) {
        std::vector<int> displs(mpi_size);

        {
            int offset = 0;

            for (int i = 0; i < mpi_size; i++)
            {
                displs[i] = offset;
                offset += recvcounts[i];
            }
        }

        MPI_Allgatherv(        // MPI_Allgatherv
            sendBuff.data(),   //   void         *sendbuf,
            sendBuff.size(),   //   int          sendcount,
            mpi_type,          //   MPI_Datatype sendtype,
            rcvBuff.data(),    //   void         *recvbuf,
            recvcounts.data(), //   int          *recvcounts,
            displs.data(),     //   int          *displs,
            mpi_type,          //   MPI_Datatype recvtype,
            MPI_COMM_WORLD     //   MPI_Comm     comm
        );
    });
}

template<typename Vector>
std::vector<Vector> collectVector(Vector const& sendBuff, int mpi_size = 0)
{
    using Data = typename Vector::value_type;

    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    auto const perMPISize = collect(static_cast<int>(sendBuff.size()), mpi_size);
    std::vector<Data> rcvBuff(std::accumulate(perMPISize.begin(), perMPISize.end(), 0));
    _collect_vector<Data>(sendBuff, rcvBuff, perMPISize, mpi_size);

    size_t offset = 0;
    std::vector<Vector> collected;
    for (int i = 0; i < mpi_size; i++)
    {
        collected.emplace_back(&rcvBuff[offset], &rcvBuff[offset] + perMPISize[i]);
        offset += perMPISize[i];
    }
    return collected;
}

template<typename T, typename Vector>
decltype(auto) collectSplitVector(Vector const& sendBuff, int mpi_size = 0)
{
    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    kul::SplitVector<T, int> rcvBuff{collect(static_cast<int>(sendBuff.size()), mpi_size)};
    _collect_vector<T>(sendBuff, rcvBuff, rcvBuff.sizes, mpi_size);

    return rcvBuff;
}



template<typename Container>
auto collectArrays(Container const& data, int mpi_size)
{
    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    std::size_t maxMPISize = 0;
    if constexpr (core::is_std_vector_v<Container>)
        maxMPISize = max(data.size(), mpi_size);
    else if constexpr (core::is_std_array_v<typename Container::value_type, 1>)
        maxMPISize = data.size();

    auto perMPISize = collect(data.size(), mpi_size);

    std::vector<typename Container::value_type> datas(maxMPISize * mpi_size);
    _collect(data.data(), datas, data.size(), maxMPISize);
    if constexpr (core::is_std_vector_v<Container>)
    {
        std::vector<Container> values;
        for (int i = 0; i < mpi_size; i++)
        {
            auto const* const array = &datas[maxMPISize * i];
            values.emplace_back(array, array + perMPISize[i]);
        }
        return values;
    }
    else
    {
        std::vector<Container> values(mpi_size);
        for (int i = 0; i < mpi_size; i++)
        {
            auto* array = &datas[maxMPISize * i];
            auto* valp  = &values[i];
            std::memcpy(valp, array, maxMPISize);
        }
        return values;
    }
}


template<typename T>
decltype(auto) collect_raw(std::vector<T> const& data, int mpi_size)
{
    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    return collectSplitVector<T>(data, mpi_size);
}


template<typename Data>
std::vector<Data> collect(Data const& data, int mpi_size)
{
    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    if constexpr (std::is_same_v<std::string, Data> or core::is_std_vector_v<Data>)
        return collectVector(data, mpi_size);
    else if constexpr (core::is_std_array_v<Data, 1>)
        return collectArrays(data, mpi_size);
    else
    {
        std::vector<Data> values(mpi_size);
        _collect(&data, values);
        return values;
    }
}
} // namespace PHARE::core::mpi


#endif /* PHARE_CORE_UTILITIES_MPI_H */
