
#ifndef PHARE_CORE_UTILITIES_MPI_H
#define PHARE_CORE_UTILITIES_MPI_H

#include <vector>
#include <string>
#include <cassert>
#include <cstring>

// clang-format off
#include "initializer/pragma_disable.h"
DISABLE_WARNING(cast-function-type, bad-function-cast, 42)
#include "mpi.h"
ENABLE_WARNING(cast-function-type, bad-function-cast, 42)
// clang-format on

#include "core/utilities/span.h"
#include "core/utilities/types.h"

namespace PHARE::core::mpi
{
template<typename Data>
std::vector<Data> collect(Data const& data, int mpi_size = 0);

std::size_t max(std::size_t const local, int mpi_size = 0);

int size();

int rank();

void barrier();

template<typename Data>
auto mpi_type_for()
{
    if constexpr (std::is_same_v<double, Data>)
        return MPI_DOUBLE;
    else if constexpr (std::is_same_v<float, Data>)
        return MPI_FLOAT;
    else if constexpr (std::is_same_v<int, Data>)
        return MPI_INT;
    else if constexpr (std::is_same_v<std::uint32_t, Data>)
        return MPI_UNSIGNED;
    else if constexpr (std::is_same_v<uint8_t, Data>)
        return MPI_UNSIGNED_SHORT;
    else if constexpr (std::is_same_v<std::size_t, Data>)
        return MPI_UINT64_T;
    else if constexpr (std::is_same_v<char, Data>)
        return MPI_CHAR;
    
    // don't return anything = compile failure if tried to use this function
}

template<typename Data>
void _collect(Data const* const sendbuf, std::vector<Data>& rcvBuff,
              std::size_t const sendcount = 1, std::size_t const recvcount = 1)
{
    auto mpi_type = mpi_type_for<Data>();

    MPI_Allgather(      // MPI_Allgather
        sendbuf,        //   void         *sendbuf,
        sendcount,      //   int          sendcount,
        mpi_type,       //   MPI_Datatype sendtype,
        rcvBuff.data(), //   void         *recvbuf,
        recvcount,      //   int          recvcount,
        mpi_type,       //   MPI_Datatype recvtype,
        MPI_COMM_WORLD  //   MPI_Comm     comm
    );
}



template<typename Data, typename SendBuff, typename RcvBuff>
void _collect_vector(SendBuff const& sendBuff, RcvBuff& rcvBuff, std::vector<int> const& recvcounts,
                     std::vector<int> const& displs, int const mpi_size)
{
    auto mpi_type = mpi_type_for<Data>();
    
    assert(recvcounts.size() == displs.size() and static_cast<int>(displs.size()) == mpi_size);
    
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
}

template<typename Vector>
std::vector<Vector> collectVector(Vector const& sendBuff, int mpi_size = 0)
{
    using Data = typename Vector::value_type;

    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    std::vector<int> const perMPISize = collect(static_cast<int>(sendBuff.size()), mpi_size);
    std::vector<int> const displs     = core::displacementFrom(perMPISize);
    std::vector<Data> rcvBuff(std::accumulate(perMPISize.begin(), perMPISize.end(), 0));
    _collect_vector<Data>(sendBuff, rcvBuff, perMPISize, displs, mpi_size);

    std::size_t offset = 0;
    std::vector<Vector> collected;
    for (int i = 0; i < mpi_size; i++)
    {
        collected.emplace_back(&rcvBuff[offset], &rcvBuff[offset] + perMPISize[i]);
        offset += perMPISize[i];
    }
    return collected;
}

template<typename T, typename Vector>
SpanSet<T, int> collectSpanSet(Vector const& sendBuff, int mpi_size = 0)
{
    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    SpanSet<T, int> rcvBuff{collect(static_cast<int>(sendBuff.size()), mpi_size)};
    _collect_vector<T>(sendBuff, rcvBuff, rcvBuff.sizes, rcvBuff.displs, mpi_size);

    return rcvBuff;
}



template<typename T, std::size_t size>
auto collectArrays(std::array<T, size> const& arr, int mpi_size)
{
    using Array = std::array<T, size>;
    using Data  = typename Array::value_type;

    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    std::size_t maxMPISize = arr.size();
    std::vector<Data> datas(maxMPISize * mpi_size);
    _collect(arr.data(), datas, arr.size(), maxMPISize);

    std::vector<Array> values(mpi_size);
    for (int i = 0; i < mpi_size; i++)
        std::memcpy(&values[i], &datas[maxMPISize * i], maxMPISize);

    return values;
}


template<typename Vector>
SpanSet<typename Vector::value_type, int> collect_raw(Vector const& data, int mpi_size)
{
    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    return collectSpanSet<typename Vector::value_type>(data, mpi_size);
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
