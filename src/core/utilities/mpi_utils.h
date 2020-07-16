
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


namespace PHARE::core::mpi
{
template<typename Data>
std::vector<Data> collect(Data const& data, int mpi_size = 0);

size_t max(size_t local, int mpi_size = 0);

int size();

template<typename Data, typename GatherFunc>
void _gather(GatherFunc&& gather)
{
    if constexpr (std::is_same_v<double, Data>)
        gather(MPI_DOUBLE);
    else if constexpr (std::is_same_v<float, Data>)
        gather(MPI_FLOAT);
    else if constexpr (std::is_same_v<int, Data>)
        gather(MPI_INT);
    else if constexpr (std::is_same_v<uint32_t, Data>)
        gather(MPI_UNSIGNED);
    else if constexpr (std::is_same_v<uint8_t, Data>)
        gather(MPI_UNSIGNED_SHORT);
    else if constexpr (std::is_same_v<size_t, Data>)
        gather(MPI_UINT64_T);
    else if constexpr (std::is_same_v<char, Data>)
        gather(MPI_CHAR);
    else
        throw std::runtime_error("Unhandled MPI data type collection");
}

template<typename Data>
void _collect(Data const* const data, std::vector<Data>& values, size_t send = 1,
              size_t receive = 1)
{
    _gather<Data>([&](auto mpi_type) {
        MPI_Allgather(     // MPI_Allgather
            data,          //   void         *sendbuf,
            send,          //   int          sendcount,
            mpi_type,      //   MPI_Datatype sendtype,
            values.data(), //   void         *recvbuf,
            receive,       //   int          recvcount,
            mpi_type,      //   MPI_Datatype recvtype,
            MPI_COMM_WORLD //   MPI_Comm     comm

        );
    });
}



template<typename Send, typename Data>
void _collect_vector(Send const& send, std::vector<Data>& receive, int mpi_size)
{
    _gather<Data>([&](auto mpi_type) {
        auto displs     = std::vector<int>(mpi_size);
        auto recvcounts = std::vector<int>(mpi_size);

        {
            int offset = 0;
            auto sizes = core::mpi::collect(send.size(), mpi_size);

            for (int i = 0; i < mpi_size; i++)
            {
                displs[i]     = offset;
                recvcounts[i] = sizes[i];
                offset += sizes[i];
            }
        }

        MPI_Allgatherv(        // MPI_Allgatherv
            send.data(),       //   void         *sendbuf,
            send.size(),       //   int          sendcount,
            mpi_type,          //   MPI_Datatype sendtype,
            receive.data(),    //   void         *recvbuf,
            recvcounts.data(), //   int          *recvcounts,
            displs.data(),     //   int          *displs,
            mpi_type,          //   MPI_Datatype recvtype,
            MPI_COMM_WORLD     //   MPI_Comm     comm
        );
    });
}


template<typename Vector>
std::vector<Vector> collectVector(Vector const& send, int mpi_size = 0)
{
    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    auto perMPI     = collect(send.size(), mpi_size);
    auto maxMPISize = *std::max_element(perMPI.begin(), perMPI.end());

    std::vector<typename Vector::value_type> receive(maxMPISize * mpi_size);
    _collect_vector(send, receive, mpi_size);

    std::vector<Vector> collected;
    for (int i = 0; i < mpi_size; i++)
    {
        auto* data = &receive[maxMPISize * i];
        collected.emplace_back(data, data + perMPI[i]);
    }
    return collected;
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

    auto perMPI = collect(data.size(), mpi_size);

    std::vector<typename Container::value_type> datas(maxMPISize * mpi_size);
    _collect(data.data(), datas, data.size(), maxMPISize);
    if constexpr (core::is_std_vector_v<Container>)
    {
        std::vector<Container> values;
        for (int i = 0; i < mpi_size; i++)
        {
            auto* array = &datas[maxMPISize * i];
            values.emplace_back(array, array + perMPI[i]);
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



template<typename Data>
std::vector<Data> collect(Data const& data, int mpi_size)
{
    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    std::vector<Data> values(mpi_size);

    if constexpr (std::is_same_v<std::string, Data> or core::is_std_vector_v<Data>)
        values = collectVector(data, mpi_size);
    else if constexpr (core::is_std_array_v<Data, 1>)
        values = collectArrays(data, mpi_size);
    else
        _collect(&data, values);
    return values;
}
} // namespace PHARE::core::mpi


#endif /* PHARE_CORE_UTILITIES_MPI_H */
