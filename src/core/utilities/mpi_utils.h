
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


template<typename Data>
void _collect(Data const* const data, std::vector<Data>& values, size_t send = 1,
              size_t receive = 1)
{
    auto gather = [&](auto mpi_type) {
        MPI_Allgather(data, send, mpi_type, values.data(), receive, mpi_type, MPI_COMM_WORLD);
    };

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

std::vector<std::string> collectStrings(std::string str, int mpi_size = 0,
                                        std::string null_str = "null");


template<typename Data>
std::vector<std::vector<Data>> collectVector(std::vector<Data> const& data, int mpi_size,
                                             bool empty = false)
{ // can't send empty vectors
    if (data.size() == 0)
    {
        std::vector<Data> data_(1);
        return collectVector(data_, mpi_size, true);
    }

    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    auto maxMPISize = max(data.size(), mpi_size);
    auto perMPI     = collect(empty ? 0 : data.size(), mpi_size);

    std::vector<Data> datas(maxMPISize * mpi_size);
    _collect(data.data(), datas, data.size(), maxMPISize);

    std::vector<std::vector<Data>> values;
    for (int i = 0; i < mpi_size; i++)
    {
        auto* array = &datas[maxMPISize * i];
        values.emplace_back(array, array + perMPI[i]);
    }
    return values;
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

    if constexpr (std::is_same_v<std::string, Data>)
        values = collectStrings(data, mpi_size);
    else if constexpr (core::is_std_vector_v<Data>)
        values = collectVector(data, mpi_size);
    else if constexpr (core::is_std_array_v<Data, 1>)
        values = collectArrays(data, mpi_size);
    else
        _collect(&data, values);
    return values;
}
} // namespace PHARE::core::mpi


#endif /* PHARE_CORE_UTILITIES_MPI_H */
