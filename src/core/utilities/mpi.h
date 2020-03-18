
#ifndef PHARE_CORE_UTILITIES_MPI_H
#define PHARE_CORE_UTILITIES_MPI_H

namespace PHARE::mpi
{
template<typename Data>
std::vector<Data> collect(Data const& data, int mpi_size = 0);

size_t max(size_t local, int mpi_size = 0)
{
    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    auto perMPI = collect(local, mpi_size);
    return *std::max_element(std::begin(perMPI), std::end(perMPI));
}

template<typename Data>
void _collect(Data const* const data, std::vector<Data>& values, size_t send = 1,
              size_t receive = 1)
{
    if constexpr (std::is_same_v<double, Data>)
        MPI_Allgather(data, send, MPI_DOUBLE, values.data(), receive, MPI_DOUBLE, MPI_COMM_WORLD);
    else if constexpr (std::is_same_v<float, Data>)
        MPI_Allgather(data, send, MPI_FLOAT, values.data(), receive, MPI_FLOAT, MPI_COMM_WORLD);
    else if constexpr (std::is_same_v<size_t, Data>)
        MPI_Allgather(data, send, MPI_UINT64_T, values.data(), receive, MPI_UINT64_T,
                      MPI_COMM_WORLD);
    else if constexpr (std::is_same_v<char, Data>)
        MPI_Allgather(data, send, MPI_CHAR, values.data(), receive, MPI_CHAR, MPI_COMM_WORLD);
    else
        static_assert("Unhandled MPI data type collection");
}

std::vector<std::string> collectStrings(std::string str, int mpi_size = 0,
                                        std::string null_str = "null")
{
    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    str             = str.empty() ? null_str : str;
    auto maxMPISize = max(str.size(), mpi_size);
    auto perMPI     = collect(str.size(), mpi_size);
    std::vector<char> chars(maxMPISize * mpi_size);
    _collect(str.c_str(), chars, str.size(), maxMPISize);
    std::vector<std::string> values;
    for (int i = 0; i < mpi_size; i++)
    {
        std::string data{&chars[maxMPISize * i], perMPI[i]};
        data = data == null_str ? "" : data;
        values.emplace_back(data);
    }
    return values;
}


template<typename Data>
std::vector<std::vector<Data>> collectVector(std::vector<Data> const& data, int mpi_size)
{
    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    auto maxMPISize = max(data.size(), mpi_size);
    auto perMPI     = collect(data.size(), mpi_size);
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
    else
        _collect(&data, values);
    return values;
}


} // namespace PHARE::mpi


#endif /* PHARE_CORE_UTILITIES_MPI_H */
