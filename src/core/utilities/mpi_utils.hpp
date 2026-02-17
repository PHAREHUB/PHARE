#ifndef PHARE_CORE_UTILITIES_MPI_HPP
#define PHARE_CORE_UTILITIES_MPI_HPP

#include "core/def.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/utilities/span.hpp"
#include "core/utilities/types.hpp"

#include <vector>
#include <string>
#include <cassert>
#include <cstring>

namespace PHARE::core::mpi
{
template<typename Data>
NO_DISCARD std::vector<Data> collect(Data const& data, int mpi_size = 0);

NO_DISCARD auto max(auto const local);
NO_DISCARD auto min(auto const local);

NO_DISCARD bool any(bool);

NO_DISCARD int size();

NO_DISCARD int rank();

void barrier();

NO_DISCARD std::string date_time(std::string format = "%Y-%m-%d-%H:%M:%S");

NO_DISCARD std::int64_t unix_timestamp_now();

auto inline sum_on_rank_0(std::size_t const s)
{
    std::size_t localsum[1]  = {s};
    std::size_t globalsum[1] = {0};
    MPI_Reduce(localsum, globalsum, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    return static_cast<std::size_t>(globalsum[0]);
}

inline bool is_init()
{
    int flag = 0;
    MPI_Initialized(&flag);
    return flag > 0;
}

template<typename Data>
NO_DISCARD auto mpi_type_for()
{
    if constexpr (std::is_same_v<double, Data>)
        return MPI_DOUBLE;
    else if constexpr (std::is_same_v<float, Data>)
        return MPI_FLOAT;
    else if constexpr (std::is_same_v<int, Data>)
        return MPI_INT;
    else if constexpr (std::is_same_v<std::uint32_t, Data>)
        return MPI_UNSIGNED;
    else if constexpr (std::is_same_v<std::int64_t, Data>)
        return MPI_INT64_T;
    else if constexpr (std::is_same_v<std::uint8_t, Data>)
        return MPI_UNSIGNED_SHORT;
    else if constexpr (std::is_same_v<std::size_t, Data>)
        return MPI_UINT64_T;
    else if constexpr (std::is_same_v<char, Data>)
        return MPI_CHAR;

    // don't return anything = compile failure if tried to use this function
}

auto mpi_type_for(auto const& val)
{
    return mpi_type_for<std::decay_t<decltype(val)>>();
}


template<typename Fn, typename... Args>
auto all_get_from(int const& rank_, Fn&& fn, Args&&... args)
{
    using Data = std::decay_t<std::invoke_result_t<Fn&, Args&...>>;

    Data var;
    auto local_rank = rank();
    if (local_rank == rank_)
        var = fn(args...);
    void* data = &var;

    int count = 1; // default
    MPI_Datatype sendtype;
    if constexpr (std::is_same_v<std::string, Data> or core::is_std_vector_v<Data>)
    {
        sendtype = mpi_type_for<typename Data::value_type>();
        count    = all_get_from(rank_, [&]() { return var.size(); });
        if (local_rank != rank_)
            var.reserve(count);
        data = var.data();
    }
    else
        sendtype = mpi_type_for<Data>();

    MPI_Bcast(         // MPI_Bcast
        data,          //   void *buffer
        count,         //   int count
        sendtype,      //   MPI_Datatype sendtype
        0,             //   int root == 0
        MPI_COMM_WORLD //   MPI_Comm comm
    );
    return var;
}


template<typename Fn, typename... Args>
auto all_get_from_rank_0(Fn&& fn, Args&&... args)
{
    return all_get_from(0, std::forward<Fn>(fn), std::forward<Args>(args)...);
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
NO_DISCARD std::vector<Vector> collectVector(Vector const& sendBuff, int mpi_size = 0)
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
    for (int i = 0; i < mpi_size; ++i)
    {
        if (perMPISize[i] == 0) // NO OUT OF BOUNDS ON LAST RANK
            collected.emplace_back();
        else
            collected.emplace_back(&rcvBuff[offset], &rcvBuff[offset] + perMPISize[i]);
        offset += perMPISize[i];
    }
    return collected;
}

template<typename T, typename Vector>
NO_DISCARD SpanSet<T, int> collectSpanSet(Vector const& sendBuff, int mpi_size = 0)
{
    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    SpanSet<T, int> rcvBuff{collect(static_cast<int>(sendBuff.size()), mpi_size)};
    _collect_vector<T>(sendBuff, rcvBuff, rcvBuff.sizes, rcvBuff.displs, mpi_size);

    return rcvBuff;
}



template<typename T, std::size_t size>
NO_DISCARD auto collectArrays(std::array<T, size> const& arr, int mpi_size)
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
NO_DISCARD SpanSet<typename Vector::value_type, int> collect_raw(Vector const& data, int mpi_size)
{
    if (mpi_size == 0)
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    return collectSpanSet<typename Vector::value_type>(data, mpi_size);
}


template<typename Data>
NO_DISCARD std::vector<Data> collect(Data const& data, int mpi_size)
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

auto max_on_rank0(auto const local)
{
    auto global = local;
    MPI_Reduce(&local, &global, 1, mpi_type_for(local), MPI_MAX, 0, MPI_COMM_WORLD);
    return global;
}


auto min_on_rank0(auto const local)
{
    auto global = local;
    MPI_Reduce(&local, &global, 1, mpi_type_for(local), MPI_MIN, 0, MPI_COMM_WORLD);
    return global;
}


auto max(auto const local)
{
    auto const mpi_size = size();
    auto const perMPI   = collect(local, mpi_size);
    return *std::max_element(std::begin(perMPI), std::end(perMPI));
}


auto min(auto const local)
{
    auto global = local;
    MPI_Allreduce(&local, &global, 1, mpi_type_for(local), MPI_MIN, MPI_COMM_WORLD);
    return global;
}


} // namespace PHARE::core::mpi


#endif /* PHARE_CORE_UTILITIES_MPI_H */
