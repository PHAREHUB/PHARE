#ifndef H5FILE_H
#define H5FILE_H

#include "highfive/H5File.hpp"
#include "highfive/H5Easy.hpp"
#include "diagnostic/diagnostic_manager.hpp"

namespace PHARE::diagnostic::h5
{
using HiFile = HighFive::File;


/*
  Highfive cannot accept a single flat array into >= 2d shaped datasets
*/
template<std::size_t dim, typename Data>
auto pointer_dim_caster(Data* data)
{
    if constexpr (dim == 1)
        return data;
    if constexpr (dim == 2)
        return reinterpret_cast<Data const** const>(data);
    if constexpr (dim == 3)
        return reinterpret_cast<Data const*** const>(data);
}


template<std::size_t dim, typename Data>
auto decay_to_pointer(Data& data)
{
    if constexpr (dim == 1)
        return data.data();
    if constexpr (dim == 2)
        return data[0].data();
    if constexpr (dim == 3)
        return data[0][0].data();
}

template<typename Data, std::size_t dim>
auto vector_for_dim()
{
    if constexpr (dim == 1)
        return std::vector<Data>();
    if constexpr (dim == 2)
        return std::vector<std::vector<Data>>();
    if constexpr (dim == 3)
        return std::vector<std::vector<std::vector<Data>>>();
}

struct HighFiveFile
{
    HiFile file_;
    PHARE::diagnostic::Mode mode_ = PHARE::diagnostic::Mode::LIGHT;

    static auto createHighFiveFile(std::string const path, unsigned flags)
    {
        return HiFile
        {
            path, flags
#if defined(H5_HAVE_PARALLEL)
                ,
                HighFive::MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL)
#endif
        };
    }

    HighFiveFile(std::string const path, unsigned flags,
                 PHARE::diagnostic::Mode mode = PHARE::diagnostic::Mode::LIGHT)
        : file_{createHighFiveFile(path, flags)}
        , mode_{mode}
    {
    }
    ~HighFiveFile() {}

    HiFile& file() { return file_; }

    HighFiveFile(const HighFiveFile&)  = delete;
    HighFiveFile(const HighFiveFile&&) = delete;
    HighFiveFile& operator=(const HighFiveFile&) = delete;
    HighFiveFile& operator=(const HighFiveFile&&) = delete;


    template<typename T, std::size_t dim = 1>
    auto read_data_set_flat(std::string path) const
    {
        std::vector<T> data(H5Easy::getSize(file_, path));
        file_.getDataSet(path).read(pointer_dim_caster<dim>(data.data()));
        return data;
    }

    template<typename T, std::size_t dim = 1>
    auto read_data_set(std::string path) const
    {
        auto data = vector_for_dim<T, dim>();
        file_.getDataSet(path).read(data);
        return data;
    }


    template<std::size_t dim = 1, typename Data>
    auto& write_data_set(std::string path, Data const& data)
    {
        file_.getDataSet(path).write(data);
        return *this;
    }

    template<std::size_t dim = 1, typename Data>
    auto& write_data_set_flat(std::string path, Data const& data)
    {
        file_.getDataSet(path).write(pointer_dim_caster<dim>(data));
        return *this;
    }
};




} // namespace PHARE::diagnostic::h5


#endif // H5FILE_H
