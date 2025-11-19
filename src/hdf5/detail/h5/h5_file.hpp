#ifndef PHARE_HDF5_H5FILE_HPP
#define PHARE_HDF5_H5FILE_HPP

#include "core/def.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/utilities/types.hpp"
#include "core/utilities/mpi_utils.hpp"
#include "core/utilities/meta/meta_utilities.hpp"


#include "highfive/H5File.hpp"
#include "highfive/H5Easy.hpp"



namespace PHARE::hdf5::h5::detail
{

// https://support.hdfgroup.org/documentation/hdf5/latest/hdf5_chunking.html
static inline auto const CHUNK_SIZE = core::get_env_as("PHARE_H5_CHUNK_SIZE", std::size_t{1024});

} // namespace PHARE::hdf5::h5::detail

namespace PHARE::hdf5::h5
{
using HiFile = HighFive::File;
using FileOp = HighFive::File::AccessMode;


template<typename Data, std::size_t dim>
NO_DISCARD auto vector_for_dim()
{
    if constexpr (dim == 1)
        return std::vector<Data>();
    if constexpr (dim == 2)
        return std::vector<std::vector<Data>>();
    if constexpr (dim == 3)
        return std::vector<std::vector<std::vector<Data>>>();
}

class HighFiveFile
{
public:
    template<typename FileAccessProps>
    static auto createHighFiveFile(std::string const path, FileOp flags, bool para,
                                   FileAccessProps& fapl)
    {
        if (para)
        {
#if defined(H5_HAVE_PARALLEL)
            fapl.add(HighFive::MPIOFileAccess{MPI_COMM_WORLD, MPI_INFO_NULL});
#else
            std::cout << "WARNING: PARALLEL HDF5 not available" << std::endl;
            if (core::mpi::size() > 1)
            {
                throw std::runtime_error("HDF5 NOT PARALLEL!");
            }
#endif
        }
        return HiFile{path, flags, fapl};
    }

    HighFiveFile(std::string const path, FileOp flags = HiFile::ReadWrite, bool para = true)
        : fapl_{}
        , h5file_{createHighFiveFile(path, flags, para, fapl_)}
    {
    }

    ~HighFiveFile() {}

    NO_DISCARD HiFile& file() { return h5file_; }


    template<typename T, std::size_t dim = 1>
    NO_DISCARD auto read_data_set_flat(std::string path) const
    {
        std::vector<T> data(H5Easy::getSize(h5file_, path));
        h5file_.getDataSet(path).read_raw(data.data());
        return data;
    }

    template<typename T, std::size_t dim = 1>
    NO_DISCARD auto read_data_set(std::string path) const
    {
        auto data = vector_for_dim<T, dim>();
        h5file_.getDataSet(path).read(data);
        return data;
    }


    template<std::size_t dim = 1, typename Data>
    auto& write_data_set(std::string path, Data const& data)
    {
        h5file_.getDataSet(path).write(data);
        return *this;
    }

    template<std::size_t dim = 1, typename Data>
    auto& write_data_set_flat(std::string path, Data const& data)
    {
        h5file_.getDataSet(path).write_raw(data);
        return *this;
    }


    template<typename Type, typename Size>
    auto create_data_set(std::string const& path, Size const& dataSetSize)
    {
        if (exist(path))
            return h5file_.getDataSet(path);
        createGroupsToDataSet(path);
        return h5file_.createDataSet<Type>(path, HighFive::DataSpace(dataSetSize));
    }


    template<typename Type>
    auto create_chunked_data_set(auto const& path, auto const chunk, auto const& dataspace)
    {
        if (exist(path))
            return h5file_.getDataSet(path);
        createGroupsToDataSet(path);
        HighFive::DataSetCreateProps props;
        props.add(HighFive::Chunking{chunk});
        return h5file_.createDataSet(path, dataspace, HighFive::create_datatype<Type>(), props);
    }

    template<typename Type, std::size_t cols>
    auto create_resizable_2d_data_set(auto const& path)
    {
        return create_chunked_data_set<Type>(
            path, std::vector<hsize_t>{detail::CHUNK_SIZE, cols},
            HighFive::DataSpace({0, cols}, {HighFive::DataSpace::UNLIMITED, cols}));
    }

    template<typename Type>
    auto create_resizable_1d_data_set(auto const& path)
    {
        return create_chunked_data_set<Type>(
            path, std::vector<hsize_t>{detail::CHUNK_SIZE},
            HighFive::DataSpace({0}, {HighFive::DataSpace::UNLIMITED}));
    }


    /*
     * Communicate all dataset paths and sizes to all MPI process to allow each to create all
     * datasets independently. This is a requirement of HDF5.
     * in the case of a disparate number of datasets per MPI process, dataSetSize may be 0
     * such that the current process creates datasets for all other processes with non-zero
     * sizes. Recommended to use similar sized paths, if possible.
     */
    template<typename Type, typename Size>
    void create_data_set_per_mpi(std::string const& path, Size const& dataSetSize)
    {
        auto const mpi_size = core::mpi::size();
        auto const sizes    = core::mpi::collect(dataSetSize, mpi_size);
        auto const paths    = core::mpi::collect(path, mpi_size);
        for (int i = 0; i < mpi_size; i++)
        { // in the case an mpi node lacks something to write
            if (is_zero(sizes[i]) || paths[i] == "")
                continue;

            create_data_set<Type>(paths[i], sizes[i]);
        }
    }


    /*
     * Write attribute on all mpi cores, considered global, always the same path/key and value
     */
    template<typename Data>
    void write_attribute(std::string const& keyPath, std::string const& key, Data const& data)
    {
        // assumes all keyPaths and values are identical, and no null patches
        // clang-format off
        PHARE_DEBUG_DO(
            auto const paths = core::mpi::collect(keyPath, core::mpi::size());
            if (!core::all(paths, [&](auto const& path) { return path == paths[0]; }))
                throw std::runtime_error("Function does not support different paths per mpi core");
        )
        // clang-format on

        constexpr bool data_is_vector = core::is_std_vector_v<Data>;

        auto doAttribute = [&](auto node, auto const& _key, auto const& value) {
            if constexpr (data_is_vector)
                node.template createAttribute<typename Data::value_type>(
                        _key, HighFive::DataSpace::From(value))
                    .write(value);
            else
                node.template createAttribute<Data>(_key, HighFive::DataSpace::From(value))
                    .write(value);
        };

        if (h5file_.exist(keyPath)
            && h5file_.getObjectType(keyPath) == HighFive::ObjectType::Dataset)
        {
            if (!h5file_.getDataSet(keyPath).hasAttribute(key))
                doAttribute(h5file_.getDataSet(keyPath), key, data);
        }
        else // group
        {
            createGroupsToDataSet(keyPath + "/dataset");
            if (!h5file_.getGroup(keyPath).hasAttribute(key))
                doAttribute(h5file_.getGroup(keyPath), key, data);
        }
    }



    /*
     * Communicate all attribute paths and values to all MPI process to allow each to create all
     * attributes independently. This is a requirement of HDF5.
     * in the case of a disparate number of attributes per MPI process, path may be an empty string
     * such that the current process creates attributes for all other processes with non-zero
     * sizes. Recommended to use similar sized paths, if possible. key is always assumed to the be
     * the same
     */
    template<typename Data>
    void write_attributes_per_mpi(std::string path, std::string key, Data const& data)
    {
        constexpr bool data_is_vector = core::is_std_vector_v<Data>;

        auto doAttribute = [&](auto node, auto const& _key, auto const& value) {
            if constexpr (data_is_vector)
            {
                if (value.size())
                    node.template createAttribute<typename Data::value_type>(
                            _key, HighFive::DataSpace(value.size()))
                        .write(value.data());
            }
            else
                node.template createAttribute<Data>(_key, HighFive::DataSpace::From(value))
                    .write(value);
        };

        int const mpi_size = core::mpi::size();
        auto const values  = [&]() {
            if constexpr (data_is_vector)
                return core::mpi::collect_raw(data, mpi_size);
            else
                return core::mpi::collect(data, mpi_size);
        }();
        auto const paths = core::mpi::collect(path, mpi_size);

        for (int i = 0; i < mpi_size; i++)
        {
            std::string const keyPath = paths[i] == "null" ? "" : paths[i];
            if (keyPath.empty())
                continue;

            if (h5file_.exist(keyPath)
                && h5file_.getObjectType(keyPath) == HighFive::ObjectType::Dataset)
            {
                if (!h5file_.getDataSet(keyPath).hasAttribute(key))
                    doAttribute(h5file_.getDataSet(keyPath), key, values[i]);
            }
            else // group
            {
                createGroupsToDataSet(keyPath + "/dataset");
                if (!h5file_.getGroup(keyPath).hasAttribute(key))
                    doAttribute(h5file_.getGroup(keyPath), key, values[i]);
            }
        }
    }

    template<typename Attr = std::string>
    auto read_attribute(std::string const& path, std::string const& key)
    {
        auto at = h5file_.getGroup(path).getAttribute(key);
        Attr attr;
        at.read(attr);
        return attr;
    }

    bool exist(std::string const& s) const { return h5file_.exist(s); }
    auto getDataSet(std::string const& s)
    {
        if (!exist(s))
            throw std::runtime_error("Dataset does not exist: " + s);
        return h5file_.getDataSet(s);
    }

    HighFiveFile(HighFiveFile const&)            = delete;
    HighFiveFile(HighFiveFile&&)                 = delete;
    HighFiveFile& operator=(HighFiveFile const&) = delete;
    HighFiveFile& operator=(HighFiveFile&&)      = delete;

private:
    HighFive::FileAccessProps fapl_;
    HiFile h5file_;


    // during attribute/dataset creation, we currently don't require the parents of the group to
    // exist, but create all that are missing - this used to exist in highfive
    inline std::string getParentName(std::string const& path)
    {
        std::size_t idx = path.find_last_of("/");
        if (idx == std::string::npos or idx == 0)
            return "/";
        return path.substr(0, idx);
    }

    inline void createGroupsToDataSet(std::string const& path)
    {
        std::string group_name = getParentName(path);
        if (!h5file_.exist(group_name))
            h5file_.createGroup(group_name);
    }



    template<typename Size>
    static bool is_zero(Size size)
    {
        if constexpr (core::is_iterable_v<Size>)
            return std::all_of(size.begin(), size.end(), [](auto const& val) { return val == 0; });

        else
            return size == 0;
    }
};




} // namespace PHARE::hdf5::h5


#endif /* PHARE_HDF5_H5FILE_HPP */
