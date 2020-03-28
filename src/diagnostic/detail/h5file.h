#ifndef H5FILE_H
#define H5FILE_H

#include "highfive/H5File.hpp"
#include "diagnostic/diagnostic_manager.h"

namespace PHARE::diagnostic::h5
{
using HiFile = HighFive::File;

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
};




} // namespace PHARE::diagnostic::h5


#endif // H5FILE_H
