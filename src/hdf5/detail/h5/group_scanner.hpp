#ifndef _PHARE_HDF5_DETAIL_H5_GROUP_SCANNER_HPP_
#define _PHARE_HDF5_DETAIL_H5_GROUP_SCANNER_HPP_

#include "core/def.hpp"
#include "core/def/phare_mpi.hpp"

#include "highfive/H5File.hpp"
#include "highfive/H5Easy.hpp"

#include "core/utilities/types.hpp"
#include "core/utilities/mpi_utils.hpp"
#include "core/utilities/meta/meta_utilities.hpp"

#include <unordered_set>

namespace PHARE::hdf5::h5
{

template<typename HighFiveFile_t>
struct GroupScanner
{
    GroupScanner(HighFiveFile_t const& h5_, std::vector<std::string> const& contains_ = {})
        : h5{h5_}
        , contains{contains_}
    {
    }


    auto& scan(std::string const& from = "/")
    {
        // prevent double / at start of full path
        scan(h5.file().getGroup(from), from == "/" ? "" : from);
        return groups;
    }
    void scan(HighFive::Group const& group, std::string const& path)
    {
        for (auto const& node : group.listObjectNames())
        {
            if (group.getObjectType(node) == HighFive::ObjectType::Group)
                scan(group.getGroup(node), path + "/" + node);
            else
            {
                auto fpath = path + "/" + node;
                if (contains.size() == 0)
                    groups.insert(fpath);
                else
                {
                    bool cont = false;
                    for (auto const& c : contains)
                        if (fpath.find(c) == std::string::npos)
                        {
                            cont = true;
                            break;
                        }
                    if (cont) // next node in listObjectNames
                        continue;
                    groups.insert(fpath);
                }
            }
        }
    }


    HighFiveFile_t const& h5;
    std::vector<std::string> contains;
    std::unordered_set<std::string> groups{};
};


} // namespace PHARE::hdf5::h5


#endif /* _PHARE_HDF5_DETAIL_H5_GROUP_SCANNER_HPP_ */
