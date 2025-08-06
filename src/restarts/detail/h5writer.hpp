#ifndef PHARE_DETAIL_RESTART_HIGHFIVE_HPP
#define PHARE_DETAIL_RESTART_HIGHFIVE_HPP


#include "core/utilities/mpi_utils.hpp"

#include "restarts/restarts_props.hpp"
#include "restarts/restarts_manager.hpp"

#include "initializer/data_provider.hpp"

#include "hdf5/detail/h5/h5_file.hpp"

namespace PHARE::restarts::h5
{
template<typename Hierarchy, typename ResourceManager_t>
class Writer
{
public:
    using This = Writer<Hierarchy, ResourceManager_t>;


    Writer(Hierarchy const& hier, ResourceManager_t const& resman, std::string const filePath)
        : hierarchy_{hier}
        , resman_{resman}
        , path_{filePath}
    {
    }

    ~Writer() {}



    static auto make_unique(Hierarchy const& hier, ResourceManager_t const& resman,
                            initializer::PHAREDict const& dict)
    {
        std::string filePath = dict["filePath"].template to<std::string>();
        return std::make_unique<This>(hier, resman, filePath);
    }


    void dump(RestartsProperties const& properties, double timestamp)
    {
        auto restart_file = writeRestartFile(restartFilePathForTime(path_, timestamp));

        // write model patch_data_ids to file with highfive
        // SAMRAI restart files are PER RANK
        PHARE::hdf5::h5::HighFiveFile h5File{restart_file, HighFive::File::ReadWrite,
                                             /*para=*/false};

        auto const& patch_ids = patch_data_ids();
        h5File.create_data_set<int>("/phare/patch/ids", patch_ids.size());
        h5File.write_data_set("/phare/patch/ids", patch_ids);

        h5File.write_attribute(
            "/phare", "serialized_simulation",
            properties.fileAttributes["serialized_simulation"].template to<std::string>());

        core::mpi::barrier();
    }



private:
    NO_DISCARD auto writeRestartFile(std::string const& path) const
    {
        return hierarchy_.writeRestartFile(path);
    }

    NO_DISCARD auto static restartFilePathForTime(std::string path, double timestamp)
    {
        return Hierarchy::restartFilePathForTime(path, timestamp);
    }


    NO_DISCARD auto patch_data_ids() const { return resman_.restart_patch_data_ids(); }



    Hierarchy const& hierarchy_;
    ResourceManager_t const& resman_;
    std::string const path_;
};



} // namespace PHARE::restarts::h5

#endif /* PHARE_DETAIL_RESTART_HIGHFIVE_H */
