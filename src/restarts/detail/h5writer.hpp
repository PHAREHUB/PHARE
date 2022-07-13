
#ifndef PHARE_DETAIL_RESTART_HIGHFIVE_HPP
#define PHARE_DETAIL_RESTART_HIGHFIVE_HPP

#include "core/logger.hpp"
#include "core/utilities/types.hpp"

#include "restarts/restarts_props.hpp"

#include "hdf5/detail/h5/h5_file.hpp"

namespace PHARE::restarts::h5
{
template<typename ModelView>
class Writer
{
public:
    using This = Writer<ModelView>;

    template<typename Hierarchy, typename Model>
    Writer(Hierarchy const& hier, Model const& model, std::string const filePath)
        : path_{filePath}
        , modelView_{hier, model}
    {
    }

    ~Writer() {}


    template<typename Hierarchy, typename Model>
    static auto make_unique(Hierarchy const& hier, Model const& model,
                            initializer::PHAREDict const& dict)
    {
        std::string filePath = dict["filePath"].template to<std::string>();
        return std::make_unique<This>(hier, model, filePath);
    }


    void dump(RestartsProperties const& properties, double timestamp)
    {
        auto restart_file
            = modelView_.writeRestartFile(ModelView::restartFilePathForTime(path_, timestamp));

        // write model patch_data_ids to file with highfive
        // SAMRAI restart files are PER RANK
        PHARE::hdf5::h5::HighFiveFile h5File{restart_file, HighFive::File::ReadWrite,
                                             /*para=*/false};

        auto patch_ids = modelView_.patch_data_ids();
        h5File.create_data_set<int>("/phare/patch/ids", patch_ids.size());
        h5File.write_data_set("/phare/patch/ids", patch_ids);

        h5File.write_attribute(
            "/phare", "serialized_simulation",
            properties.fileAttributes["serialized_simulation"].template to<std::string>());

        core::mpi::barrier();
    }

    auto& modelView() { return modelView_; }


private:
    std::string const path_;
    ModelView const modelView_;
};



} // namespace PHARE::restarts::h5

#endif /* PHARE_DETAIL_RESTART_HIGHFIVE_H */
