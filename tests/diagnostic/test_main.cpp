
#include <filesystem>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Easy.hpp>

#include "diagnostic/include.hpp"
#include "diagnostic/samrai_lifecycle.hpp"
#include "diagnostic/defs.hpp"
using namespace PHARE_test::_1d;
#include "diagnostic/integrator.hpp"
#include "diagnostic/tag_strat.hpp"
#include "diagnostic/hierarchy.hpp"

struct Hi5Diagnostic : public AfullHybridBasicHierarchy
{
    HighFive::File file{"/tmp/new_diagnostic.file.h5", HighFive::File::ReadWrite
                                                           | HighFive::File::Create
                                                           | HighFive::File::Truncate};
};

TEST_F(Hi5Diagnostic, initializesFieldsOnRefinedLevels)
{
    using FIELD_ref = decltype(hybridModel->state.electromag.E.getComponent(Component::X));
    using FIELD     = std::remove_reference_t<FIELD_ref>;
    auto& hifile    = this->file;
    auto& hierarchy = basicHierarchy->getHierarchy();
    auto& rm        = hybridModel->resourcesManager;

    hifile.createGroup("/timestep");
    size_t ts_idx = 0;
    [](HighFive::File& file, size_t ts_idx) {
        const auto ts_idx_str = std::to_string(ts_idx);
        std::stringstream ss;
        ss << "/timestep/" << ts_idx_str;
        auto ts_group = file.createGroup(ss.str());
        ts_group.createAttribute<std::string>("index", HighFive::DataSpace::From(ts_idx_str))
            .write(ts_idx_str);
        ss << "/patchlevel/";
        file.createGroup(ss.str());
    }(hifile, ts_idx);

    std::string _0_path("/timestep/" + std::to_string(ts_idx) + "/patchlevel/");
    H5Easy::detail::createGroupsToDataSet(hifile, _0_path);
    for (int iLevel = 0; iLevel < hierarchy.getNumberOfLevels(); ++iLevel)
    {
        std::stringstream levelPath;
        levelPath << _0_path << iLevel;
        file.createGroup(levelPath.str());

        auto const& level = hierarchy.getPatchLevel(iLevel);

        size_t patch_idx = 0;
        for (auto& patch : *level)
        {
            std::stringstream patchPath;
            patchPath << levelPath.str() << "/patch/" << patch_idx;
            file.createGroup(patchPath.str());

            auto onPatch
                = rm->setOnPatch(*patch, hybridModel->state.electromag, hybridModel->state.ions);

            auto layout = layoutFromPatch<typename HybridModelT::gridLayout_type>(*patch);

            auto writeMyField = [&layout, &patchPath, &hifile](auto const& field, std::string id) {
                using T     = typename FIELD::data_type;
                auto iStart = layout.physicalStartIndex(field, Direction::X);
                auto iEnd   = layout.physicalEndIndex(field, Direction::X);
                std::vector<T> data(iEnd - iStart + 1);
                size_t field_idx = 0;
                for (auto ix = iStart; ix <= iEnd; ++ix)
                {
                    data[field_idx++] = field(ix);
                }
                std::stringstream fieldPath;
                fieldPath << patchPath.str() << "/" << id;
                HighFive::DataSet dataset
                    = hifile.createDataSet<T>(fieldPath.str(), HighFive::DataSpace::From(data));
                dataset.write(data);
            };
            std::map<std::string, std::reference_wrapper<FIELD>> fields
                = {{{"Ex", hybridModel->state.electromag.E.getComponent(Component::X)},
                    {"Ey", hybridModel->state.electromag.E.getComponent(Component::Y)},
                    {"Ez", hybridModel->state.electromag.E.getComponent(Component::Z)},
                    {"Bx", hybridModel->state.electromag.B.getComponent(Component::X)},
                    {"By", hybridModel->state.electromag.B.getComponent(Component::Y)},
                    {"Bz", hybridModel->state.electromag.B.getComponent(Component::Z)}}};

            for (auto [id, field] : fields)
            {
                writeMyField(field.get(), id);
            }

            patch_idx++;
        }
    }
}

int main(int argc, char** argv)
{
    PHARE_test::SamraiLifeCycle samrai_lifecycle(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
