#ifndef PHARE_TEST_DIAGNOSTIC_INCLUDE_H
#define PHARE_TEST_DIAGNOSTIC_INCLUDE_H

#include "tests/simulator/per_test.h"

#include "diagnostic/diagnostic_model_view.h"
#include "diagnostic/detail/h5writer.h"
#include "diagnostic/detail/types/electromag.h"
#include "diagnostic/detail/types/particle.h"
#include "diagnostic/detail/types/fluid.h"

using namespace PHARE;
using namespace PHARE::diagnostic;
using namespace PHARE::diagnostic::h5;

constexpr unsigned NEW_HI5_FILE
    = HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate;

template<typename FieldFilter, typename F, typename GridLayout>
std::array<uint32_t, GridLayout::dimension> fieldIndices(FieldFilter ff, F&& func,
                                                         GridLayout& layout)
{
    constexpr auto dim = GridLayout::dimension;
    static_assert(dim >= 1 and dim <= 3, "Invalid dimension.");

    auto get = [&](auto dir) { return static_cast<uint32_t>(((ff).*(func))(layout, dir)); };

    if constexpr (dim == 1)
        return {get(PHARE::core::Direction::X)};
    if constexpr (dim == 2)
        return {get(PHARE::core::Direction::X), get(PHARE::core::Direction::Y)};
    if constexpr (dim == 3)
        return {get(PHARE::core::Direction::X), get(PHARE::core::Direction::Y),
                get(PHARE::core::Direction::Z)};
}

template<typename GridLayout, typename Field, typename FieldFilter = PHARE::FieldNullFilter>
void checkField(HighFive::File& file, GridLayout& layout, Field& field, std::string path,
                FieldFilter ff = FieldFilter{})
{
    constexpr auto dim = GridLayout::dimension;
    static_assert(dim >= 1 and dim <= 3, "Invalid dimension.");

    std::vector<float> fieldV;
    file.getDataSet(path).read(fieldV);
    EXPECT_EQ(fieldV.size(), field.size());

    auto siz = fieldIndices(FieldNullFilter{}, &FieldNullFilter::template end<Field, GridLayout>,
                            layout);

    auto beg = fieldIndices(ff, &FieldFilter::template start<Field, GridLayout>, layout);
    auto end = fieldIndices(ff, &FieldFilter::template end<Field, GridLayout>, layout);

    if constexpr (dim == 1)
    {
        PHARE::core::NdArrayView<1, float> view{fieldV, siz};
        for (size_t i = beg[0]; i < end[0]; i++)
        {
            if (std::isnan(view(i)) || std::isnan(field(i)))
                throw std::runtime_error("This field should not be NaN");
            EXPECT_FLOAT_EQ(view(i), field(i));
        }
    }
    else if constexpr (dim == 2)
    {
        PHARE::core::NdArrayView<2, float> view{fieldV, siz};
        for (size_t i = beg[0]; i < end[0]; i++)
        {
            for (size_t j = beg[1]; j < end[1]; j++)
            {
                if (std::isnan(view(i, j)) || std::isnan(field(i, j)))
                    throw std::runtime_error("This field should not be NaN");
                EXPECT_FLOAT_EQ(view(i, j), field(i, j));
            }
        }
    }
    else
        throw std::runtime_error("Unhandled dimension");
}

template<typename GridLayout, typename VecField, typename FieldFilter = PHARE::FieldNullFilter>
void checkVecField(HighFive::File& file, GridLayout& layout, VecField& vecField, std::string path,
                   FieldFilter ff = FieldFilter{})
{
    for (auto& [id, type] : PHARE::core::Components::componentMap)
    {
        checkField(file, layout, vecField.getComponent(type), path + "_" + id, ff);
    }
}


template<typename Hierarchy, typename HybridModel>
struct Hi5Diagnostic
{
    using ModelView_t = ModelView<Hierarchy, HybridModel>;
    using Writer_t    = Writer<ModelView_t>;

    Hi5Diagnostic(Hierarchy& hierarchy, HybridModel& hybridModel, std::string out, unsigned flags)
        : hierarchy_{hierarchy}
        , model_{hybridModel}
        , out_{out}
        , flags_{flags}
        , dMan{std::make_unique<Writer_t>(hierarchy_, model_, "phare_outputs", flags_)}
        , writer{dMan.writer()}
        , modelView{writer.modelView()}
    {
    }
    ~Hi5Diagnostic() {}

    auto dict(std::string&& type, std::string& quantity)
    {
        PHARE::initializer::PHAREDict dict;
        dict["name"]      = quantity;
        dict["type"]      = type;
        dict["quantity"]  = quantity;
        dict["time_step"] = double{1};

        dict["write_timestamps"]   = std::vector<double>{0, 1, 2};
        dict["compute_timestamps"] = std::vector<double>{0, 1, 2};

        return dict;
    }
    auto electromag(std::string&& type) { return dict("electromag", type); }
    auto particles(std::string&& type) { return dict("particle", type); }
    auto fluid(std::string&& type) { return dict("fluid", type); }

    std::string getPatchPath(int level, std::string patch, std::string timestamp = "0.000000")
    {
        return Writer_t::getFullPatchPath(timestamp, level, patch);
    }

    Hierarchy& hierarchy_;
    HybridModel& model_;
    std::string out_;
    unsigned flags_;

    DiagnosticsManager<Writer_t> dMan;
    Writer_t& writer;
    ModelView_t const& modelView;
};



template<typename Simulator, typename Hi5Diagnostic>
void validateFluidDump(Simulator& sim, Hi5Diagnostic& hi5)
{
    using namespace std::string_literals;

    using GridLayout  = typename Simulator::PHARETypes::GridLayout_t;
    auto& hybridModel = *sim.getHybridModel();

    auto checkF = [&](auto& layout, auto& path, auto tree, auto name, auto& val) {
        auto hifile = hi5.writer.makeFile(hi5.writer.fileString(tree + name));
        checkField(hifile->file(), layout, val, path + name, FieldDomainFilter{});
    };
    auto checkVF = [&](auto& layout, auto& path, auto tree, auto name, auto& val) {
        auto hifile = hi5.writer.makeFile(hi5.writer.fileString(tree + name));
        checkVecField(hifile->file(), layout, val, path + name, FieldDomainFilter{});
    };

    auto visit = [&](GridLayout& layout, std::string patchID, size_t iLevel) {
        auto path  = hi5.getPatchPath(iLevel, patchID);
        auto& ions = hi5.modelView.getIons();
        for (auto& pop : ions)
        {
            checkF(layout, path, "/ions/pop/" + pop.name(), "/density"s, pop.density());
            checkVF(layout, path, "/ions/pop/" + pop.name(), "/flux"s, pop.flux());
        }
        checkF(layout, path, "/ions"s, "/density"s, ions.density());

        std::string tree{"/ions"}, var{"/bulkVelocity"};
        auto hifile = hi5.writer.makeFile(hi5.writer.fileString(tree + var));
        checkVecField(hifile->file(), layout, ions.velocity(), path + var, FieldDomainFilter{});
    };

    PHARE::amr::visitHierarchy<GridLayout>(*sim.hierarchy, *hybridModel.resourcesManager, visit, 0,
                                           sim.hierarchy->getNumberOfLevels(), hybridModel);
}



template<typename Simulator, typename Hi5Diagnostic>
void validateElectromagDump(Simulator& sim, Hi5Diagnostic& hi5)
{
    using GridLayout = typename Simulator::PHARETypes::GridLayout_t;

    auto& hybridModel = *sim.getHybridModel();

    auto checkVF = [&](auto& layout, auto& path, auto tree, auto& val) {
        auto hifile = hi5.writer.makeFile(hi5.writer.fileString(tree));
        checkVecField(hifile->file(), layout, val, path + tree);
    };

    auto visit = [&](GridLayout& layout, std::string patchID, size_t iLevel) {
        auto path = hi5.getPatchPath(iLevel, patchID) + "/";
        checkVF(layout, path, "/EM_B", hybridModel.state.electromag.B);
        checkVF(layout, path, "/EM_E", hybridModel.state.electromag.E);
    };

    PHARE::amr::visitHierarchy<GridLayout>(*sim.hierarchy, *hybridModel.resourcesManager, visit, 0,
                                           sim.hierarchy->getNumberOfLevels(), hybridModel);
}


template<typename Simulator, typename Hi5Diagnostic>
void validateParticleDump(Simulator& sim, Hi5Diagnostic& hi5)
{
    using GridLayout = typename Simulator::PHARETypes::GridLayout_t;

    auto& hybridModel = *sim.getHybridModel();

    auto checkParticles = [&](auto& file, auto& particles, auto path) {
        if (!particles.size())
            return;
        std::vector<float> weightV, chargeV, vV;
        file.getDataSet(path + "weight").read(weightV);
        file.getDataSet(path + "charge").read(chargeV);
        file.getDataSet(path + "v").read(vV);
        std::vector<int> iCellV;
        file.getDataSet(path + "iCell").read(iCellV);
        std::vector<float> deltaV;
        file.getDataSet(path + "delta").read(deltaV);

        core::ParticlePacker packer{particles};

        auto first       = packer.empty();
        size_t iCellSize = std::get<2>(first).size();
        size_t deltaSize = std::get<3>(first).size();
        size_t vSize     = std::get<4>(first).size();
        size_t part_idx  = 0;
        while (packer.hasNext())
        {
            auto next = packer.next();

            for (size_t i = 0; i < iCellSize; i++)
                EXPECT_EQ(iCellV[(part_idx * iCellSize) + i], std::get<2>(next)[i]);

            for (size_t i = 0; i < deltaSize; i++)
                EXPECT_FLOAT_EQ(deltaV[(part_idx * deltaSize) + i], std::get<3>(next)[i]);

            for (size_t i = 0; i < vSize; i++)
                EXPECT_FLOAT_EQ(vV[(part_idx * vSize) + i], std::get<4>(next)[i]);

            part_idx++;
        }
    };

    auto checkFile = [&](auto& path, auto tree, auto& particles) {
        auto hifile = hi5.writer.makeFile(hi5.writer.fileString(tree));
        checkParticles(hifile->file(), particles, path + "/");
    };

    auto visit = [&](GridLayout&, std::string patchID, size_t iLevel) {
        auto path = hi5.getPatchPath(iLevel, patchID);
        for (auto& pop : hybridModel.state.ions)
        {
            checkFile(path, "/ions/pop/" + pop.name() + "/domain", pop.domainParticles());
            checkFile(path, "/ions/pop/" + pop.name() + "/levelGhost", pop.levelGhostParticles());
        }
    };

    PHARE::amr::visitHierarchy<GridLayout>(*sim.hierarchy, *hybridModel.resourcesManager, visit, 0,
                                           sim.hierarchy->getNumberOfLevels(), hybridModel);
}


template<typename Simulator, typename Hi5Diagnostic>
void validateAttributes(Simulator& sim, Hi5Diagnostic& hi5)
{
    using GridLayout = typename Simulator::PHARETypes::GridLayout_t;

    auto& hybridModel = *sim.getHybridModel();
    auto hifile       = hi5.writer.makeFile(hi5.writer.fileString("/EM_B"));

    auto visit = [&](GridLayout& gridLayout, std::string patchID, size_t iLevel) {
        std::string patchPath = hi5.getPatchPath(iLevel, patchID), origin;
        hifile->file().getGroup(patchPath).getAttribute("origin").read(origin);
        EXPECT_EQ(gridLayout.origin().str(), origin);
    };

    PHARE::amr::visitHierarchy<GridLayout>(*sim.hierarchy, *hybridModel.resourcesManager, visit, 0,
                                           sim.hierarchy->getNumberOfLevels(), hybridModel);
}


#endif /*PHARE_TEST_DIAGNOSTIC_INCLUDE_H*/
