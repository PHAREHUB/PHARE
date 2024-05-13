#ifndef PHARE_TEST_DIAGNOSTIC_INCLUDE_HPP
#define PHARE_TEST_DIAGNOSTIC_INCLUDE_HPP

#include "tests/simulator/per_test.hpp"

#include "diagnostic/diagnostic_model_view.hpp"
#include "diagnostic/detail/h5writer.hpp"
#include "diagnostic/detail/types/electromag.hpp"
#include "diagnostic/detail/types/particle.hpp"
#include "diagnostic/detail/types/fluid.hpp"

#include <functional>


using namespace PHARE;
using namespace PHARE::diagnostic;
using namespace PHARE::diagnostic::h5;

constexpr unsigned NEW_HI5_FILE
    = HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate;


template<typename GridLayout, typename Field, typename FieldFilter = PHARE::FieldNullFilter>
auto checkField(HighFiveFile const& hifile, GridLayout const& layout, Field const& field,
                std::string const path, FieldFilter const ff = FieldFilter{})
{
    constexpr auto dim = GridLayout::dimension;
    static_assert(dim >= 1 and dim <= 3, "Invalid dimension.");

    auto fieldV = hifile.read_data_set_flat<float, dim>(path);
    PHARE::core::test(layout, field, fieldV, ff);
    return fieldV; // possibly unused
}

template<typename GridLayout, typename VecField, typename FieldFilter = PHARE::FieldNullFilter>
void checkVecField(HighFiveFile const& file, GridLayout const& layout, VecField const& vecField,
                   std::string const path, FieldFilter const ff = FieldFilter{})
{
    for (auto& [id, type] : core::Components::componentMap())
        checkField(file, layout, vecField.getComponent(type), path + "_" + id, ff);
}


template<typename Hierarchy, typename HybridModel>
struct Hi5Diagnostic
{
    using ModelView_t = ModelView<Hierarchy, HybridModel>;
    using Writer_t    = H5Writer<ModelView_t>;

    Hi5Diagnostic(Hierarchy& hierarchy, HybridModel& hybridModel, std::string out,
                  unsigned flags = NEW_HI5_FILE)
        : hierarchy_{hierarchy}
        , model_{hybridModel}
        , out_{out}
        , flags_{flags}
        , dMan{std::make_unique<Writer_t>(hierarchy_, model_, out, flags_)}
        , writer{dMan.writer()}
        , modelView{writer.modelView()}
    {
    }
    ~Hi5Diagnostic() {}

    auto dict(std::string&& type, std::string& quantity)
    {
        PHARE::initializer::PHAREDict dict;
        dict["name"]        = quantity;
        dict["type"]        = type;
        dict["quantity"]    = quantity;
        dict["time_step"]   = double{1};
        dict["flush_every"] = Writer_t::flush_never;

        dict["write_timestamps"]   = std::vector<double>{0, 1, 2};
        dict["compute_timestamps"] = std::vector<double>{0, 1, 2};
        dict["n_attributes"]       = std::size_t{0}; // real diag attrs are loaded from python

        return dict;
    }
    auto electromag(std::string&& type) { return dict("electromag", type); }
    auto particles(std::string&& type) { return dict("particle", type); }
    auto fluid(std::string&& type) { return dict("fluid", type); }

    // timestamp is constant precision of 10 places
    std::string getPatchPath(int level, std::string patch, std::string timestamp = "0.0000000000")
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
    using GridLayout = typename Simulator::PHARETypes::GridLayout_t;

    auto& hybridModel = *sim.getHybridModel();

    auto checkF = [&](auto& layout, auto& path, auto tree, auto name, auto& field) {
        auto hifile = hi5.writer.makeFile(hi5.writer.fileString(tree + name), hi5.flags_);
        auto&& data = checkField(*hifile, layout, field, path + name, FieldDomainFilter{});
    };

    auto checkVF = [&](auto& layout, auto& path, auto tree, auto name, auto& val) {
        auto hifile = hi5.writer.makeFile(hi5.writer.fileString(tree + name), hi5.flags_);
        checkVecField(*hifile, layout, val, path + name, FieldDomainFilter{});
    };

    auto visit = [&](GridLayout& layout, std::string patchID, std::size_t iLevel) {
        auto path  = hi5.getPatchPath(iLevel, patchID);
        auto& ions = hi5.modelView.getIons();
        for (auto& pop : ions)
        {
            checkF(layout, path, "/ions/pop/" + pop.name(), "/density"s, pop.density());
            checkVF(layout, path, "/ions/pop/" + pop.name(), "/flux"s, pop.flux());
        }
        checkF(layout, path, "/ions"s, "/density"s, ions.density());

        std::string tree{"/ions"}, var{"/bulkVelocity"};
        auto hifile = hi5.writer.makeFile(hi5.writer.fileString(tree + var), hi5.flags_);
        checkVecField(*hifile, layout, ions.velocity(), path + var, FieldDomainFilter{});
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
        auto hifile = hi5.writer.makeFile(hi5.writer.fileString(tree), hi5.flags_);
        checkVecField(*hifile, layout, val, path + tree);
    };

    auto visit = [&](GridLayout& layout, std::string patchID, std::size_t iLevel) {
        auto path = hi5.getPatchPath(iLevel, patchID) + "/";
        checkVF(layout, path, "EM_B", hybridModel.state.electromag.B);
        checkVF(layout, path, "EM_E", hybridModel.state.electromag.E);
    };

    PHARE::amr::visitHierarchy<GridLayout>(*sim.hierarchy, *hybridModel.resourcesManager, visit, 0,
                                           sim.hierarchy->getNumberOfLevels(), hybridModel);
}


template<typename Simulator, typename Hi5Diagnostic>
void validateParticleDump(Simulator& sim, Hi5Diagnostic& hi5)
{
    using GridLayout = typename Simulator::PHARETypes::GridLayout_t;

    auto& hybridModel = *sim.getHybridModel();

    auto checkParticles = [&](auto& hifile, auto& particles, auto path) {
        if (!particles.size())
            return;
        auto weightV = hifile.template read_data_set_flat<float, 2>(path + "weight");
        auto chargeV = hifile.template read_data_set_flat<float, 2>(path + "charge");
        auto vV      = hifile.template read_data_set_flat<float, 2>(path + "v");
        auto iCellV  = hifile.template read_data_set_flat<float, 2>(path + "iCell");
        auto deltaV  = hifile.template read_data_set_flat<float, 2>(path + "delta");

        core::ParticlePacker packer{particles};

        auto first            = packer.empty();
        std::size_t iCellSize = std::get<2>(first).size();
        std::size_t deltaSize = std::get<3>(first).size();
        std::size_t vSize     = std::get<4>(first).size();
        std::size_t part_idx  = 0;
        while (packer.hasNext())
        {
            auto next = packer.next();

            for (std::size_t i = 0; i < iCellSize; i++)
                EXPECT_EQ(iCellV[(part_idx * iCellSize) + i], std::get<2>(next)[i]);

            for (std::size_t i = 0; i < deltaSize; i++)
                EXPECT_FLOAT_EQ(deltaV[(part_idx * deltaSize) + i], std::get<3>(next)[i]);

            for (std::size_t i = 0; i < vSize; i++)
                EXPECT_FLOAT_EQ(vV[(part_idx * vSize) + i], std::get<4>(next)[i]);

            part_idx++;
        }
    };

    auto checkFile = [&](auto& path, auto tree, auto& particles) {
        auto hifile = hi5.writer.makeFile(hi5.writer.fileString(tree), hi5.flags_);
        checkParticles(*hifile, particles, path + "/");
    };

    auto visit = [&](GridLayout&, std::string patchID, std::size_t iLevel) {
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
    using GridLayout                           = typename Simulator::PHARETypes::GridLayout_t;
    constexpr auto dimension                   = Simulator::dimension;
    constexpr std::size_t expectedPopNbr       = 2;
    constexpr std::size_t expectedPopAttrFiles = 5;

    std::string const ionsPopPath = "/ions/pop/";

    auto& hybridModel = *sim.getHybridModel();
    auto& dict        = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();

    auto nbrPop = dict["simulation"]["ions"]["nbrPopulations"].template to<std::size_t>();
    EXPECT_EQ(nbrPop, expectedPopNbr);

    std::vector<std::string> h5FileTypes{"/EM_B", "/EM_E", "/ions/density", "/ions/bulkVelocity"};

    for (std::size_t i = 0; i < nbrPop; ++i)
    {
        std::string popName = dict["simulation"]["ions"]["pop" + std::to_string(i)]["name"]
                                  .template to<std::string>();

        h5FileTypes.emplace_back(ionsPopPath + popName + "/domain");
        h5FileTypes.emplace_back(ionsPopPath + popName + "/levelGhost");
        h5FileTypes.emplace_back(ionsPopPath + popName + "/patchGhost");
        h5FileTypes.emplace_back(ionsPopPath + popName + "/density");
        h5FileTypes.emplace_back(ionsPopPath + popName + "/flux");
    }

    auto _check_equal = [](auto& group, auto expected, auto key) {
        std::vector<typename decltype(expected)::value_type> attr;
        group.getAttribute(key).read(attr);
        EXPECT_EQ(expected, attr);
    };

    std::vector<std::string> const boundaryTypes(dimension, "periodic");
    std::size_t popAttrChecks = 0;
    for (auto const& fileType : h5FileTypes)
    {
        auto hifile = hi5.writer.makeFile(hi5.writer.fileString(fileType), hi5.flags_);

        auto visit = [&](GridLayout& grid, std::string patchID, std::size_t iLevel) {
            auto group = hifile->file().getGroup(hi5.getPatchPath(iLevel, patchID));

            _check_equal(group, grid.origin().toVector(), "origin");
            _check_equal(group, core::Point<std::uint32_t, dimension>{grid.nbrCells()}.toVector(),
                         "nbrCells");
            _check_equal(group, grid.AMRBox().lower.toVector(), "lower");
            _check_equal(group, grid.AMRBox().upper.toVector(), "upper");
        };

        PHARE::amr::visitHierarchy<GridLayout>(*sim.hierarchy, *hybridModel.resourcesManager, visit,
                                               0, sim.hierarchy->getNumberOfLevels(), hybridModel);

        auto rootGroup = hifile->file().getGroup("/");
        if (fileType.find(ionsPopPath) == 0)
        {
            ++popAttrChecks;
            double pop_mass = 0;
            rootGroup.getAttribute("pop_mass").read(pop_mass);
            EXPECT_DOUBLE_EQ(pop_mass, 1.0);
            _check_equal(rootGroup, boundaryTypes, "boundary_conditions");
        }
    }
    EXPECT_EQ(popAttrChecks, expectedPopNbr * expectedPopAttrFiles);
}


#endif /*PHARE_TEST_DIAGNOSTIC_INCLUDE_H*/
