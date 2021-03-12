#ifndef PHARE_TEST_DIAGNOSTIC_INCLUDE_H
#define PHARE_TEST_DIAGNOSTIC_INCLUDE_H

#include "tests/simulator/per_test.h"

#include "diagnostic/diagnostic_model_view.h"
#include "diagnostic/detail/h5writer.h"
#include "diagnostic/detail/types/electromag.h"
#include "diagnostic/detail/types/particle.h"
#include "diagnostic/detail/types/fluid.h"

#include <functional>


using namespace PHARE;
using namespace PHARE::diagnostic;
using namespace PHARE::diagnostic::h5;

constexpr unsigned NEW_HI5_FILE
    = HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate;

template<typename FieldFilter, typename Func, typename GridLayout, typename Field>
std::array<std::uint32_t, GridLayout::dimension> fieldIndices(FieldFilter ff, Func&& func,
                                                              GridLayout& layout, Field& field)
{
    constexpr auto dim = GridLayout::dimension;
    static_assert(dim >= 1 and dim <= 3, "Invalid dimension.");

    auto get = [&](auto dir) { return static_cast<std::uint32_t>(func(ff, layout, field, dir)); };

    if constexpr (dim == 1)
        return {get(core::Direction::X)};
    if constexpr (dim == 2)
        return {get(core::Direction::X), get(core::Direction::Y)};
    if constexpr (dim == 3)
        return {get(core::Direction::X), get(core::Direction::Y), get(core::Direction::Z)};
}


template<std::size_t dim, typename T, typename GridLayout, typename Field>
void validateFluidGhosts(std::vector<T> const& data, GridLayout const& layout, Field const& field)
{
    using Filter = FieldDomainPlusNFilter;
    auto filter  = Filter{1}; // include one ghost on each side

    auto beg = fieldIndices(filter, std::mem_fn(&Filter::template start<Field, GridLayout>), layout,
                            field);
    auto end = fieldIndices(filter, std::mem_fn(&Filter::template end<Field, GridLayout>), layout,
                            field);
    auto siz = fieldIndices(FieldNullFilter{},
                            std::mem_fn(&FieldNullFilter::template size<Field, GridLayout>), layout,
                            field);

    for (std::size_t d = 0; d < dim; d++)
    {
        ASSERT_TRUE(beg[d] > 0);
        ASSERT_TRUE(end[d] > beg[d]);
        ASSERT_TRUE(siz[d] > end[d]);
    }

    core::NdArrayView<dim, float> view{data, siz};

    {
        std::size_t nans = 0;
        for (auto const& v : data)
            if (std::isnan(v))
                nans++;

        auto phyStartIdx = layout.physicalStartIndex(field.physicalQuantity(), core::Direction::X);
        core::NdArrayMask mask{0, phyStartIdx - 2};
        EXPECT_EQ(mask.nCells(field), nans);
    }

    if constexpr (dim == 1)
    {
        EXPECT_FLOAT_EQ(data[beg[0]], data[beg[0] + 1]); // left
        EXPECT_FLOAT_EQ(data[end[0] - 1], data[end[0]]); // right
    }
    else if constexpr (dim == 2)
    {
        for (std::size_t j = beg[1] + 1; j < end[1]; j++)
        {
            EXPECT_FLOAT_EQ(view(beg[0], j), view(beg[0] + 1, j)); // bottom
            EXPECT_FLOAT_EQ(view(end[0], j), view(end[0] - 1, j)); // top
        }

        for (std::size_t i = beg[0] + 1; i < end[0]; i++)
        {
            EXPECT_FLOAT_EQ(view(i, beg[1]), view(i, beg[1] + 1)); // left
            EXPECT_FLOAT_EQ(view(i, end[1]), view(i, end[1] - 1)); // right
        }

        // bottom left
        EXPECT_FLOAT_EQ(view(beg[0], beg[1]), view(beg[0] + 1, beg[1] + 1));

        // bottom right
        EXPECT_FLOAT_EQ(view(beg[0], end[1]), view(beg[0] + 1, end[1] - 1));

        // top left square
        EXPECT_FLOAT_EQ(view(end[0], beg[1]), view(end[0] - 1, beg[1] + 1));

        // top right square
        EXPECT_FLOAT_EQ(view(end[0], end[1]), view(end[0] - 1, end[1] - 1));
    }
    else if constexpr (dim == 3)
    {
        throw std::runtime_error("not implemented");
    }
}


template<typename GridLayout, typename Field, typename FieldFilter = PHARE::FieldNullFilter>
auto checkField(HighFive::File const& file, GridLayout const& layout, Field const& field,
                std::string const path, FieldFilter const ff = FieldFilter{})
{
    constexpr auto dim = GridLayout::dimension;
    static_assert(dim >= 1 and dim <= 3, "Invalid dimension.");

    std::vector<float> fieldV;
    file.getDataSet(path).read(fieldV);
    EXPECT_EQ(fieldV.size(), field.size());

    auto siz = fieldIndices(FieldNullFilter{},
                            std::mem_fn(&FieldNullFilter::template size<Field, GridLayout>), layout,
                            field);

    std::size_t items
        = std::accumulate(std::begin(siz), std::end(siz), 1., std::multiplies<std::size_t>());
    EXPECT_EQ(items, field.size());

    auto beg = fieldIndices(ff, std::mem_fn(&FieldFilter::template start<Field, GridLayout>),
                            layout, field);
    auto end = fieldIndices(ff, std::mem_fn(&FieldFilter::template end<Field, GridLayout>), layout,
                            field);

    core::NdArrayView<dim, float> view{fieldV, siz};
    if constexpr (dim == 1)
    {
        for (std::size_t i = beg[0]; i < end[0]; i++)
        {
            if (std::isnan(view(i)) || std::isnan(field(i)))
                throw std::runtime_error("This field should not be NaN");
            EXPECT_FLOAT_EQ(view(i), field(i));
        }
    }
    else if constexpr (dim == 2)
    {
        for (std::size_t i = beg[0]; i < end[0]; i++)
        {
            for (std::size_t j = beg[1]; j < end[1]; j++)
            {
                if (std::isnan(view(i, j)) || std::isnan(field(i, j)))
                    throw std::runtime_error("This field should not be NaN");
                EXPECT_FLOAT_EQ(view(i, j), field(i, j));
            }
        }
    }
    else if constexpr (dim == 3)
    {
        for (std::size_t i = beg[0]; i < end[0]; i++)
        {
            for (std::size_t j = beg[1]; j < end[1]; j++)
            {
                for (std::size_t k = beg[2]; k < end[2]; k++)
                {
                    if (std::isnan(view(i, j, k)) || std::isnan(field(i, j, k)))
                        throw std::runtime_error("This field should not be NaN");
                    EXPECT_FLOAT_EQ(view(i, j, k), field(i, j, k));
                }
            }
        }
    }

    return fieldV; // possibly unused
}

template<typename GridLayout, typename VecField, typename FieldFilter = PHARE::FieldNullFilter>
void checkVecField(HighFive::File const& file, GridLayout const& layout, VecField const& vecField,
                   std::string const path, FieldFilter const ff = FieldFilter{})
{
    for (auto& [id, type] : core::Components::componentMap)
    {
        checkField(file, layout, vecField.getComponent(type), path + "_" + id, ff);
    }
}


template<typename Hierarchy, typename HybridModel>
struct Hi5Diagnostic
{
    using ModelView_t = ModelView<Hierarchy, HybridModel>;
    using Writer_t    = Writer<ModelView_t>;

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
        auto&& data = checkField(hifile->file(), layout, field, path + name, FieldDomainFilter{});
        /*
          Validate ghost of first border node is equal to border node
          see fixMomentGhosts() in src/core/data/ions/ions.h
        */
        validateFluidGhosts<GridLayout::dimension>(data, layout, field);
    };
    auto checkVF = [&](auto& layout, auto& path, auto tree, auto name, auto& val) {
        auto hifile = hi5.writer.makeFile(hi5.writer.fileString(tree + name), hi5.flags_);
        checkVecField(hifile->file(), layout, val, path + name, FieldDomainFilter{});
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
        auto hifile = hi5.writer.makeFile(hi5.writer.fileString(tree), hi5.flags_);
        checkVecField(hifile->file(), layout, val, path + tree);
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
        checkParticles(hifile->file(), particles, path + "/");
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

    int nbrPop = dict["simulation"]["ions"]["nbrPopulations"].template to<int>();
    EXPECT_EQ(nbrPop, expectedPopNbr);

    std::vector<std::string> h5FileTypes{"/EM_B", "/EM_E", "/ions/density", "/ions/bulkVelocity"};

    for (int i = 0; i < nbrPop; i++)
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
        }
    }
    EXPECT_EQ(popAttrChecks, expectedPopNbr * expectedPopAttrFiles);
}


#endif /*PHARE_TEST_DIAGNOSTIC_INCLUDE_H*/
