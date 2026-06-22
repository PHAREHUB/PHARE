


#include "simulator/simulator.hpp"
#include "amr/tagging/concrete_tagger.hpp"
#include "amr/tagging/tagging_criteria.hpp"

#include "tests/core/data/gridlayout/gridlayout_test.hpp"
#include "tests/core/data/vecfield/test_vecfield_fixtures.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cmath>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>

using namespace PHARE::amr;
using namespace PHARE::core;


// runtime dict contract: method + nbr_quantities + Q{i}/{name,threshold}
PHARE::initializer::PHAREDict taggingDict(std::string const& method,
                                          std::vector<std::pair<std::string, double>> const& qtys)
{
    PHARE::initializer::PHAREDict dict;
    dict["method"]         = method;
    dict["nbr_quantities"] = static_cast<int>(qtys.size());
    for (std::size_t i = 0; i < qtys.size(); ++i)
    {
        auto const path         = "Q" + std::to_string(i);
        dict[path]["name"]      = qtys[i].first;
        dict[path]["threshold"] = qtys[i].second;
    }
    return dict;
}


//-----------------------------------------------------------------------------
//  Construction (uses a real model type; ctor does not touch the model)
//-----------------------------------------------------------------------------

TEST(test_tagger, constructsWithValidMethodAndQuantity)
{
    auto static constexpr opts = PHARE::SimOpts{1ul, 1ul, 2ul};
    using phare_types          = PHARE::solver::PHARE_Types<opts>;
    using hybrid_model         = phare_types::HybridModel_t;
    auto dict                  = taggingDict("default", {{"B", 0.2}});
    EXPECT_NO_THROW((ConcreteTagger<hybrid_model>{dict}));
    auto dictLohner = taggingDict("lohner", {{"B", 0.2}});
    EXPECT_NO_THROW((ConcreteTagger<hybrid_model>{dictLohner}));
}

TEST(test_tagger, throwsOnUnknownMethod)
{
    auto static constexpr opts = PHARE::SimOpts{1ul, 1ul, 2ul};
    using phare_types          = PHARE::solver::PHARE_Types<opts>;
    using hybrid_model         = phare_types::HybridModel_t;
    auto dict                  = taggingDict("invalidStrat", {{"B", 0.2}});
    EXPECT_THROW((ConcreteTagger<hybrid_model>{dict}), std::runtime_error);
}


//-----------------------------------------------------------------------------
//  Pure criteria (tagging_criteria.hpp)
//-----------------------------------------------------------------------------

// minimal 1D field for testing the pure criteria functions directly
struct MockField1D
{
    std::vector<double> data;
    double operator()(std::uint32_t i) const { return data[i]; }
};

TEST(test_criteria, lohnerSeparatesStepFromFlat)
{
    // tanh step centered at cell 20 over 40 cells
    std::size_t const n = 40;
    MockField1D f;
    f.data.resize(n);
    for (std::size_t i = 0; i < n; ++i)
        f.data[i] = std::tanh((static_cast<double>(i) - 20.) / 1.5);

    std::vector<MockField1D const*> comps{&f};

    // Sample a flank of the step (cell 18), NOT the center: lohner is a second-
    // difference estimator and vanishes at the inflection point of a symmetric
    // tanh (a_p = -a_m, a_0 = 0). Curvature, and thus the indicator, is maximal
    // on the flanks.
    auto const atStep = lohnerIndicator<1>(comps, {18u});
    auto const atFlat = lohnerIndicator<1>(comps, {5u});
    EXPECT_GT(atStep, atFlat);
    EXPECT_GT(atStep, 10. * atFlat);
}

TEST(test_criteria, defaultSeparatesStepFromFlat)
{
    std::size_t const n = 40;
    MockField1D f;
    f.data.resize(n);
    for (std::size_t i = 0; i < n; ++i)
        f.data[i] = std::tanh((static_cast<double>(i) - 20.) / 1.5);

    std::vector<MockField1D const*> comps{&f};

    auto const atStep = defaultIndicator<1>(comps, {20u});
    auto const atFlat = defaultIndicator<1>(comps, {5u});
    EXPECT_GT(atStep, atFlat);
}

TEST(test_criteria, parseTaggingMethod)
{
    EXPECT_EQ(parseTaggingMethod("default"), TaggingMethod::Default);
    EXPECT_EQ(parseTaggingMethod("lohner"), TaggingMethod::Lohner);
    EXPECT_THROW(parseTaggingMethod("nope"), std::runtime_error);
}

// A field PRIMAL in x sampled at the dual cell index sits half a cell off, which
// breaks the symmetry of the centered stencil. CellCenteredSampler projects it onto
// the cell center (averaging bracketing nodes) and restores symmetry. Here a tent
// peaked on node 10 is symmetric about cell center 9.5: the projected indicator must
// satisfy ind(9)==ind(10); the raw co-indexed read does NOT (that is the bug).
TEST(test_criteria, cellCenteredProjectionRestoresSymmetry)
{
    MockField1D f;
    f.data.resize(21);
    for (std::size_t i = 0; i < f.data.size(); ++i)
        f.data[i] = -std::abs(static_cast<double>(i) - 10.0); // tent peaked at node 10

    // primal-in-x sampler -> projects to cell center; dual sampler -> raw passthrough
    CellCenteredSampler<1, MockField1D> primal{f, {true}};
    CellCenteredSampler<1, MockField1D> dual{f, {false}};
    std::vector<CellCenteredSampler<1, MockField1D> const*> projected{&primal};
    std::vector<CellCenteredSampler<1, MockField1D> const*> raw{&dual};

    auto const dProj9  = defaultIndicator<1>(projected, {9u});
    auto const dProj10 = defaultIndicator<1>(projected, {10u});
    EXPECT_NEAR(dProj9, dProj10, 1e-12); // projected: symmetric about cell 9.5

    auto const lProj9  = lohnerIndicator<1>(projected, {9u});
    auto const lProj10 = lohnerIndicator<1>(projected, {10u});
    EXPECT_NEAR(lProj9, lProj10, 1e-12);

    // raw co-indexing is asymmetric (regression guard: the projection is what fixes it)
    auto const dRaw9  = defaultIndicator<1>(raw, {9u});
    auto const dRaw10 = defaultIndicator<1>(raw, {10u});
    EXPECT_GT(std::abs(dRaw9 - dRaw10), 0.1);
}


//-----------------------------------------------------------------------------
//  Per-quantity field selection (generic resource-tree resolver) + union
//-----------------------------------------------------------------------------

namespace
{
constexpr std::size_t ib_tagger_dim    = 2;
constexpr std::size_t ib_tagger_interp = 1;
using IBTaggerGridLayout = GridLayout<GridLayoutImplYee<ib_tagger_dim, ib_tagger_interp>>;
} // namespace

namespace
{
struct TagFieldsMockState
{
    UsableVecField<ib_tagger_dim>& B;
    UsableVecField<ib_tagger_dim>& E;
    auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(B, E); }
    auto getCompileTimeResourcesViewList() const { return std::forward_as_tuple(B, E); }
};

struct TagFieldsMockModel
{
    using gridlayout_type                  = IBTaggerGridLayout;
    static constexpr std::size_t dimension = ib_tagger_dim;

    TagFieldsMockState state;
    UsableVecField<ib_tagger_dim>* B_ptr = nullptr;

    auto& get_B() { return *B_ptr; }
};

// fill component (ix,iy) with `ix` -> gradient large enough to tag every cell at threshold=0.1
template<typename Field, typename Layout, typename Qty>
void fillRamp(Field& f, Layout const& layout, Qty qty)
{
    auto const alloc = layout.allocSize(qty);
    for (std::size_t ix = 0; ix < alloc[0]; ++ix)
        for (std::size_t iy = 0; iy < alloc[1]; ++iy)
            f(ix, iy) = static_cast<double>(ix);
}

template<typename Field, typename Layout, typename Qty>
void fillZero(Field& f, Layout const& layout, Qty qty)
{
    auto const alloc = layout.allocSize(qty);
    for (std::size_t ix = 0; ix < alloc[0]; ++ix)
        for (std::size_t iy = 0; iy < alloc[1]; ++iy)
            f(ix, iy) = 0.0;
}
} // namespace


TEST(TagFields, EquantityTagsOnEOnly)
{
    constexpr std::size_t dim = ib_tagger_dim;
    auto const layout         = TestGridLayout<IBTaggerGridLayout>::make(20);

    // B identically zero, E_y a ramp -> "B" tags nothing, "E" tags everything.
    UsableVecField<dim> B{"B", layout, HybridQuantity::Vector::B};
    UsableVecField<dim> E{"E", layout, HybridQuantity::Vector::E};
    fillZero(B.getComponent(PHARE::core::Component::Y), layout, HybridQuantity::Scalar::By);
    fillRamp(E.getComponent(PHARE::core::Component::Y), layout, HybridQuantity::Scalar::Ey);

    TagFieldsMockModel model{TagFieldsMockState{B, E}, &B};
    auto const ncells = layout.nbrCells();

    std::vector<int> tagsB(ncells[0] * ncells[1], 0);
    ConcreteTaggerKernel<TagFieldsMockModel>{taggingDict("default", {{"B", 0.1}})}.tagFields(
        model, layout, tagsB.data());

    std::vector<int> tagsE(ncells[0] * ncells[1], 0);
    ConcreteTaggerKernel<TagFieldsMockModel>{taggingDict("default", {{"E", 0.1}})}.tagFields(
        model, layout, tagsE.data());

    bool constexpr fortran = false;
    auto tagsBv            = NdArrayView<dim, int, fortran>(tagsB.data(), ncells);
    auto tagsEv            = NdArrayView<dim, int, fortran>(tagsE.data(), ncells);
    for (auto const& p : boxFromNbrCells(ncells))
    {
        EXPECT_EQ(tagsBv(p.toArray()), 0) << "B should be 0 with zero B";
        EXPECT_EQ(tagsEv(p.toArray()), 1) << "E ramp should tag";
    }
}


TEST(TagFields, BxCompactNameSelectsSingleComponent)
{
    constexpr std::size_t dim = ib_tagger_dim;
    auto const layout         = TestGridLayout<IBTaggerGridLayout>::make(20);

    UsableVecField<dim> B{"B", layout, HybridQuantity::Vector::B};
    UsableVecField<dim> E{"E", layout, HybridQuantity::Vector::E};
    // Bx has a ramp, By/Bz zero -> "Bx" tags everywhere, "By" tags nothing.
    fillRamp(B.getComponent(PHARE::core::Component::X), layout, HybridQuantity::Scalar::Bx);

    TagFieldsMockModel model{TagFieldsMockState{B, E}, &B};
    auto const ncells = layout.nbrCells();

    std::vector<int> tagsBx(ncells[0] * ncells[1], 0);
    ConcreteTaggerKernel<TagFieldsMockModel>{taggingDict("default", {{"Bx", 0.1}})}.tagFields(
        model, layout, tagsBx.data());

    std::vector<int> tagsBy(ncells[0] * ncells[1], 0);
    ConcreteTaggerKernel<TagFieldsMockModel>{taggingDict("default", {{"By", 0.1}})}.tagFields(
        model, layout, tagsBy.data());

    bool constexpr fortran = false;
    auto tagsBxv           = NdArrayView<dim, int, fortran>(tagsBx.data(), ncells);
    auto tagsByv           = NdArrayView<dim, int, fortran>(tagsBy.data(), ncells);
    for (auto const& p : boxFromNbrCells(ncells))
    {
        EXPECT_EQ(tagsBxv(p.toArray()), 1);
        EXPECT_EQ(tagsByv(p.toArray()), 0);
    }
}


TEST(TagFields, UnionTagsIfAnyQuantityExceeds)
{
    constexpr std::size_t dim = ib_tagger_dim;
    auto const layout         = TestGridLayout<IBTaggerGridLayout>::make(20);

    // B zero, E ramp. {("B",0.1)} alone tags nothing; {("B",0.1),("E",0.1)} tags everywhere
    // because the union accepts a cell if ANY quantity's indicator exceeds its threshold.
    UsableVecField<dim> B{"B", layout, HybridQuantity::Vector::B};
    UsableVecField<dim> E{"E", layout, HybridQuantity::Vector::E};
    fillZero(B.getComponent(PHARE::core::Component::Y), layout, HybridQuantity::Scalar::By);
    fillRamp(E.getComponent(PHARE::core::Component::Y), layout, HybridQuantity::Scalar::Ey);

    TagFieldsMockModel model{TagFieldsMockState{B, E}, &B};
    auto const ncells = layout.nbrCells();

    std::vector<int> tags(ncells[0] * ncells[1], 0);
    ConcreteTaggerKernel<TagFieldsMockModel>{taggingDict("default", {{"B", 0.1}, {"E", 0.1}})}.tagFields(
        model, layout, tags.data());

    bool constexpr fortran = false;
    auto tagsv             = NdArrayView<dim, int, fortran>(tags.data(), ncells);
    for (auto const& p : boxFromNbrCells(ncells))
        EXPECT_EQ(tagsv(p.toArray()), 1) << "union should tag via E even though B is flat";
}


TEST(TagFields, UnknownQuantityNameThrows)
{
    constexpr std::size_t dim = ib_tagger_dim;
    auto const layout         = TestGridLayout<IBTaggerGridLayout>::make(20);

    UsableVecField<dim> B{"B", layout, HybridQuantity::Vector::B};
    UsableVecField<dim> E{"E", layout, HybridQuantity::Vector::E};
    TagFieldsMockModel model{TagFieldsMockState{B, E}, &B};

    ConcreteTaggerKernel<TagFieldsMockModel> tagger{
        taggingDict("default", {{"bogus_does_not_exist", 0.1}})};
    auto const ncells = layout.nbrCells();
    std::vector<int> tags(ncells[0] * ncells[1], 0);

    EXPECT_THROW(tagger.tagFields(model, layout, tags.data()), std::runtime_error);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    int testResult = RUN_ALL_TESTS();
    return testResult;
}
