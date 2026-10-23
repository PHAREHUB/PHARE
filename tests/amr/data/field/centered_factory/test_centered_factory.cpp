#include "phare_mpi.hpp"
#include "phare_core.hpp"
#include "simulator/simulator_def.hpp"

#include "core/data/grid/grid.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"

#include "amr/data/field/field_variable.hpp"
#include "amr/data/field/field_variable_fill_pattern.hpp"
#include "amr/data/tensorfield/tensor_field_variable.hpp"
#include "amr/messengers/refiner.hpp"

#include <SAMRAI/hier/VariableDatabase.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>

#include "gtest/gtest.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

using namespace PHARE::core;
using namespace PHARE::amr;

// field factories have a distinct type per centering, so that SAMRAI's equivalence classes
// (factory type + ghost width) never mix centerings: see centered_factory.hpp

namespace
{
constexpr std::size_t dim    = 2;
constexpr std::size_t interp = 1;

using GridLayout_t = PHARE_Types<PHARE::SimOpts{dim, interp}>::Hybrid::GridLayout_t;
using Grid_t       = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;

using ScalarVariable = FieldVariable<GridLayout_t, Grid_t>;
using VectorVariable = TensorFieldVariable<1, GridLayout_t, Grid_t, HybridQuantity>;


// SAMRAI variable names are global: each call gives a new one
std::string uniqueName(std::string const& name)
{
    static int counter = 0;
    return name + "_" + std::to_string(counter++);
}

template<typename Variable, typename Quantity>
int registerVariable(std::string const& name, Quantity qty)
{
    auto* db = SAMRAI::hier::VariableDatabase::getDatabase();
    return db->registerVariableAndContext(
        std::make_shared<Variable>(uniqueName(name), qty), db->getContext("centered_factory"),
        SAMRAI::hier::IntVector::getZero(SAMRAI::tbox::Dimension{dim}));
}

auto const& factoryOf(int const id)
{
    return *SAMRAI::hier::VariableDatabase::getDatabase()
                ->getPatchDescriptor()
                ->getPatchDataFactory(id);
}

struct Item
{
    int id;
    std::shared_ptr<SAMRAI::xfer::VariableFillPattern> fillPattern = nullptr;
};

SAMRAI::xfer::RefineClasses& classesOf(SAMRAI::xfer::RefineAlgorithm& algo,
                                       std::vector<Item> const& items)
{
    for (auto const& item : items)
        algo.registerRefine(item.id, item.id, item.id, nullptr, item.fillPattern);
    return *algo.getEquivalenceClasses();
}
} // namespace



TEST(CenteredFactory, differentCenteringsGetDifferentFactoryTypesAndClasses)
{
    auto const B         = registerVariable<VectorVariable>("B", HybridQuantity::Vector::B);
    auto const E         = registerVariable<VectorVariable>("E", HybridQuantity::Vector::E);
    auto const& BFactory = factoryOf(B);
    auto const& EFactory = factoryOf(E);
    EXPECT_NE(typeid(BFactory), typeid(EFactory));

    SAMRAI::xfer::RefineAlgorithm vectors;
    EXPECT_EQ(2, classesOf(vectors, {{B}, {E}}).getNumberOfEquivalenceClasses());

    auto const Bx = registerVariable<ScalarVariable>("Bx", HybridQuantity::Scalar::Bx);
    auto const By = registerVariable<ScalarVariable>("By", HybridQuantity::Scalar::By);

    SAMRAI::xfer::RefineAlgorithm scalars;
    EXPECT_EQ(2, classesOf(scalars, {{Bx}, {By}}).getNumberOfEquivalenceClasses());
}


TEST(CenteredFactory, sameCenteringSharesFactoryTypeAndClass)
{
    auto const E         = registerVariable<VectorVariable>("E", HybridQuantity::Vector::E);
    auto const J         = registerVariable<VectorVariable>("J", HybridQuantity::Vector::J);
    auto const& EFactory = factoryOf(E);
    auto const& JFactory = factoryOf(J);
    EXPECT_EQ(typeid(EFactory), typeid(JFactory));

    SAMRAI::xfer::RefineAlgorithm vectors;
    EXPECT_EQ(1, classesOf(vectors, {{E}, {J}}).getNumberOfEquivalenceClasses());

    auto const rho = registerVariable<ScalarVariable>("rho", HybridQuantity::Scalar::rho);
    auto const Vx  = registerVariable<ScalarVariable>("Vx", HybridQuantity::Scalar::Vx);

    SAMRAI::xfer::RefineAlgorithm scalars;
    EXPECT_EQ(1, classesOf(scalars, {{rho}, {Vx}}).getNumberOfEquivalenceClasses());
}


TEST(CenteredFactory, guardThrowsOnDifferentOverwriteInteriorInOneClass)
{
    auto const rho = registerVariable<ScalarVariable>("rho", HybridQuantity::Scalar::rho);
    auto const Vx  = registerVariable<ScalarVariable>("Vx", HybridQuantity::Scalar::Vx);

    auto const overwrite    = std::make_shared<FieldFillPattern<dim>>(true);
    auto const notOverwrite = std::make_shared<FieldFillPattern<dim>>(false);

    SAMRAI::xfer::RefineAlgorithm agreeing;
    auto& agreeingClasses = classesOf(agreeing, {{rho, notOverwrite}, {Vx, notOverwrite}});
    EXPECT_EQ(1, agreeingClasses.getNumberOfEquivalenceClasses());
    EXPECT_NO_THROW(checkOverwriteInteriorAgrees(agreeingClasses));

    // fill patterns of one type are not told apart by SAMRAI: Vx would get rho's overlaps
    SAMRAI::xfer::RefineAlgorithm mismatched;
    auto& mismatchedClasses = classesOf(mismatched, {{rho, notOverwrite}, {Vx, overwrite}});
    EXPECT_EQ(1, mismatchedClasses.getNumberOfEquivalenceClasses());
    EXPECT_THROW(checkOverwriteInteriorAgrees(mismatchedClasses), std::runtime_error);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();

    int testResult = RUN_ALL_TESTS();

    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();

    return testResult;
}
