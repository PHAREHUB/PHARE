#ifndef PHARE_TEST_SIMULATOR_PER_TEST_H
#define PHARE_TEST_SIMULATOR_PER_TEST_H

#include "phare/include.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

// Having multiple PythonDataProvider per binary execution doesn't work so well

struct __attribute__((visibility("hidden"))) StaticIntepreter
{
    static std::shared_ptr<PHARE::initializer::PythonDataProvider> input;

    StaticIntepreter()
    {
        input = std::make_shared<PHARE::initializer::PythonDataProvider>();
        input->read();
    }

    static StaticIntepreter& INSTANCE()
    {
        static StaticIntepreter i;
        return i;
    }
};
std::shared_ptr<PHARE::initializer::PythonDataProvider> StaticIntepreter::input = 0;


template<size_t _dim>
struct HierarchyMaker
{
    HierarchyMaker(PHARE::initializer::PHAREDict& dict)
        : hierarchy{std::make_shared<PHARE::amr::DimHierarchy<_dim>>(dict)}
    {
    }
    std::shared_ptr<PHARE::amr::Hierarchy> hierarchy;
};



template<size_t _dim, size_t _interp, size_t _nbRefinePart = 2>
struct TestSimulator : public HierarchyMaker<_dim>,
                       public PHARE::Simulator<_dim, _interp, _nbRefinePart>
{
    static constexpr size_t dim          = _dim;
    static constexpr size_t interp       = _interp;
    static constexpr size_t nbRefinePart = _nbRefinePart;

    using Simulator  = PHARE::Simulator<dim, interp, nbRefinePart>;
    using PHARETypes = PHARE::PHARE_Types<dim, interp, nbRefinePart>;
    using Hierarchy  = PHARE::amr::Hierarchy;

    using HybridModel = typename PHARETypes::HybridModel_t;
    using MHDModel    = typename PHARETypes::MHDModel_t;

    using ModelView        = PHARE::diagnostic::ModelView<Hierarchy, HybridModel>;
    using DiagnosticWriter = PHARE::diagnostic::h5::Writer<ModelView>;

    std::unique_ptr<ModelView> modelView;
    std::unique_ptr<DiagnosticWriter> writer;
    std::unique_ptr<PHARE::diagnostic::DiagnosticsManager<DiagnosticWriter>> dMan;


    auto& dict()
    {
        StaticIntepreter::INSTANCE();
        auto st = SAMRAI::hier::VariableDatabase::getDatabase();
        return PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
    }

    TestSimulator()
        : HierarchyMaker<_dim>{dict()}
        , Simulator{dict(), this->hierarchy}
    {
        Simulator::initialize();

        if (dict()["simulation"].contains("diagnostics"))
        {
            modelView = std::make_unique<ModelView>(*this->hierarchy, *this->getHybridModel());
            writer    = DiagnosticWriter::from(*modelView, dict()["simulation"]["diagnostics"]);
            dMan      = PHARE::diagnostic::DiagnosticsManager<DiagnosticWriter>::from(
                *writer, dict()["simulation"]["diagnostics"]);
        }
    }

    template<typename DMan>
    void dump(DMan& dman, double current_timestamp = 0, double current_timestep = 1)
    {
        dman.dump(current_timestamp, current_timestep);
    }
};

template<typename Simulator>
struct SimulatorTest : public ::testing::Test
{
};

using Simulators = testing::Types<TestSimulator<1, 1>, TestSimulator<1, 2>, TestSimulator<1, 3>/*,   // dim 1
                                  TestSimulator<2, 1>, TestSimulator<2, 2>, TestSimulator<2, 3>,     // dim 2
                                  TestSimulator<3, 1>, TestSimulator<3, 2>, TestSimulator<3, 3>*/>;  // dim 3
TYPED_TEST_SUITE(SimulatorTest, Simulators);

/*
int main(int argc, char** argv)
{
    int testResult = RUN_ALL_TESTS();

    StaticIntepreter::INSTANCE().kill(); // <-- mandatory

    return testResult;
}
*/


namespace PHARE
{
class FieldNullFilter
{
public:
    template<typename Field, typename GridLayout>
    size_t start(GridLayout& layout, core::Direction direction)
    {
        return layout.ghostStartIndex(typename Field::physical_quantity_type{}, direction);
    }
    template<typename Field, typename GridLayout>
    size_t end(GridLayout& layout, core::Direction direction)
    {
        return layout.ghostEndIndex(typename Field::physical_quantity_type{}, direction);
    }
};

class FieldDomainFilter
{
public:
    template<typename Field, typename GridLayout>
    size_t start(GridLayout& layout, core::Direction direction)
    {
        return layout.physicalStartIndex(typename Field::physical_quantity_type{}, direction);
    }
    template<typename Field, typename GridLayout>
    size_t end(GridLayout& layout, core::Direction direction)
    {
        return layout.physicalStartIndex(typename Field::physical_quantity_type{}, direction);
    }
};

class FieldDomainPlusNFilter
{
public:
    FieldDomainPlusNFilter(size_t n)
        : n_{n}
    {
    }
    template<typename Field, typename GridLayout>
    size_t start(GridLayout& layout, core::Direction direction)
    {
        return layout.physicalStartIndex(typename Field::physical_quantity_type{}, direction) - n_;
    }
    template<typename Field, typename GridLayout>
    size_t end(GridLayout& layout, core::Direction direction)
    {
        return layout.physicalStartIndex(typename Field::physical_quantity_type{}, direction) + n_;
    }

private:
    size_t n_;
};

} // namespace PHARE

#endif /* PHARE_TEST_SIMULATOR_PER_TEST_H */
