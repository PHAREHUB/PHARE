#ifndef PHARE_TEST_SIMULATOR_PER_TEST_H
#define PHARE_TEST_SIMULATOR_PER_TEST_H

#include "phare/phare.h"
#include "initializer/python_data_provider.h"
#include "tests/core/data/field/test_field.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


struct __attribute__((visibility("hidden"))) StaticIntepreter
{
    static std::shared_ptr<PHARE::initializer::PythonDataProvider> input;

    static StaticIntepreter& INSTANCE()
    {
        static StaticIntepreter i;
        return i;
    }
};
std::shared_ptr<PHARE::initializer::PythonDataProvider> StaticIntepreter::input{nullptr};


template<std::size_t _dim>
struct HierarchyMaker
{
    HierarchyMaker(PHARE::initializer::PHAREDict& dict)
        : hierarchy{std::make_shared<PHARE::amr::DimHierarchy<_dim>>(dict)}
    {
    }
    std::shared_ptr<PHARE::amr::Hierarchy> hierarchy;
};



template<std::size_t _dim, std::size_t _interp, std::size_t _nbRefinePart>
struct SimulatorTestParam : private HierarchyMaker<_dim>,
                            public PHARE::Simulator<_dim, _interp, _nbRefinePart>
{
    static constexpr std::size_t dim          = _dim;
    static constexpr std::size_t interp       = _interp;
    static constexpr std::size_t nbRefinePart = _nbRefinePart;

    using Simulator   = PHARE::Simulator<dim, interp, nbRefinePart>;
    using PHARETypes  = PHARE::PHARE_Types<dim, interp, nbRefinePart>;
    using Hierarchy   = PHARE::amr::Hierarchy;
    using HybridModel = typename PHARETypes::HybridModel_t;
    using MHDModel    = typename PHARETypes::MHDModel_t;
    using HierarchyMaker<dim>::hierarchy;

    std::unique_ptr<PHARE::diagnostic::IDiagnosticsManager> dMan;

    auto& dict(std::string job_py)
    {
        auto& input = StaticIntepreter::INSTANCE().input;
        if (!input)
        {
            input = std::make_shared<PHARE::initializer::PythonDataProvider>(job_py);
            input->read();
        }
        SAMRAI::hier::VariableDatabase::getDatabase();
        return PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
    }

    SimulatorTestParam(std::string job_py = "job")
        : HierarchyMaker<_dim>{dict(job_py)}
        , Simulator{dict(job_py), this->hierarchy}
    {
        Simulator::initialize();

        if (dict(job_py)["simulation"].contains("diagnostics"))
        {
            dMan = PHARE::diagnostic::DiagnosticsManagerResolver::make_unique(
                *this->hierarchy, *this->getHybridModel(),
                dict(job_py)["simulation"]["diagnostics"]);
        }
    }


    template<typename DMan>
    void dump(DMan& dman, double current_timestamp = 0, double current_timestep = 1)
    {
        dman.dump(current_timestamp, current_timestep);
    }
};

template<typename SimulatorParam>
struct SimulatorTest : public ::testing::Test
{
};

using Simulators = testing::Types<SimulatorTestParam<1, 1, 2>>;
TYPED_TEST_SUITE(SimulatorTest, Simulators);


template<typename SimulatorParam>
struct Simulator1dTest : public ::testing::Test
{
};


using Simulators1d = testing::Types<
    SimulatorTestParam<1, 1, 2>, SimulatorTestParam<1, 1, 3>, SimulatorTestParam<1, 2, 2>,
    SimulatorTestParam<1, 2, 3>, SimulatorTestParam<1, 2, 4>, SimulatorTestParam<1, 3, 2>,
    SimulatorTestParam<1, 3, 3>, SimulatorTestParam<1, 3, 4>, SimulatorTestParam<1, 3, 5>>;
TYPED_TEST_SUITE(Simulator1dTest, Simulators1d);


template<typename SimulatorParam>
struct Simulator2dTest : public ::testing::Test
{
};


using Simulators2d = testing::Types<
    SimulatorTestParam<2, 1, 4>, SimulatorTestParam<2, 1, 5>, SimulatorTestParam<2, 1, 8>,
    SimulatorTestParam<2, 1, 9>, SimulatorTestParam<2, 2, 4>, SimulatorTestParam<2, 2, 5>,
    SimulatorTestParam<2, 2, 8>, SimulatorTestParam<2, 2, 9>, SimulatorTestParam<2, 2, 16>,
    SimulatorTestParam<2, 3, 4>, SimulatorTestParam<2, 3, 5>, SimulatorTestParam<2, 3, 8>,
    SimulatorTestParam<2, 3, 25>>;
TYPED_TEST_SUITE(Simulator2dTest, Simulators2d);




#endif /* PHARE_TEST_SIMULATOR_PER_TEST_H */
