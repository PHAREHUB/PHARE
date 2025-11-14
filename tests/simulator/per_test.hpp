#ifndef PHARE_TEST_SIMULATOR_PER_TEST_HPP
#define PHARE_TEST_SIMULATOR_PER_TEST_HPP

#include "phare/phare.hpp"
#include "initializer/python_data_provider.hpp"
#include "tests/core/data/field/test_field.hpp"
#include "python3/mhd_defaults/default_mhd_time_stepper.hpp"

#include "gtest/gtest.h"

using SimOpts = PHARE::SimOpts<>;

struct __attribute__((visibility("hidden"))) StaticIntepreter
{
    static inline std::shared_ptr<PHARE::initializer::PythonDataProvider> input{nullptr};

    static StaticIntepreter& INSTANCE()
    {
        static StaticIntepreter i;
        return i;
    }
};


template<std::size_t _dim>
struct HierarchyMaker
{
    HierarchyMaker(PHARE::initializer::PHAREDict& dict)
        : hierarchy{std::make_shared<PHARE::amr::DimHierarchy<_dim>>(dict)}
    {
    }
    std::shared_ptr<PHARE::amr::Hierarchy> hierarchy;
};



template<auto opts>
struct SimulatorTestParam : private HierarchyMaker<opts.dimension>, public PHARE::Simulator<opts>
{
    static constexpr std::size_t dim = opts.dimension;

    using Simulator   = PHARE::Simulator<opts>;
    using PHARETypes  = PHARE::PHARE_Types<opts>;
    using Hierarchy   = PHARE::amr::Hierarchy;
    using HybridModel = PHARETypes::HybridModel_t;
    using MHDModel    = PHARETypes::MHDModel_t;
    using HierarchyMaker<dim>::hierarchy;

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
        : HierarchyMaker<dim>{dict(job_py)}
        , Simulator{dict(job_py), this->hierarchy}
    {
        Simulator::initialize();
    }

    void reset_dman() { Simulator::reset_dman(); }
};


template<typename SimulatorParam>
struct SimulatorTest : public ::testing::Test
{
};

using Simulators = testing::Types<SimulatorTestParam<SimOpts{1, 1, 2}>>;
TYPED_TEST_SUITE(SimulatorTest, Simulators);


template<typename SimulatorParam>
struct Simulator1dTest : public ::testing::Test
{
};

// clang-format off
using Simulators1d = testing::Types<
    SimulatorTestParam<SimOpts{1, 1, 2}>, SimulatorTestParam<SimOpts{1, 1, 3}>,
    SimulatorTestParam<SimOpts{1, 2, 2}>, SimulatorTestParam<SimOpts{1, 2, 3}>,
    SimulatorTestParam<SimOpts{1, 2, 4}>, SimulatorTestParam<SimOpts{1, 3, 2}>,
    SimulatorTestParam<SimOpts{1, 3, 3}>, SimulatorTestParam<SimOpts{1, 3, 4}>,
    SimulatorTestParam<SimOpts{1, 3, 5}>
>;

TYPED_TEST_SUITE(Simulator1dTest, Simulators1d);


template<typename SimulatorParam>
struct Simulator2dTest : public ::testing::Test
{
};


using Simulators2d = testing::Types<
    SimulatorTestParam<SimOpts{2, 1, 4}>, SimulatorTestParam<SimOpts{2, 1, 5}>,
    SimulatorTestParam<SimOpts{2, 1, 8}>, SimulatorTestParam<SimOpts{2, 1, 9}>,
    SimulatorTestParam<SimOpts{2, 2, 4}>, SimulatorTestParam<SimOpts{2, 2, 5}>,
    SimulatorTestParam<SimOpts{2, 2, 8}>, SimulatorTestParam<SimOpts{2, 2, 9}>,
    SimulatorTestParam<SimOpts{2, 2, 16}>, SimulatorTestParam<SimOpts{2, 3, 4}>,
    SimulatorTestParam<SimOpts{2, 3, 5}>, SimulatorTestParam<SimOpts{2, 3, 8}>,
    SimulatorTestParam<SimOpts{2, 3, 25}>
>;

TYPED_TEST_SUITE(Simulator2dTest, Simulators2d);
// clang-format on



#endif /* PHARE_TEST_SIMULATOR_PER_TEST_H */
