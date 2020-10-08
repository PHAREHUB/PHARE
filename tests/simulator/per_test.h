#ifndef PHARE_TEST_SIMULATOR_PER_TEST_H
#define PHARE_TEST_SIMULATOR_PER_TEST_H

#include "phare/phare.h"
#include "initializer/python_data_provider.h"

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
struct TestSimulator : public HierarchyMaker<_dim>,
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

    TestSimulator(std::string job_py = "job")
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

template<typename Simulator>
struct /*[[deprecated]]*/ SimulatorTest : public ::testing::Test
{
};

// deprecated
using Simulators
    = testing::Types<TestSimulator<1, 1, 2>, TestSimulator<1, 2, 2>, TestSimulator<1, 3, 2>>;
TYPED_TEST_SUITE(SimulatorTest, Simulators);


template<typename Simulator>
struct Simulator1dTest : public ::testing::Test
{
};
using Simulators1d
    = testing::Types<TestSimulator<1, 1, 2>, TestSimulator<1, 2, 2>, TestSimulator<1, 3, 2>>;
TYPED_TEST_SUITE(Simulator1dTest, Simulators1d);


template<typename Simulator>
struct Simulator2dTest : public ::testing::Test
{
};
using Simulators2d
    = testing::Types<TestSimulator<2, 1, 8>, TestSimulator<2, 2, 8>, TestSimulator<2, 3, 8>>;
TYPED_TEST_SUITE(Simulator2dTest, Simulators2d);



namespace PHARE
{
class FieldNullFilter
{
public:
    template<typename Field, typename GridLayout>
    std::size_t start(GridLayout const& layout, Field const& field, core::Direction const direction)
    {
        return layout.ghostStartIndex(field, direction);
    }

    template<typename Field, typename GridLayout>
    std::size_t end(GridLayout const& layout, Field const& field, core::Direction const direction)
    {
        return layout.ghostEndIndex(field, direction);
    }

    template<typename Field, typename GridLayout>
    std::size_t size(GridLayout const& layout, Field const& field, core::Direction const direction)
    {
        return end(layout, field, direction) - start(layout, field, direction) + 1;
    }
};

class FieldDomainPlusNFilter
{
public:
    FieldDomainPlusNFilter(std::size_t n = 0)
        : n_{n}
    {
    }

    template<typename Field, typename GridLayout>
    std::size_t start(GridLayout const& layout, Field const& field, core::Direction const direction)
    {
        return layout.physicalStartIndex(field, direction) - n_;
    }

    template<typename Field, typename GridLayout>
    std::size_t end(GridLayout const& layout, Field const& field, core::Direction const direction)
    {
        return layout.physicalEndIndex(field, direction) + n_;
    }

private:
    std::size_t n_;
};


struct FieldDomainFilter : public FieldDomainPlusNFilter
{
};


} // namespace PHARE

#endif /* PHARE_TEST_SIMULATOR_PER_TEST_H */
