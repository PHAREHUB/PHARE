#ifndef PHARE_TEST_SIMULATOR_PER_TEST_H
#define PHARE_TEST_SIMULATOR_PER_TEST_H

#include "simulator/simulator.h"

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


template<size_t _dim, size_t _interp>
struct TestSimulator : public PHARE::Simulator<_dim, _interp>
{
    static constexpr size_t dim    = _dim;
    static constexpr size_t interp = _interp;

    using Simulator  = PHARE::Simulator<dim, interp>;
    using PHARETypes = PHARE::PHARE_Types<dim, interp>;
    using Hierarchy  = typename PHARE::amr::SAMRAI_Types::hierarchy_t;

    using HybridModel = typename PHARETypes::HybridModel_t;
    using MHDModel    = typename PHARETypes::MHDModel_t;

    using DiagnosticModelView = PHARE::diagnostic::AMRDiagnosticModelView<Hierarchy, HybridModel>;
    using DiagnosticWriter = PHARE::diagnostic::h5::HighFiveDiagnosticWriter<DiagnosticModelView>;

    std::unique_ptr<DiagnosticModelView> modelView;
    std::unique_ptr<DiagnosticWriter> writer;
    std::unique_ptr<PHARE::diagnostic::DiagnosticsManager<DiagnosticWriter>> dMan;

    auto& dict()
    {
        StaticIntepreter::INSTANCE();
        return PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
    }

    TestSimulator()
        : Simulator{dict()}
    {
        Simulator::initialize();

        if (dict()["simulation"].contains("diagnostics"))
        {
            modelView = std::make_unique<DiagnosticModelView>(*this->getPrivateHierarchy(),
                                                              *this->getHybridModel());
            writer    = DiagnosticWriter::from(*modelView, dict()["simulation"]["diagnostics"]);
            dMan      = PHARE::diagnostic::DiagnosticsManager<DiagnosticWriter>::from(
                *writer, dict()["simulation"]["diagnostics"]);
        }
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

namespace PHARE_test
{
class StreamAppender : public SAMRAI::tbox::Logger::Appender
{
public:
    StreamAppender(std::ostream* stream) { d_stream = stream; }
    void logMessage(const std::string& message, const std::string& filename, const int line)
    {
        (*d_stream) << "At :" << filename << " line :" << line << " message: " << message
                    << std::endl;
    }

private:
    std::ostream* d_stream;
};

class SamraiLifeCycle
{
public:
    SamraiLifeCycle(int argc, char** argv,
                    std::vector<std::function<void()>> funcs = std::vector<std::function<void()>>{})
    {
        SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
        SAMRAI::tbox::SAMRAIManager::initialize();
        SAMRAI::tbox::SAMRAIManager::startup();

        std::shared_ptr<SAMRAI::tbox::Logger::Appender> appender
            = std::make_shared<StreamAppender>(StreamAppender{&std::cout});
        SAMRAI::tbox::Logger::getInstance()->setWarningAppender(appender);
        for (auto& func : funcs)
            func();
    }
    ~SamraiLifeCycle()
    {
        for (auto& func : funcs_)
            func();
        PHARE::initializer::PHAREDictHandler::INSTANCE().stop();
        StaticIntepreter::input.reset();
        SAMRAI::tbox::SAMRAIManager::shutdown();
        SAMRAI::tbox::SAMRAIManager::finalize();
        SAMRAI::tbox::SAMRAI_MPI::finalize();
    }

    SamraiLifeCycle& doOnDestruct(std::vector<std::function<void()>> funcs)
    {
        funcs_ = funcs;
        return *this;
    }

private:
    std::vector<std::function<void()>> funcs_;
};

} // namespace PHARE_test

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
