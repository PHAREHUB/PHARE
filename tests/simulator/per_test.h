#ifndef PHARE_TEST_SIMULATOR_PER_TEST_H
#define PHARE_TEST_SIMULATOR_PER_TEST_H

#include "simulator/simulator.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

// Having multiple PythonDataProvider per binary execution doesn't work so well

struct __attribute__((visibility("hidden"))) StaticIntepreter
{
    std::shared_ptr<PHARE::initializer::PythonDataProvider> input;

    StaticIntepreter()
        : input{std::make_shared<PHARE::initializer::PythonDataProvider>()}
    {
        input->read();
    }

    void kill()
    {
        PHARE::initializer::PHAREDictHandler::INSTANCE().stop();
        input.reset();
    }

    static StaticIntepreter& INSTANCE()
    {
        static StaticIntepreter i;
        return i;
    }
};

template<size_t dim, size_t interp>
struct TestSimulator : public PHARE::Simulator<dim, interp>
{
    using Simulator                 = PHARE::Simulator<dim, interp>;
    static constexpr size_t dim_    = dim;
    static constexpr size_t interp_ = interp;

    auto& dict()
    {
        StaticIntepreter::INSTANCE();
        return PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
    }

    TestSimulator()
        : Simulator{dict()}
    {
        Simulator::initialize();
    }
};

template<typename Simulator>
struct SimulatorTest : public ::testing::Test
{
};

using Simulators = testing::Types<TestSimulator<1, 1>/*, TestSimulator<1, 2>, TestSimulator<1, 3>*//*,   // dim 1
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
        StaticIntepreter::INSTANCE().kill();
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

#endif /* PHARE_TEST_SIMULATOR_PER_TEST_H */