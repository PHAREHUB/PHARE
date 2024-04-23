
#ifndef PHARE_PHARE_INCLUDE_HPP
#define PHARE_PHARE_INCLUDE_HPP

#include "simulator/simulator.hpp"
#include "core/utilities/algorithm.hpp"
#include "core/utilities/run_timer.hpp"


#include <memory>
#include <iostream>

namespace PHARE
{
class StreamAppender : public SAMRAI::tbox::Logger::Appender
{
public:
    StreamAppender(std::ostream* stream) { d_stream = stream; }
    void logMessage(std::string const& message, std::string const& filename, const int line)
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
    SamraiLifeCycle(int argc = 0, char** argv = nullptr)
    {
        SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
        SAMRAI::tbox::SAMRAIManager::initialize();
        SAMRAI::tbox::SAMRAIManager::startup();
        // SAMRAI::tbox::SAMRAI_MPI::setCallAbortInParallelInsteadOfMPIAbort();
        std::shared_ptr<SAMRAI::tbox::Logger::Appender> appender
            = std::make_shared<StreamAppender>(StreamAppender{&std::cout});
        SAMRAI::tbox::Logger::getInstance()->setWarningAppender(appender);
    }
    ~SamraiLifeCycle()
    {
        core::StaticRunTimer::reset();
        SAMRAI::tbox::SAMRAIManager::shutdown();
        SAMRAI::tbox::SAMRAIManager::finalize();
        SAMRAI::tbox::SAMRAI_MPI::finalize();
    }

    static void reset()
    {
        core::StaticRunTimer::reset();
        PHARE::initializer::PHAREDictHandler::INSTANCE().stop();
        SAMRAI::tbox::SAMRAIManager::shutdown();
        SAMRAI::tbox::SAMRAIManager::startup();
    }
};


} // namespace PHARE


#endif /*PHARE_PHARE_INCLUDE_H*/
