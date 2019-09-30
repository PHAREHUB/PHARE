#ifndef PHARE_TEST_DIAGNOSTIC_SAMRAI_LIFECYCLE_H
#define PHARE_TEST_DIAGNOSTIC_SAMRAI_LIFECYCLE_H

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
    SamraiLifeCycle(int argc, char** argv)
    {
        SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
        SAMRAI::tbox::SAMRAIManager::initialize();
        SAMRAI::tbox::SAMRAIManager::startup();

        std::shared_ptr<SAMRAI::tbox::Logger::Appender> appender
            = std::make_shared<StreamAppender>(StreamAppender{&std::cout});
        SAMRAI::tbox::Logger::getInstance()->setWarningAppender(appender);
    }
    ~SamraiLifeCycle()
    {
        SAMRAI::tbox::SAMRAIManager::shutdown();
        SAMRAI::tbox::SAMRAIManager::finalize();
        SAMRAI::tbox::SAMRAI_MPI::finalize();
    }
};
} // namespace PHARE_test

#endif /*PHARE_TEST_DIAGNOSTIC_SAMRAI_LIFECYCLE_H*/