
#ifndef PHARE_AMR_SAMRAI_HPP
#define PHARE_AMR_SAMRAI_HPP

#include "core/def/phlop.hpp" // scope timing

#include "amr/wrappers/integrator.hpp"
#include "core/utilities/mpi_utils.hpp"

#include <SAMRAI/hier/VariableDatabase.h>
#include <SAMRAI/tbox/RestartManager.h>
#include <SAMRAI/hier/PatchDataRestartManager.h>

#include <memory>
#include <iostream>

namespace PHARE
{
class StreamAppender : public SAMRAI::tbox::Logger::Appender
{
public:
    StreamAppender(std::ostream* stream) { d_stream = stream; }
    void logMessage(std::string const& message, std::string const& filename, int const line)
    {
        (*d_stream) << "At :" << filename << " line :" << line << " message: " << message
                    << std::endl;
    }

private:
    std::ostream* d_stream;
};

class SamraiLifeCycle //
{
public:
    SamraiLifeCycle(int argc = 0, char** argv = nullptr);

    ~SamraiLifeCycle();

    static void reset();

    static SAMRAI::hier::VariableDatabase* getDatabase();

    static SAMRAI::hier::PatchDataRestartManager* getPatchDataRestartManager();

    static SAMRAI::tbox::RestartManager* getRestartManager();
};


} // namespace PHARE


#endif /*PHARE_AMR_SAMRAI_HPP*/
