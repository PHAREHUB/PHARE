
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
/**
 * @brief Custom logger appender for SAMRAI that redirects log messages to a stream.
 * 
 * This class allows SAMRAI library logging to be captured and redirected to
 * a custom output stream (e.g., std::cout, std::cerr, or a file).
 */
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

/**
 * @brief Manages the lifecycle of SAMRAI library initialization and shutdown.
 * 
 * This RAII class ensures SAMRAI is properly initialized before use and
 * cleanly shut down when done. It also provides access to SAMRAI's singleton
 * managers for variables, restart data, and restart operations.
 * 
 * Note: Only one instance should exist per process. The instance can be reset
 * via the reset() static method if needed (e.g., when switching between
 * different simulator configurations).
 */
class SamraiLifeCycle //
{
public:
    SamraiLifeCycle(int argc = 0, char** argv = nullptr);

    ~SamraiLifeCycle();

    /**
     * @brief Reset SAMRAI library state.
     * 
     * This should be called when switching between different simulator instances
     * or when reinitializing the library is required.
     */
    static void reset();

    /**
     * @brief Get SAMRAI's variable database singleton.
     * 
     * The variable database manages all SAMRAI variables (field data, particle data, etc.)
     * and their associated patch data indices.
     */
    static SAMRAI::hier::VariableDatabase* getDatabase();

    /**
     * @brief Get SAMRAI's patch data restart manager singleton.
     * 
     * This manager handles registration and unregistration of patch data that
     * should be included in restart files.
     */
    static SAMRAI::hier::PatchDataRestartManager* getPatchDataRestartManager();

    /**
     * @brief Get SAMRAI's restart manager singleton.
     * 
     * This manager handles reading and writing restart files for checkpointing
     * and simulation restarts.
     */
    static SAMRAI::tbox::RestartManager* getRestartManager();
};


} // namespace PHARE


#endif /*PHARE_AMR_SAMRAI_HPP*/
