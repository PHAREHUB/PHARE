
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
 * @brief Custom logger appender for SAMRAI to redirect log messages to a stream
 * 
 * This class allows SAMRAI warning and error messages to be directed to 
 * a specific output stream (e.g., std::cout, std::cerr, or a file stream).
 */
class StreamAppender : public SAMRAI::tbox::Logger::Appender
{
public:
    /**
     * @brief Construct appender with target output stream
     * @param stream Pointer to the output stream for log messages
     */
    StreamAppender(std::ostream* stream) { d_stream = stream; }
    
    /**
     * @brief Format and output a log message
     * @param message The log message content
     * @param filename Source file where the message originated
     * @param line Line number where the message originated
     */
    void logMessage(std::string const& message, std::string const& filename, int const line)
    {
        (*d_stream) << "At :" << filename << " line :" << line << " message: " << message
                    << std::endl;
    }

private:
    std::ostream* d_stream;
};

/**
 * @brief Manages SAMRAI library lifecycle and provides centralized access to SAMRAI singletons
 * 
 * This class encapsulates SAMRAI initialization, shutdown, and access to global SAMRAI managers.
 * It ensures proper initialization order and cleanup of SAMRAI resources.
 * 
 * Key responsibilities:
 * - Initialize/finalize SAMRAI MPI and manager subsystems
 * - Configure SAMRAI logging appenders
 * - Initialize optional performance timing (PHLOP)
 * - Provide static access to SAMRAI singleton managers (avoiding multiple static instances)
 * 
 * @note This class exists to solve the static singleton problem when SAMRAI is built as
 *       shared libraries and loaded by multiple Python modules. Each module gets the same
 *       SAMRAI instance through these static accessors.
 */
class SamraiLifeCycle
{
public:
    /**
     * @brief Initialize SAMRAI subsystems
     * @param argc Command line argument count (default: 0)
     * @param argv Command line arguments (default: nullptr)
     */
    SamraiLifeCycle(int argc = 0, char** argv = nullptr);

    /**
     * @brief Cleanup SAMRAI subsystems in proper order
     */
    ~SamraiLifeCycle();

    /**
     * @brief Reset SAMRAI to initial state (for simulation reruns)
     * 
     * Shuts down and restarts SAMRAI manager, clears restart items.
     * Used when running multiple simulations in the same process.
     */
    static void reset();

    /**
     * @brief Get SAMRAI's global variable database
     * @return Pointer to the singleton VariableDatabase
     */
    static SAMRAI::hier::VariableDatabase* getDatabase();

    /**
     * @brief Get SAMRAI's patch data restart manager
     * @return Pointer to the singleton PatchDataRestartManager
     */
    static SAMRAI::hier::PatchDataRestartManager* getPatchDataRestartManager();

    /**
     * @brief Get SAMRAI's restart manager
     * @return Pointer to the singleton RestartManager
     */
    static SAMRAI::tbox::RestartManager* getRestartManager();
};


} // namespace PHARE


#endif /*PHARE_AMR_SAMRAI_HPP*/
