
#include "core/def/phlop.hpp" // IWYU pragma: keep // scope timing
#include "core/utilities/mpi_utils.hpp"

#include "initializer/data_provider.hpp"

#include "samrai.hpp"


#include <SAMRAI/tbox/SAMRAIManager.h>


#include "H5pubconf.h" // may define H5_HAVE_SUBFILING_VFD

#if !defined(H5_HAVE_SUBFILING_VFD)
#define H5_HAVE_SUBFILING_VFD 0
#endif


namespace PHARE
{

SamraiLifeCycle::SamraiLifeCycle(int argc, char** argv)
{
#if H5_HAVE_SUBFILING_VFD
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided < MPI_THREAD_MULTIPLE)
        throw std::runtime_error(
            "MPI_THREAD_MULTIPLE required for HDF5 subfiling but not provided");
    SAMRAI::tbox::SAMRAI_MPI::init(MPI_COMM_WORLD);

#else  // normal way
    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
#endif // H5_HAVE_SUBFILING_VFD

    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();
    // uncomment next line for debugging samrai issues
    // SAMRAI::tbox::SAMRAI_MPI::setCallAbortInParallelInsteadOfMPIAbort();
    std::shared_ptr<SAMRAI::tbox::Logger::Appender> appender
        = std::make_shared<StreamAppender>(StreamAppender{&std::cout});
    SAMRAI::tbox::Logger::getInstance()->setWarningAppender(appender);
    PHARE_WITH_PHLOP( //
        if (auto e = core::get_env("PHARE_SCOPE_TIMING", "false"); e == "1" || e == "true")
            phlop::ScopeTimerMan::INSTANCE()
                .file_name(".phare/timings/rank." + std::to_string(core::mpi::rank()) + ".txt")
                .init(); //
    )
}

SamraiLifeCycle::~SamraiLifeCycle()
{
    PHARE_WITH_PHLOP(phlop::ScopeTimerMan::reset());
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();

    SAMRAI::tbox::SAMRAI_MPI::finalize();
#if H5_HAVE_SUBFILING_VFD
    MPI_Finalize();
#endif // H5_HAVE_SUBFILING_VFD
}

void SamraiLifeCycle::reset()
{
    PHARE::initializer::PHAREDictHandler::INSTANCE().stop();
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::startup();
    getRestartManager()->clearRestartItems();
}


SAMRAI::hier::VariableDatabase* SamraiLifeCycle::getDatabase()
{
    return SAMRAI::hier::VariableDatabase::getDatabase();
}

SAMRAI::hier::PatchDataRestartManager* SamraiLifeCycle::getPatchDataRestartManager()
{
    return SAMRAI::hier::PatchDataRestartManager::getManager();
}


SAMRAI::tbox::RestartManager* SamraiLifeCycle::getRestartManager()
{
    return SAMRAI::tbox::RestartManager::getManager();
}


} // namespace PHARE
