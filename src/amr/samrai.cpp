
#include "core/def/pragma_disable.hpp"
#include "core/def/phlop.hpp" // IWYU pragma: keep // scope timing
#include "mpi/mpi_utils.hpp"

#include "initializer/data_provider.hpp"

#include "samrai.hpp"

#include <SAMRAI/tbox/SAMRAIManager.h>

namespace PHARE
{

// GCC (observed on 16.x, -O3) mis-attributes the template argument of
// _Sp_counted_ptr_inplace<Tp,...>::_M_dispose() for this shared_ptr<Appender>, reporting it as
// instantiated for std::string instead of StreamAppender -- a false positive.
DISABLE_WARNING(array-bounds, array-bounds, 42)
SamraiLifeCycle::SamraiLifeCycle(int argc, char** argv)
{
    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();
    // uncomment next line for debugging samrai issues
    // SAMRAI::tbox::SAMRAI_MPI::setCallAbortInParallelInsteadOfMPIAbort();
    std::shared_ptr<SAMRAI::tbox::Logger::Appender> appender
        = std::make_shared<StreamAppender>(&std::cout);
    SAMRAI::tbox::Logger::getInstance()->setWarningAppender(appender);
    PHARE_WITH_PHLOP({
        if (auto e = core::get_env("PHARE_SCOPE_TIMING", "false"); e == "1" || e == "true")
            phlop::scope_timer()
                .file_name(".phare/timings/rank." + std::to_string(mpi::rank()) + ".txt")
                .init();
    })
}
ENABLE_WARNING(array-bounds, array-bounds, 42)

SamraiLifeCycle::~SamraiLifeCycle()
{
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();
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
