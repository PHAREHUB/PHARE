#ifndef PHARE_TEST_AMR_AMR_HPP
#define PHARE_TEST_AMR_AMR_HPP

#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

namespace PHARE::test::amr
{
class SamraiLifeCycle
{
public:
    SamraiLifeCycle(int argc = 0, char** argv = nullptr)
    {
        SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
        SAMRAI::tbox::SAMRAIManager::initialize();
        SAMRAI::tbox::SAMRAIManager::startup();
    }
    ~SamraiLifeCycle()
    {
        SAMRAI::tbox::SAMRAIManager::shutdown();
        SAMRAI::tbox::SAMRAIManager::finalize();
        SAMRAI::tbox::SAMRAI_MPI::finalize();
    }

    static void reset()
    {
        SAMRAI::tbox::SAMRAIManager::shutdown();
        SAMRAI::tbox::SAMRAIManager::startup();
    }
};

} // namespace PHARE::test::amr

#endif /* PHARE_TEST_AMR_AMR_HPP */
