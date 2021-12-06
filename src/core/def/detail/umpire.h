

#ifndef PHARE_CORE_DEF_DETAIL_UMPIRE_H
#define PHARE_CORE_DEF_DETAIL_UMPIRE_H


#if defined(HAVE_SYS_TIMES_H)
#undef HAVE_SYS_TIMES_H // https://github.com/LLNL/SAMRAI/issues/93
#undef HAVE_UNISTD_H
#endif

// #include "SAMRAI/SAMRAI_config.h"

// #if defined(HAVE_UMPIRE)
// #include "kul/log.hpp" // TORM
// #endif

#if !defined(HAVE_UMPIRE)
#define PHARE_HAVE_UMPIRE 0
#else
#define PHARE_HAVE_UMPIRE 1
#endif


#if PHARE_HAVE_UMPIRE
#define PHARE_WITH_UMPIRE(...) __VA_ARGS__
#else
#define PHARE_WITH_UMPIRE(...)
#endif


#if PHARE_HAVE_UMPIRE
#include "SAMRAI/tbox/Collectives.h" // tbox::parallel_synchronize();
#include "umpire/Allocator.hpp"
#include "umpire/ResourceManager.hpp"
#include "umpire/TypedAllocator.hpp"
#include "umpire/strategy/AllocationAdvisor.hpp"
#include "umpire/strategy/QuickPool.hpp"

namespace PHARE::_umpire_
{
inline auto static_init()
{
    auto& rm = umpire::ResourceManager::getInstance();

    if (rm.isAllocator("samrai::data_allocator")
        and rm.getAllocator("samrai::data_allocator").getPlatform() != umpire::Platform::host)
        throw std::runtime_error("Invalid samrai allocator state");

    if (!rm.isAllocator("samrai::data_allocator"))
    { // initialize samrai internals allocator to HOST
        [[maybe_unused]] auto samrai_allocator = rm.makeAllocator<umpire::strategy::QuickPool>(
            "samrai::data_allocator", rm.getAllocator(umpire::resource::Host));
        assert(samrai_allocator.getPlatform() == umpire::Platform::host);
    }
    assert(rm.isAllocator("samrai::data_allocator"));
    assert(rm.getAllocator("samrai::data_allocator").getPlatform() == umpire::Platform::host);


    if (!rm.isAllocator("samrai::tag_allocator"))
    { // initialize samrai internals allocator to HOST
        [[maybe_unused]] auto samrai_allocator = rm.makeAllocator<umpire::strategy::QuickPool>(
            "samrai::tag_allocator", rm.getAllocator(umpire::resource::Unified));
        assert(samrai_allocator.getPlatform() != umpire::Platform::host);
    }
    assert(rm.isAllocator("samrai::tag_allocator"));
    assert(rm.getAllocator("samrai::tag_allocator").getPlatform() != umpire::Platform::host);

    if (!rm.isAllocator("PHARE::data_allocator"))
    {
        auto allocator = rm.makeAllocator<umpire::strategy::AllocationAdvisor>(
            "internal::PHARE::um_allocation_advisor", rm.getAllocator(umpire::resource::Unified),
            // Set preferred location to GPU
            "PREFERRED_LOCATION");

        [[maybe_unused]] auto phare_allocator
            = rm.makeAllocator<umpire::strategy::QuickPool>("PHARE::data_allocator", allocator);
        assert(phare_allocator.getPlatform() != umpire::Platform::host);
        assert(rm.isAllocator("PHARE::data_allocator"));
    }

    return true;
}
inline static auto _static_sim_init = static_init();
} // namespace PHARE::_umpire_




#endif /* PHARE_HAVE_UMPIRE */

#endif /* PHARE_CORE_DEF_DETAIL_UMPIRE_H */
