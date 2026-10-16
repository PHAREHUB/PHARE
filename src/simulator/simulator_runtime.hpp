#ifndef PHARE_SIMULATOR_SIMULATOR_RUNTIME_HPP
#define PHARE_SIMULATOR_SIMULATOR_RUNTIME_HPP

#include "core/logger.hpp"

#include "initializer/data_provider.hpp"

#include "amr/wrappers/hierarchy.hpp"

#include "simulator/simulator_def.hpp"

namespace PHARE
{

template<typename Maker>
std::unique_ptr<ISimulator> makeSimulatorAtRuntime(SimOpts const& r_opts, Maker&& maker)
{
    std::unique_ptr<ISimulator> p = nullptr;

    auto constexpr sims = phare_exe_default_simulators();
    core::for_N<sims.size()>([&](auto ic) {
        constexpr auto c_opts = sims[ic()];
        if (!p)
            p = maker.template operator()<c_opts>(r_opts);
    });

    return p;
}

template<typename Maker> // used from PHARE::amr::Hierarchy
auto makeHierarchyAtRuntime(std::size_t dim, Maker&& maker)
{
    using Ptr_t = decltype(maker(dim, 1));
    Ptr_t p{};

    auto constexpr sims = phare_exe_default_simulators();
    core::for_N<sims.size()>([&](auto ic) {
        constexpr auto opts = sims[ic()];

        if (!p)
            p = maker(dim, core::DimConst<opts.dimension>{});
    });

    return p;
}

inline auto make_hierarchy()
{
    PHARE::initializer::PHAREDict const& theDict
        = PHARE::initializer::PHAREDictHandler::INSTANCE().dict();
    auto dim  = theDict["simulation"]["dimension"].template to<int>();
    auto hier = makeHierarchyAtRuntime<amr::HierarchyMaker>(dim, amr::HierarchyMaker{theDict});
    if (hier)
        return hier;
    PHARE_LOG_LINE_SS("hierarchy not found for params:\n"
                      << dim << " " << theDict["simulation"]["interp_order"].template to<int>()
                      << " " << theDict["simulation"]["refined_particle_nbr"].template to<int>());
    throw std::runtime_error("Likely unsupported template parameters");
}

} // namespace PHARE

#endif /*PHARE_SIMULATOR_SIMULATOR_DEF_HPP*/
