#ifndef PHARE_PHARE_HPPPP
#define PHARE_PHARE_HPPPP

#include "amr/wrappers/hierarchy.hpp"
#include "initializer/data_provider.hpp"
#include "simulator/simulator.hpp"

namespace PHARE
{

struct SimulatorMaker
{
    SimulatorMaker(std::shared_ptr<PHARE::amr::Hierarchy>& hierarchy)
        : hierarchy_{hierarchy}
    {
    }

    std::shared_ptr<PHARE::amr::Hierarchy>& hierarchy_;

    template<SimOpts opts>
    std::unique_ptr<ISimulator> operator()(SimOpts const& userOpts)
    {
        if (userOpts == opts)
            return std::make_unique<Simulator<opts>>(
                initializer::PHAREDictHandler::INSTANCE().dict(), hierarchy_);

        return nullptr;
    }
};

} // namespace PHARE

#endif /*PHARE_PHARE_HPP*/
