#ifndef PHARE_PYTHON_DATA_WRANGLER_HPP
#define PHARE_PYTHON_DATA_WRANGLER_HPP


#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <vector>

#include "simulator/simulator.hpp"
#include "amr/wrappers/hierarchy.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "core/utilities/mpi_utils.hpp"
#include "core/utilities/point/point.hpp"
#include "cppdict/include/dict.hpp"
#include "initializer/data_provider.hpp"
#include "python3/patch_data.hpp"
#include "python3/patch_level.hpp"

namespace PHARE::pydata
{
template<std::size_t dimension, std::size_t interp_order, std::size_t nbRefinedPart,
         template<typename> typename MHDTimeStepper>
class SimulatorCaster
{
public:
    using Simulator_t = Simulator<dimension, interp_order, nbRefinedPart, MHDTimeStepper>;

    SimulatorCaster(std::shared_ptr<ISimulator> const& _simulator)
        : simulator{_simulator}
    {
    }

    template<typename Dimension, typename InterpOrder, typename NbRefinedPart>
    Simulator_t* operator()(std::size_t userDim, std::size_t userInterpOrder,
                            std::size_t userNbRefinedPart, Dimension dimension_fn,
                            InterpOrder interp_order_fn, NbRefinedPart nbRefinedPart_fn)
    {
        if (userDim == dimension_fn() and userInterpOrder == interp_order_fn()
            and userNbRefinedPart == nbRefinedPart_fn())
        {
            std::size_t constexpr d  = dimension_fn();
            std::size_t constexpr io = interp_order_fn();
            std::size_t constexpr nb = nbRefinedPart_fn();

            // extra if constexpr as cast is templated and not generic interface
            if constexpr (d == dimension and io == interp_order and nb == nbRefinedPart)
                return dynamic_cast<Simulator_t*>(simulator.get());
        }
        return nullptr;
    }

private:
    std::shared_ptr<ISimulator> const& simulator;
};



template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart,
         template<typename> typename MHDTimeStepper>
class DataWrangler
{
public:
    static constexpr std::size_t dimension     = _dimension;
    static constexpr std::size_t interp_order  = _interp_order;
    static constexpr std::size_t nbRefinedPart = _nbRefinedPart;

    using Simulator   = PHARE::Simulator<dimension, interp_order, nbRefinedPart, MHDTimeStepper>;
    using HybridModel = typename Simulator::HybridModel;

    DataWrangler(std::shared_ptr<ISimulator> const& simulator,
                 std::shared_ptr<amr::Hierarchy> const& hierarchy)
        : simulator_{cast_simulator(simulator)}
        , hierarchy_{hierarchy}

    {
    }


    auto getNumberOfLevels() const { return hierarchy_->getNumberOfLevels(); }

    auto getPatchLevel(size_t lvl)
    {
        return PatchLevel<_dimension, _interp_order, _nbRefinedPart, MHDTimeStepper>{
            *hierarchy_, *simulator_.getHybridModel(), lvl};
    }

    auto sort_merge_1d(std::vector<PatchData<std::vector<double>, dimension>> const&& input,
                       bool shared_patch_border = false)
    {
        std::vector<std::pair<double, PatchData<std::vector<double>, dimension> const*>> sorted;
        for (auto const& data : input)
            sorted.emplace_back(core::Point<double, 1>::fromString(data.origin)[0], &data);
        std::sort(sorted.begin(), sorted.end(), [](auto& a, auto& b) { return a.first < b.first; });
        std::vector<double> ret;
        for (size_t i = 0; i < sorted.size(); i++)
        { // skip empty patches in case of unequal patches across MPI domains
            if (!sorted[i].second->data.size())
                continue;
            auto& data   = sorted[i].second->data;
            auto& ghosts = sorted[i].second->nGhosts;
            auto end     = ghosts;
            // primal nodes share a cell wall when patches touch so drop duplicate value if so
            if (shared_patch_border)
                end = i == sorted.size() - 1 ? end : end + 1;
            ret.insert(std::end(ret), std::begin(data) + ghosts, std::end(data) - end);
        }
        return ret;
    }

    auto sync(std::vector<PatchData<std::vector<double>, dimension>> const& input)
    {
        int mpi_size = core::mpi::size();
        std::vector<PatchData<std::vector<double>, dimension>> collected;

        auto reinterpret_array = [&](auto& py_array) {
            return reinterpret_cast<std::array<std::size_t, dimension>&>(
                *static_cast<std::size_t*>(py_array.request().ptr));
        };

        auto collect = [&](PatchData<std::vector<double>, dimension> const& patch_data) {
            auto patchIDs = core::mpi::collect(patch_data.patchID, mpi_size);
            auto origins  = core::mpi::collect(patch_data.origin, mpi_size);
            auto lower    = core::mpi::collect_raw(makeSpan(patch_data.lower), mpi_size);
            auto upper    = core::mpi::collect_raw(makeSpan(patch_data.upper), mpi_size);
            auto ghosts   = core::mpi::collect(patch_data.nGhosts, mpi_size);
            auto datas    = core::mpi::collect(patch_data.data, mpi_size);

            for (int i = 0; i < mpi_size; i++)
            {
                auto& data = collected.emplace_back();
                setPatchData(data, patchIDs[i], origins[i], lower[i], upper[i]);
                data.nGhosts = ghosts[i];
                data.data    = std::move(datas[i]);
            }
        };

        std::size_t max = core::mpi::max(input.size(), mpi_size);

        PatchData<std::vector<double>, dimension> empty;

        for (size_t i = 0; i < max; i++)
        {
            if (i < input.size())
                collect(input[i]);
            else
                collect(empty);
        }
        return collected;
    }

    auto sync_merge(std::vector<PatchData<std::vector<double>, dimension>> const& input,
                    [[maybe_unused]] bool primal)
    {
        if constexpr (dimension == 1)
            return sort_merge_1d(sync(input), primal);

        throw std::runtime_error("Not handled for >1 dim");
    }

private:
    Simulator& simulator_;
    std::shared_ptr<amr::Hierarchy> hierarchy_;




    static Simulator& cast_simulator(std::shared_ptr<ISimulator> const& simulator)
    {
        using SimulatorCaster
            = SimulatorCaster<dimension, interp_order, nbRefinedPart, MHDTimeStepper>;

        auto const& simDict = initializer::PHAREDictHandler::INSTANCE().dict()["simulation"];

        Simulator* simulator_ptr = core::makeAtRuntime<SimulatorCaster>(
            simDict["dimension"].template to<int>(), simDict["interp_order"].template to<int>(),
            simDict["refined_particle_nbr"].template to<int>(), SimulatorCaster{simulator});
        if (!simulator_ptr)
            throw std::runtime_error("Data Wranger creation error: failed to cast Simulator");

        return *simulator_ptr;
    }
};
} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_DATA_WRANGLER_H*/
