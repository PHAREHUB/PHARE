#ifndef PHARE_PYTHON_DATA_WRANGLER_H
#define PHARE_PYTHON_DATA_WRANGLER_H


#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <vector>
#include "amr/wrappers/hierarchy.h"
#include "core/utilities/meta/meta_utilities.h"
#include "core/utilities/mpi_utils.h"
#include "core/utilities/point/point.h"
#include "cppdict/include/dict.hpp"
#include "initializer/data_provider.h"
#include "python3/patch_data.h"
#include "python3/patch_level.h"
#include "simulator/simulator.h"

namespace PHARE::pydata
{
template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart>
class DataWrangler
{
public:
    static constexpr std::size_t dimension     = _dimension;
    static constexpr std::size_t interp_order  = _interp_order;
    static constexpr std::size_t nbRefinedPart = _nbRefinedPart;

    using Simulator   = PHARE::Simulator<dimension, interp_order, nbRefinedPart>;
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
        return PatchLevel<_dimension, _interp_order, _nbRefinedPart>{
            *hierarchy_, *simulator_.getHybridModel(), lvl};
    }

    auto sort_merge_1d(std::vector<PatchData<std::vector<double>, dimension>> const&& input,
                       bool shared_patch_border = false)
    {
        std::vector<std::pair<double, const PatchData<std::vector<double>, dimension>*>> sorted;
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

        auto collect = [&](auto& patch_data) {
            auto patchIDs = core::mpi::collect(patch_data.patchID, mpi_size);
            auto origins  = core::mpi::collect(patch_data.origin, mpi_size);
            auto lower    = core::mpi::collect(reinterpret_array(patch_data.lower), mpi_size);
            auto upper    = core::mpi::collect(reinterpret_array(patch_data.upper), mpi_size);
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

        auto max = core::mpi::max(input.size(), mpi_size);

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
                    bool primal)
    {
        if constexpr (dimension == 1)
            return sort_merge_1d(sync(input), primal);

        static_assert("Unhandled dimension: sort_merge_fluid");
    }

private:
    Simulator& simulator_;
    std::shared_ptr<amr::Hierarchy> hierarchy_;


    static Simulator& cast_simulator(std::shared_ptr<ISimulator> const& simulator)
    {
        using SimulatorCaster = PHARE::SimulatorCaster<dimension, interp_order, nbRefinedPart>;

        auto simDict = initializer::PHAREDictHandler::INSTANCE().dict()["simulation"];

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
