#ifndef PHARE_PYTHON_DATA_WRANGLER_HPP
#define PHARE_PYTHON_DATA_WRANGLER_HPP


#include "core/utilities/mpi_utils.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/meta/meta_utilities.hpp"

#include "amr/wrappers/hierarchy.hpp"

#include "core/utilities/types.hpp"
#include "initializer/data_provider.hpp"

#include "python3/patch_data.hpp"
#include "python3/patch_level.hpp"

#include "python3/pybind_def.hpp"
#include "simulator/simulator.hpp"

#include "dict.hpp"


#include <array>
#include <memory>
#include <vector>
#include <cstddef>
#include <iterator>
#include <algorithm>
#include <stdexcept>



namespace PHARE::pydata
{
template<auto opts>
class __attribute__((visibility("hidden"))) SimulatorCaster
{
public:
    using Simulator_t = Simulator<opts>;

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
            if constexpr (d == opts.dimension and io == opts.interp_order
                          and nb == opts.nbRefinedPart)
                return dynamic_cast<Simulator_t*>(simulator.get());
        }
        return nullptr;
    }

private:
    std::shared_ptr<ISimulator> const& simulator;
};



template<auto opts>
class __attribute__((visibility("hidden"))) DataWrangler
{
public:
    static constexpr std::size_t dimension = opts.dimension;

    using Simulator   = PHARE::Simulator<opts>;
    using HybridModel = Simulator::HybridModel;
    using MHDModel    = Simulator::MHDModel;



    DataWrangler(std::shared_ptr<Simulator> const& simulator,
                 std::shared_ptr<amr::Hierarchy> const& hierarchy)
        : simulator_{*simulator}
        , hierarchy_{hierarchy}
    {
    }

    DataWrangler(std::shared_ptr<ISimulator> const& simulator,
                 std::shared_ptr<amr::Hierarchy> const& hierarchy)
        : simulator_{cast_simulator(simulator)}
        , hierarchy_{hierarchy}
    {
    }


    auto getNumberOfLevels() const { return hierarchy_->getNumberOfLevels(); }

    auto getMHDPatchLevel(size_t lvl)
    {
        return PatchLevel<MHDModel>{*hierarchy_, *simulator_.getMHDModel(), lvl};
    }

    auto getHybridPatchLevel(size_t lvl)
    {
        return PatchLevel<HybridModel>{*hierarchy_, *simulator_.getHybridModel(), lvl};
    }


    auto sync(std::vector<PatchData<py_array_t<double>, dimension>> const& input)
    {
        // collect all data on rank 0!

        auto const mpi_size = core::mpi::size();
        std::vector<PatchData<py_array_t<double>, dimension>> collected;

        auto const collect = [&](PatchData<py_array_t<double>, dimension> const& patch_data) {
            auto patchIDs = core::mpi::collect(patch_data.patchID, mpi_size);
            auto shapes   = core::mpi::collectArrays(shape<dimension>(patch_data.data), mpi_size);
            auto origins  = core::mpi::collect(makeSpan(patch_data.origin), mpi_size);
            auto lower    = core::mpi::collect(makeSpan(patch_data.lower), mpi_size);
            auto upper    = core::mpi::collect(makeSpan(patch_data.upper), mpi_size);
            auto ghosts   = core::mpi::collect(patch_data.nGhosts, mpi_size);
            auto datas    = core::mpi::collect(makeSpan(patch_data.data), mpi_size);

            if (core::mpi::rank() == 0)
                for (int i = 0; i < mpi_size; ++i)
                {
                    if (datas[i].size() == 0) // missing patch on rank
                        continue;

                    auto& data = collected.emplace_back(shapes[i], strides_from<double>(shapes[i]));
                    auto const span = makeSpan(data.data);
                    std::memcpy(span.data(), datas[i].data(), span.size() * sizeof(double));
                    setPatchData(data, patchIDs[i], origins[i], lower[i], upper[i]);
                    data.nGhosts = ghosts[i];
                }
        };

        auto const max = core::mpi::max(input.size());

        PatchData<py_array_t<double>, dimension> empty{core::ConstArray<int, dimension>()};

        for (std::size_t i = 0; i < max; ++i)
            collect(i < input.size() ? input[i] : empty);

        return collected;
    }


private:
    Simulator& simulator_;
    std::shared_ptr<amr::Hierarchy> hierarchy_;


    static Simulator& cast_simulator(std::shared_ptr<ISimulator> const& simulator)
    {
        using SimulatorCaster = SimulatorCaster<opts>;

        auto const& simDict = initializer::PHAREDictHandler::INSTANCE().dict()["simulation"];

        Simulator* simulator_ptr = core::makeAtRuntime<SimulatorCaster>(
            simDict["dimension"].template to<int>(), simDict["interp_order"].template to<int>(),
            simDict["refined_particle_nbr"].template to<int>(), SimulatorCaster{simulator});
        if (!simulator_ptr)
            throw std::runtime_error("Data Wrangler creation error: failed to cast Simulator");

        return *simulator_ptr;
    }
};
} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_DATA_WRANGLER_H*/
