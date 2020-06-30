#ifndef PHARE_PYTHON_DATA_WRANGLER_H
#define PHARE_PYTHON_DATA_WRANGLER_H

#include "python3/patch_data.h"
#include "python3/patch_level.h"

namespace PHARE::pydata
{
template<std::size_t _dimension, std::size_t _interp_order, std::size_t _nbRefinedPart>
class DataWrangler
{
public:
    using This                                 = DataWrangler;
    static constexpr std::size_t dimension     = _dimension;
    static constexpr std::size_t interp_order  = _interp_order;
    static constexpr std::size_t nbRefinedPart = _nbRefinedPart;

    using PHARETypes  = PHARE_Types<dimension, interp_order, nbRefinedPart>;
    using HybridModel = typename PHARETypes::HybridModel_t;

    DataWrangler(std::shared_ptr<ISimulator> const& simulator,
                 std::shared_ptr<amr::Hierarchy> const& hierarchy)
        : simulator_{simulator}
        , hierarchy_{hierarchy}
    {
        auto simDict = initializer::PHAREDictHandler::INSTANCE().dict()["simulation"];
        if (!core::makeAtRuntime<Maker>(
                simDict["dimension"].template to<int>(), simDict["interp_order"].template to<int>(),
                simDict["refined_particle_nbr"].template to<int>(), Maker{*this}))
            throw std::runtime_error("Runtime diagnostic deduction failed");
    }


    auto getNumberOfLevels() const { return hierarchy_->getNumberOfLevels(); }

    auto getPatchLevel(size_t lvl)
    {
        return PatchLevel<_dimension, _interp_order, _nbRefinedPart>{
            *hierarchy_, *simulator_ptr_->getHybridModel(), lvl};
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
        int mpi_size = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
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
    std::shared_ptr<ISimulator> simulator_;
    std::shared_ptr<amr::Hierarchy> hierarchy_;
    Simulator<dimension, interp_order, nbRefinedPart>* simulator_ptr_ = nullptr;

    struct Maker
    {
        Maker(DataWrangler& _dw)
            : dw{_dw}
        {
        }

        template<typename Dimension, typename InterpOrder, typename NbRefinedPart>
        bool operator()(std::size_t userDim, std::size_t userInterpOrder,
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
                    return (dw.simulator_ptr_
                            = dynamic_cast<PHARE::Simulator<d, io, nb>*>(dw.simulator_.get()));
            }
            return 0;
        }

        DataWrangler& dw;
    };
};
} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_DATA_WRANGLER_H*/
