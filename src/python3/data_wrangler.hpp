#ifndef PHARE_PYTHON_DATA_WRANGLER_HPP
#define PHARE_PYTHON_DATA_WRANGLER_HPP


#include "mpi/mpi_utils.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/meta/meta_utilities.hpp"

#include "amr/wrappers/hierarchy.hpp"

#include "python3/patch_data.hpp"
#include "python3/patch_level.hpp"

#include "python3/pybind_def.hpp"

#include "simulator/simulator.hpp"

#include <array>
#include <memory>
#include <vector>
#include <cstddef>


namespace PHARE::pydata
{

template<auto opts>
class __attribute__((visibility("hidden"))) DataWrangler
{
public:
    static constexpr std::size_t dimension = opts.dimension;

    using Simulator_t = PHARE::Simulator<opts>;

    DataWrangler(std::shared_ptr<Simulator_t> const& simulator,
                 std::shared_ptr<amr::Hierarchy> const& hierarchy)
        : simulator_{*simulator}
        , hierarchy_{hierarchy}
    {
    }


    auto getNumberOfLevels() const { return hierarchy_->getNumberOfLevels(); }

    auto getMHDPatchLevel(size_t lvl)
        requires(has_mhd_v<opts>)
    {
        using MHDModel = solver::PHARE_Types<opts>::MHD::Model_t;
        return PatchLevel<MHDModel>{*hierarchy_, *simulator_.getMHDModel(), lvl};
    }

    auto getHybridPatchLevel(size_t lvl)
        requires(has_hybrid_v<opts>)
    {
        using HybridModel = solver::PHARE_Types<opts>::Hybrid::Model_t;
        return PatchLevel<HybridModel>{*hierarchy_, *simulator_.getHybridModel(), lvl};
    }


    auto sync(std::vector<PatchData<py_array_t<double>, dimension>> const& input)
    {
        // collect all data on rank 0!

        auto const mpi_size = mpi::size();
        std::vector<PatchData<py_array_t<double>, dimension>> collected;

        auto const collect = [&](PatchData<py_array_t<double>, dimension> const& patch_data) {
            auto patchIDs = mpi::collect(patch_data.patchID, mpi_size);
            auto shapes   = mpi::collectArrays(shape<dimension>(patch_data.data), mpi_size);
            auto origins  = mpi::collect(makeSpan(patch_data.origin), mpi_size);
            auto lower    = mpi::collect(makeSpan(patch_data.lower), mpi_size);
            auto upper    = mpi::collect(makeSpan(patch_data.upper), mpi_size);
            auto ghosts   = mpi::collect(patch_data.nGhosts, mpi_size);
            auto datas    = mpi::collect(makeSpan(patch_data.data), mpi_size);

            if (mpi::rank() == 0)
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

        auto const max = mpi::max(input.size());

        PatchData<py_array_t<double>, dimension> empty{core::ConstArray<int, dimension>()};

        for (std::size_t i = 0; i < max; ++i)
            collect(i < input.size() ? input[i] : empty);

        return collected;
    }


private:
    Simulator_t& simulator_;
    std::shared_ptr<amr::Hierarchy> hierarchy_;
};

} // namespace PHARE::pydata

#endif /*PHARE_PYTHON_DATA_WRANGLER_H*/
